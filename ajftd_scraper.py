"""
AJFTD scraper — extracts all accessible fuel records and selected properties.

Properties:
  - density (including density_api_by_meter as a distinct metric)
  - viscosity (including viscosity_kinematic and related metrics)
  - flash_point
  - smoke_point

Output:
  ajftd_data.csv (one row per reported measurement)
  debug/initial_results.html
  debug/detail_<fuel_id>.html

Install:
  pip install playwright beautifulsoup4
  python -m playwright install chromium

Run:
  python ajftd_scraper.py

Sign in manually in Chromium, navigate to Fuels > Search Fuels, clear filters,
then return to the terminal and press Enter.

Changes in this version (everything else is unchanged):
  * RESUME: a detail page already saved in debug/detail_<id>.html is re-used instead of re-fetching it, so a crash or a
    slow-site timeout no longer loses the whole run - just start it again.  Use `--fresh` to force re-fetching.
  * RETRIES: each detail page is attempted 3 times before being recorded as an error.
  * FUEL NAME: taken from the "Fuel Details for (...)" heading (before, every record got the site banner
    "Alternative Jet Fuels"); also written to a new fuel_label column.
  * The results table's "Showing x to y of N entries" text is printed, and you get a warning if the number of
    records discovered differs from N (a sign that pagination stopped early).
  * To get EVERY metric on the page (not only the ones in PROPERTY_GROUPS) run:  python ajftd_parse_html.py
    It re-reads the saved debug/*.html offline.
"""

from __future__ import annotations

import json
import csv
import re
import sys
import time
from pathlib import Path
from urllib.parse import urljoin, urlparse, parse_qs

from bs4 import BeautifulSoup
from playwright.sync_api import sync_playwright, TimeoutError as PlaywrightTimeout

BASE_URL = "https://altjetfuels.illinois.edu/"
SEARCH_URL = "https://altjetfuels.illinois.edu/fuels/"
SESSION_FILE = Path("ajftd_session.json")
OUTPUT_FILE = Path("ajftd_data.csv")
DEBUG_DIR = Path("debug")

# Be respectful of the service. This script visits pages sequentially.
DELAY_BETWEEN_RECORDS_SECONDS = 0.35
PAGE_LOAD_TIMEOUT_MS = 20000
MAX_RESULT_PAGES = 10000
MAX_RECORDS = 100000
RESUME = "--fresh" not in sys.argv
DETAIL_ATTEMPTS = 3
# Discover records through the site's AJAX endpoint (fast, no clicking). `--no-ajax` = old click-through paging.
USE_AJAX_DISCOVERY = "--no-ajax" not in sys.argv

# Exact metric labels seen in AJFTD detail HTML, plus common alternatives.
# Keep metrics distinct: density and API gravity are related but not identical.
# Requested properties plus composition metrics visible in AJFTD detail pages.
# The scraper retains the original metric label and reported unit.
PROPERTY_GROUPS = {
    "density": {"density", "fuel_density"},
    "density_api_by_meter": {"density_api_by_meter", "api_gravity", "api"},
    "viscosity": {"viscosity", "viscosity_dynamic", "dynamic_viscosity"},
    "kinematic_viscosity": {"viscosity_kinematic", "kinematic_viscosity"},
    "flash_point": {"flash_point", "flashpoint", "flash_temperature"},
    "pour_point": {"pour_point", "pourpoint"},
    "smoke_point": {"smoke_point", "smokepoint"},
    # Composition metrics
    "composition": {
        "acid_number", "acidity_total", "aromatics_total",
        "contaminants_existent_gum",
        "contaminants_water_reaction_interface_condition",
        "contaminants_water_reaction_separation",
        "diaromatics", "monoaromatics", "olefins_total",
        "saturates_total", "sulfur", "sulfur_mercaptan",
    },
    # Volume-related metrics. Also captures volume fields in the measurement table.
    "volume": {
        "contaminants_water_reaction_volume_change",
        "thermal_stability_spent_fuel_volume",
        "distillation", "distillation_loss", "distillation_recovery",
        "distillation_residue",
    },
}

def clean(value: str | None) -> str:
    return re.sub(r"\s+", " ", value or "").strip()

def normalize_metric(value: str) -> str:
    value = clean(value).lower()
    value = re.sub(r"\([^)]*\)", " ", value)
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")

def get_property_group(metric_label: str) -> str | None:
    normalized = normalize_metric(metric_label)
    for group, aliases in PROPERTY_GROUPS.items():
        if normalized in aliases:
            return group
    # Preserve any additional composition metric labels instead of discarding them.
    known_composition_prefixes = (
        "aromatics", "diaromatics", "monoaromatics", "olefins",
        "saturates", "sulfur", "acid", "acidity",
    )
    if normalized.startswith(known_composition_prefixes) or (
        normalized.startswith("contaminants_")
        and normalized != "contaminants_water_reaction_volume_change"
    ):
        return "composition"
    # Keep additional volume-related metrics in the output.
    if "volume" in normalized:
        return "volume"
    return None

def page_text(page) -> str:
    try:
        return page.locator("body").inner_text(timeout=5000)
    except Exception:
        return ""

def save_html(html: str, filename: str) -> None:
    DEBUG_DIR.mkdir(parents=True, exist_ok=True)
    (DEBUG_DIR / filename).write_text(html, encoding="utf-8", errors="replace")

def extract_metric_table(table, metric_label: str) -> list[dict]:
    """Extract rows from the table paired with a metric's expand/collapse button."""
    rows = table.select("tr")
    if not rows:
        return []

    headers = [
        clean(cell.get_text(" ", strip=True))
        for cell in rows[0].select("th, td")
    ]
    measurements = []

    for tr in rows[1:]:
        cells = [
            clean(cell.get_text(" ", strip=True))
            for cell in tr.select("th, td")
        ]
        if not cells:
            continue

        # Skip rows that are entirely empty.
        if not any(cells):
            continue

        row_obj = {}
        for idx, value in enumerate(cells):
            key = headers[idx] if idx < len(headers) and headers[idx] else f"column_{idx + 1}"
            # Repeated headers are unlikely, but avoid overwriting values.
            if key in row_obj:
                key = f"{key}_{idx + 1}"
            row_obj[key] = value if value else None

        # Preserve original source fields; do not coerce values or convert units.
        measurements.append(row_obj)

    return measurements

def extract_fuel_record(html: str, record_url: str) -> dict:
    soup = BeautifulSoup(html, "html.parser")
    title = soup.title.get_text(" ", strip=True) if soup.title else None

    # The record name is generally shown in the page heading or near the fuel ID.
    headings = [clean(h.get_text(" ", strip=True)) for h in soup.select("h1, h2, h3")]
    headings = [h for h in headings if h]

    parsed_url = urlparse(record_url)
    query = parse_qs(parsed_url.query)
    fuel_id = query.get("id", [None])[0]

    # The first heading is the site banner; the fuel's own name is in "Fuel Details for (...)".
    detail_heading = next((h for h in headings if h.lower().startswith("fuel details for")), None)
    fuel_label = re.sub(r"(?i)^fuel details for\s*", "", detail_heading).strip() if detail_heading else None

    # Detail page uses buttons like onclick="myFunction('tbl_15')" whose visible
    # label is the metric identifier, and tables with that exact ID.
    properties = {group: [] for group in PROPERTY_GROUPS}
    all_metrics_found = []

    for button in soup.select("button[onclick*='myFunction']"):
        onclick = button.get("onclick", "")
        match = re.search(r"""myFunction\(\s*['"]([^'"]+)['"]\s*\)""", onclick)
        if not match:
            continue

        table_id = match.group(1)
        metric_label = clean(button.get_text(" ", strip=True))
        metric_label = re.sub(r"\s+", " ", metric_label)
        # Buttons may include a trailing whitespace/icon; normalize only for matching.
        group = get_property_group(metric_label)
        if not group:
            continue

        table = soup.find(id=table_id)
        if not table:
            continue

        measurements = extract_metric_table(table, metric_label)
        if measurements:
            properties[group].append({
                "metric": metric_label,
                "table_id": table_id,
                "measurements": measurements,
            })
            all_metrics_found.append(metric_label)

    # Keep the record even when no requested property has data; missing data is
    # meaningful and must not be replaced with guessed values.
    return {
        "fuel_id": fuel_id,
        "name": fuel_label or (headings[0] if headings else None),
        "fuel_label": fuel_label,
        "headings": headings,
        "record_url": record_url,
        "page_title": title,
        "properties": properties,
        "metrics_extracted": sorted(set(all_metrics_found)),
    }

def cached_detail_html(fuel_id: str) -> str | None:
    """Return a previously saved detail page (resume support) or None."""
    if not RESUME or not fuel_id:
        return None
    path = DEBUG_DIR / f"detail_{fuel_id}.html"
    try:
        if path.exists() and path.stat().st_size > 2000:
            html = path.read_text(encoding="utf-8", errors="replace")
            if "Fuel Details for" in html:       # a complete page, not an error/login page
                return html
    except OSError:
        pass
    return None


def results_signature(page) -> str:
    """Fingerprint of the rows currently shown in the results table.

    NOTE: do not fingerprint body text - its first ~11,000 characters are the banner, menu and the filter
    sidebar (metric lists), which are identical on every results page, so a "first 2000 chars" check
    can never see the table change and reports "page did not advance" even when it did.
    """
    try:
        hrefs = page.locator('#tbl_fuels a[title="click to view this fuel"]').evaluate_all(
            "els => els.map(e => e.getAttribute('href'))")
        return "|".join(h for h in hrefs if h)
    except Exception:
        return ""


def get_result_links(page, allowed_host: str) -> list[str]:
    links = []
    # Verified from the supplied AJFTD Search Fuels HTML.
    locator = page.locator('a[title="click to view this fuel"]')
    try:
        count = locator.count()
    except Exception:
        count = 0

    for idx in range(count):
        try:
            href = locator.nth(idx).get_attribute("href")
        except Exception:
            continue
        if not href:
            continue
        absolute = urljoin(page.url, href)
        parsed = urlparse(absolute)
        if parsed.netloc != allowed_host:
            continue
        if "fuel-viewcompare.asp" not in parsed.path.lower():
            continue
        if parse_qs(parsed.query).get("id"):
            links.append(absolute)

    # Deduplicate while preserving order.
    return list(dict.fromkeys(links))

def next_page_button(page):
    # AJFTD uses a DataTables control: numbered <a class="page-link" data-dt-idx="N"> links plus a
    # #pageInput box. There is no "Next" link, so walk the numbered links from the current page.
    cur_link = page.locator('a[aria-controls="tbl_fuels"][aria-current="page"]').first
    cur = None
    try:
        if cur_link.count():
            cur = int(cur_link.get_attribute("data-dt-idx") or "")
    except Exception:
        cur = None
    if cur is None:
        try:
            cur = int(page.locator("#pageInput").input_value() or "1")
        except Exception:
            return None
    nxt = page.locator(f'a.page-link[aria-controls="tbl_fuels"][data-dt-idx="{cur + 1}"]').first
    try:
        if nxt.count() and nxt.is_visible():
            return nxt
    except Exception:
        pass
    return None

def main():
    DEBUG_DIR.mkdir(parents=True, exist_ok=True)
    all_record_urls = []
    seen_result_pages = set()
    visited_page_signatures = set()

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=False)

        if SESSION_FILE.exists():
            try:
                context = browser.new_context(storage_state=str(SESSION_FILE))
            except Exception:
                context = browser.new_context()
        else:
            context = browser.new_context()

        page = context.new_page()
        page.set_default_timeout(PAGE_LOAD_TIMEOUT_MS)

        page.goto(BASE_URL, wait_until="domcontentloaded")
        print("\nA Chromium window is open.")
        print("1. Sign in if prompted.")
        print("2. Open Fuels > Search Fuels.")
        print("3. Click 'Clear All Filters' so the full accessible list is shown.")
        print("4. If possible, choose a larger 'entries per page' value.")
        try:
            input("\nWhen the Search Fuels results table is visible, press Enter here... ")
        except EOFError:
            print("No interactive input available. Waiting 30 seconds for manual interaction...")
            time.sleep(30)

        # Save login state for the next run. This file is sensitive.
        try:
            context.storage_state(path=str(SESSION_FILE))
        except Exception as exc:
            print(f"Could not save browser session state: {exc}")

        start_url = page.url
        allowed_host = urlparse(start_url).netloc
        if not allowed_host:
            allowed_host = urlparse(BASE_URL).netloc

        # If user is on a detail or unrelated page, navigate to the search page.
        if "fuel-viewcompare.asp" in urlparse(start_url).path.lower():
            page.goto(SEARCH_URL, wait_until="domcontentloaded")
        elif "fuels" not in urlparse(start_url).path.lower():
            page.goto(SEARCH_URL, wait_until="domcontentloaded")

        try:
            page.locator('a[title="click to view this fuel"]').first.wait_for(
                state="visible", timeout=PAGE_LOAD_TIMEOUT_MS
            )
        except PlaywrightTimeout:
            print("No fuel links found yet. Confirm you are on the Search Fuels results page.")
            save_html(page.content(), "initial_results.html")
            print("Saved debug/initial_results.html for inspection.")
            context.close()
            browser.close()
            return

        save_html(page.content(), "initial_results.html")
        try:
            # searchDatatable is a top-level `let` (global lexical binding, NOT a window
            # property) - wait on the bare identifier rather than window.searchDatatable.
            page.wait_for_function(
                "() => typeof searchDatatable !== 'undefined' && searchDatatable !== null",
                timeout=20000)
            page.evaluate("() => { searchDatatable.page.len(100); searchDatatable.draw(false); }")
            page.wait_for_timeout(3000)
        except Exception as exc:
            print(f"[diag] widen via searchDatatable failed: {exc}")
        save_html(page.content(), "initial_results.html")
        try:
            info_now = page.locator("#tbl_fuels_info, .dt-info, .dataTables_info").first.inner_text(timeout=3000)
            print(f"[diag] table now says: {info_now}")
        except Exception as exc:
            print(f"[diag] could not read info text: {exc}")
        print(f"Search page: {page.url}")
        expected_total = None
        try:
            info_text = clean(page.locator("#tbl_fuels_info, .dt-info, .dataTables_info").first.inner_text(timeout=3000))
            print(f"Results table says: {info_text}")
            m_total = re.search(r"of\s+([\d,]+)\s+entries", info_text)
            expected_total = int(m_total.group(1).replace(",", "")) if m_total else None
        except Exception:
            print("(no 'Showing x to y of N entries' text found - cannot cross-check the record count)")

        # --- Discover every record through the site's own DataTables endpoint -----------------------
        # Runs inside the logged-in page with the site's jQuery and `searchDatatable.ajax.params()`, so the
        # request parameters and cookies are exactly what the site itself sends (see ajax_paginator.py).
        if USE_AJAX_DISCOVERY:
            try:
                from ajax_paginator import browser_fetcher, crawl_all, write_index
                records, ajax_total = crawl_all(browser_fetcher(page), length=100, delay=0.2)
                all_record_urls = [r["record_url"] for r in records]
                write_index(records, Path("data/processed/ajftd_fuel_index.csv"))
                expected_total = ajax_total or expected_total
                print(f"[ajax] discovered {len(all_record_urls)} fuels (server total: {ajax_total})")
            except Exception as exc:
                print(f"[ajax] discovery failed ({exc}); falling back to click-through pagination")
                all_record_urls = []

        # --- Fallback: click through the table page by page --------------------------------------------
        if not all_record_urls:
            # Crawl visible pages using the actual Next button. Stop if page contents
            # repeat or no enabled Next button exists.
            for page_num in range(1, MAX_RESULT_PAGES + 1):
                current_url = page.url
                html = page.content()
                signature = results_signature(page)

                if signature in visited_page_signatures:
                    print("Repeated results-page content detected; stopping pagination.")
                    break
                visited_page_signatures.add(signature)
                seen_result_pages.add(current_url)

                page_links = get_result_links(page, allowed_host)
                for link in page_links:
                    if link not in all_record_urls:
                        all_record_urls.append(link)

                print(
                    f"Results page {page_num}: found {len(page_links)} fuel links "
                    f"(unique total: {len(all_record_urls)})"
                )

                if len(all_record_urls) >= MAX_RECORDS:
                    print("Safety limit reached.")
                    break

                before_signature = signature
                try:
                    # Advance the way the site's own enterPageNumber does: page(N).draw(false) on the
                    # global searchDatatable (zero-based; #pageInput value is 1-based so it equals the
                    # next zero-based index). Keyboard/dispatch/link approaches didn't fire the handler.
                    zero_based_next = int(page.locator("#pageInput").input_value() or "1")
                    page.wait_for_function(
                        "() => typeof searchDatatable !== 'undefined' && searchDatatable !== null",
                        timeout=20000)
                    page.evaluate(f"() => {{ searchDatatable.page({zero_based_next}).draw(false); }}")
                    # Wait for the body text to change; the server-side redraw can take seconds.
                    deadline = time.time() + 120
                    moved = False
                    while time.time() < deadline:
                        now_sig = results_signature(page)
                        if now_sig and now_sig != before_signature:
                            moved = True
                            break
                        page.wait_for_timeout(1000)
                    if not moved:
                        print(f"Page did not advance after 120 s (tried zero-based page {zero_based_next}); falling back to Next link click.")
                        # Fallback: click the "Next" pagination link (always advances by 1 page)
                        next_link = page.locator('a.page-link[data-dt-idx="next"][aria-controls="tbl_fuels"]').first
                        if next_link.count() > 0 and next_link.is_visible():
                            next_link.click()
                            deadline2 = time.time() + 60
                            while time.time() < deadline2:
                                now_sig2 = results_signature(page)
                                if now_sig2 and now_sig2 != before_signature:
                                    moved = True
                                    break
                                page.wait_for_timeout(500)
                        if not moved:
                            print(f"Could not advance page via Next link; stopping pagination.")
                            break
                except Exception as exc:
                    print(f"Could not advance page: {exc}")
                    break

        print(f"\nUnique fuel records discovered: {len(all_record_urls)}")
        if expected_total and expected_total != len(all_record_urls):
            print(f"WARNING: the table reports {expected_total} entries but only {len(all_record_urls)} links were collected - "
                  "pagination may have stopped early (or filters are still active).")
        if not all_record_urls:
            print("No records discovered. Inspect debug/initial_results.html.")
            context.close()
            browser.close()
            return

        # Detail pages use the exact metric button/table relationship visible in
        # the supplied detail HTML.
        fuel_records = []
        for i, record_url in enumerate(all_record_urls, start=1):
            try:
                query = parse_qs(urlparse(record_url).query)
                fuel_id = query.get("id", [str(i)])[0]
                html = cached_detail_html(fuel_id)
                from_cache = html is not None
                if not from_cache:
                    last_exc = None
                    for attempt in range(1, DETAIL_ATTEMPTS + 1):
                        try:
                            page.goto(record_url, wait_until="domcontentloaded")
                            try:
                                page.wait_for_load_state("networkidle", timeout=8000)
                            except PlaywrightTimeout:
                                pass
                            page.wait_for_timeout(350)
                            html = page.content()
                            last_exc = None
                            break
                        except Exception as exc:
                            last_exc = exc
                            print(f"    attempt {attempt}/{DETAIL_ATTEMPTS} failed for id={fuel_id}: {exc}")
                            time.sleep(2 * attempt)
                    if last_exc is not None:
                        raise last_exc
                    save_html(html, f"detail_{fuel_id}.html")      # saved immediately: progress survives a crash

                record = extract_fuel_record(html, record_url)
                fuel_records.append(record)

                found = [
                    key for key, groups in record["properties"].items()
                    if groups
                ]
                print(
                    f"[{i}/{len(all_record_urls)}] "
                    f"id={fuel_id} name={record.get('name')!r} "
                    f"properties={','.join(found) if found else 'none'}"
                    f"{' (cached)' if from_cache else ''}"
                )
                if not from_cache:
                    time.sleep(DELAY_BETWEEN_RECORDS_SECONDS)

            except Exception as exc:
                print(f"[{i}/{len(all_record_urls)}] ERROR {record_url}: {exc}")
                fuel_records.append({
                    "fuel_id": parse_qs(urlparse(record_url).query).get("id", [None])[0],
                    "name": None,
                    "record_url": record_url,
                    "error": str(exc),
                    "properties": {
                        **{group: [] for group in PROPERTY_GROUPS},
                    },
                })

        # Flatten nested fuel/property/measurement data to CSV: one row per measurement.
        # This preserves multiple measurements for a fuel/property instead of losing data.
        csv_rows = []
        all_measurement_fields = set()

        for record in fuel_records:
            base = {
                "fuel_id": record.get("fuel_id"),
                "fuel_name": record.get("name"),
                "fuel_label": record.get("fuel_label"),
                "record_url": record.get("record_url"),
                "page_title": record.get("page_title"),
            }

            for group, metric_groups in (record.get("properties") or {}).items():
                for metric_group in metric_groups or []:
                    metric_label = metric_group.get("metric")
                    table_id = metric_group.get("table_id")
                    measurements = metric_group.get("measurements") or []

                    # Keep a row even if a recognized metric table exists but has no data rows.
                    if not measurements:
                        csv_rows.append({
                            **base,
                            "category": group,
                            "metric": metric_label,
                            "table_id": table_id,
                            "measurement_index": None,
                            "Value": None,
                            "Unit": None,
                        })
                        continue

                    for measurement_index, measurement in enumerate(measurements, start=1):
                        row = {
                            **base,
                            "category": group,
                            "metric": metric_label,
                            "table_id": table_id,
                            "measurement_index": measurement_index,
                        }
                        row.update(measurement or {})
                        csv_rows.append(row)
                        all_measurement_fields.update((measurement or {}).keys())

        preferred_columns = [
            "fuel_id", "fuel_name", "fuel_label", "record_url", "page_title",
            "category", "metric", "table_id", "measurement_index",
            "Value", "Unit", "temperature_value", "temperature_unit",
            "volume_evaporated_value", "volume_evaporated_unit",
            "specification_limits", "info", "Test Date", "Test Method",
        ]
        fieldnames = [c for c in preferred_columns if c in (
            set(preferred_columns) | all_measurement_fields |
            {"fuel_id", "fuel_name", "fuel_label", "record_url", "page_title", "category", "metric", "table_id", "measurement_index"}
        )]
        fieldnames += sorted(all_measurement_fields - set(fieldnames))

        with OUTPUT_FILE.open("w", newline="", encoding="utf-8-sig") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(csv_rows)

        counts = {}
        for key in PROPERTY_GROUPS:
            counts[key] = sum(
                1 for record in fuel_records
                if record.get("properties", {}).get(key)
            )

        print("\nFinished.")
        print(f"CSV: {OUTPUT_FILE.resolve()}")
        print(f"Debug HTML: {DEBUG_DIR.resolve()}")
        print(f"Records discovered / processed: {len(all_record_urls)} / {len(fuel_records)}")
        print(f"Measurement rows exported: {len(csv_rows)}")
        print("Records with each category:")
        for key, count in counts.items():
            print(f"  {key}: {count}")

        context.close()
        browser.close()

if __name__ == "__main__":
    main()
