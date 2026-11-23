"""Interactive AJFTD scraper for targeted / filtered fuels (e.g. 178 GCxGC fuels).

Launches a visible Chromium browser session:
1. Navigates to https://altjetfuels.illinois.edu/fuels/
2. Allows user to sign in if required.
3. User applies their desired filter/search on the page (e.g. the 178 fuels).
4. An on-screen floating banner with a "Start Scraping" button lets the user trigger scraping directly in the browser.
5. Queries DataTables AJAX endpoint with the current active filter to get all matching fuel IDs.
6. Downloads detail HTML pages to extras/fuelblend_data/debug/detail_<id>.html (skipping existing ones).
7. Re-runs parsing to update ajftd_fuel_details.csv and loaders_ajftd_gcxgc.py.
"""
from __future__ import annotations

import csv
import json
import re
import sys
import time
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from bs4 import BeautifulSoup
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
DEBUG_DIR = ROOT / "extras" / "fuelblend_data" / "debug"
RAW_DIR = ROOT / "extras" / "fuelblend_data" / "data" / "raw" / "ajftd"
SESSION_FILE = ROOT / "extras" / "fuelblend_data" / "ajftd_session.json"
BASE_URL = "https://altjetfuels.illinois.edu/"
SEARCH_URL = "https://altjetfuels.illinois.edu/fuels/"
DETAIL_URL = "https://altjetfuels.illinois.edu/fuels/fuel-viewcompare.asp"

_BANNER_JS = """
(() => {
  if (document.getElementById('scraper-overlay')) return;
  const div = document.createElement('div');
  div.id = 'scraper-overlay';
  div.style.cssText = 'position:fixed;top:15px;right:20px;z-index:9999999;background:#0f172a;color:#f8fafc;padding:16px 20px;border-radius:10px;font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;box-shadow:0 10px 25px rgba(0,0,0,0.5);border:2px solid #38bdf8;max-width:380px;';
  div.innerHTML = `
    <div style="font-weight:700;font-size:15px;margin-bottom:6px;color:#38bdf8;">AJFTD Target Scraper</div>
    <div id="scraper-status-text" style="font-size:13px;line-height:1.4;margin-bottom:12px;color:#cbd5e1;">
      1. Log in if prompted.<br/>
      2. Filter/search until your target fuels (e.g. 178) appear in the table.<br/>
      3. Click the button below to scrape.
    </div>
    <div style="display:flex;gap:10px;align-items:center;">
      <button id="scraper-btn-start" style="background:#0284c7;color:#fff;border:none;padding:8px 16px;border-radius:6px;font-weight:600;font-size:13px;cursor:pointer;transition:background 0.2s;">
        Start Scraping Filtered Fuels
      </button>
      <span id="scraper-count-badge" style="font-size:12px;color:#94a3b8;"></span>
    </div>
  `;
  document.body.appendChild(div);

  const btn = document.getElementById('scraper-btn-start');
  btn.addEventListener('mouseenter', () => btn.style.background = '#0369a1');
  btn.addEventListener('mouseleave', () => btn.style.background = '#0284c7');
  btn.addEventListener('click', () => {
    btn.disabled = true;
    btn.style.background = '#64748b';
    btn.textContent = 'Scraping in progress...';
    window.__scraper_requested = true;
  });

  // Periodically update fuel count badge from table info text
  setInterval(() => {
    const info = document.querySelector('#tbl_fuels_info, .dataTables_info');
    if (info) {
      const badge = document.getElementById('scraper-count-badge');
      if (badge) {
        const m = info.innerText.match(/of\\s+([\\d,]+)\\s+entries/);
        if (m) badge.textContent = `${m[1]} fuels matched`;
      }
    }
  }, 1000);
})();
"""


def get_filtered_records(page) -> list[dict]:
    """Retrieve all records matching the current DataTables filter using page's jQuery, paging until done."""
    fetch_page_js = """
    async ({start, length, draw}) => {
      if (typeof searchDatatable === 'undefined' || !searchDatatable) {
        return { error: 'searchDatatable not found' };
      }
      const params = Object.assign({}, searchDatatable.ajax.params());
      params.start = start;
      params.length = length;
      params.draw = draw;
      return await new Promise((resolve, reject) => {
        jQuery.ajax({
          url: searchDatatable.ajax.url(),
          method: 'POST',
          data: params,
          dataType: 'json',
          success: resolve,
          error: (xhr, status, err) => reject(new Error('HTTP ' + xhr.status + ' ' + err))
        });
      });
    }
    """
    records = []
    start = 0
    length = 100
    draw = 1
    total_expected = None

    while True:
        res = page.evaluate(fetch_page_js, {"start": start, "length": length, "draw": draw})
        if isinstance(res, dict) and "error" in res:
            raise RuntimeError(res["error"])
        
        if total_expected is None:
            total_expected = int(res.get("recordsFiltered") or res.get("recordsTotal") or 0)
            print(f"[Scraper] Server reports recordsFiltered={total_expected}")

        data = res.get("data", [])
        if not data:
            break

        for item in data:
            if isinstance(item, dict):
                fid = item.get("FuelID") or item.get("fuel_id") or item.get("ID")
                fname = item.get("FuelName") or item.get("fuel_name") or ""
            elif isinstance(item, list) and len(item) > 0:
                first_col = str(item[0])
                m = re.search(r"id=(\d+)", first_col)
                fid = m.group(1) if m else None
                soup = BeautifulSoup(first_col, "html.parser")
                fname = soup.get_text(strip=True)
            else:
                continue
            if fid:
                records.append({
                    "fuel_id": str(fid),
                    "fuel_name": fname,
                    "record_url": f"{DETAIL_URL}?id={fid}"
                })

        start += len(data)
        draw += 1
        if len(data) < length or (total_expected and start >= total_expected):
            break
        time.sleep(0.2)

    # De-duplicate by fuel_id while preserving order
    seen = set()
    dedup = []
    for r in records:
        if r["fuel_id"] not in seen:
            seen.add(r["fuel_id"])
            dedup.append(r)
    return dedup


def update_overlay_status(page, msg: str):
    page.evaluate(f"""
    (() => {{
      const el = document.getElementById('scraper-status-text');
      if (el) el.innerHTML = `{msg}`;
    }})()
    """)


def scrape_filtered_fuels():
    DEBUG_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 70)
    print("AJFTD TARGETED FUEL SCRAPER")
    print("=" * 70)
    print("Launching Chromium browser...")

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=False)
        context = browser.new_context(
            storage_state=str(SESSION_FILE) if SESSION_FILE.exists() else None
        )
        page = context.new_page()
        page.goto(SEARCH_URL, wait_until="domcontentloaded")

        print("\nBrowser is open at https://altjetfuels.illinois.edu/fuels/")
        print("Please log in if prompted, and set your search/filter (to show the 178 fuels).")
        print("When ready, click the blue 'Start Scraping Filtered Fuels' button on the page.")

        # Keep injecting the overlay and waiting for click
        while True:
            try:
                page.evaluate(_BANNER_JS)
                is_clicked = page.evaluate("() => Boolean(window.__scraper_requested)")
                if is_clicked:
                    break
            except Exception:
                pass
            time.sleep(1)

        print("\n[Scraper] 'Start Scraping' clicked! Reading active filter records...")
        update_overlay_status(page, "Reading records from filtered table...")

        records = get_filtered_records(page)
        print(f"[Scraper] Discovered {len(records)} fuels matching current filter.")
        update_overlay_status(page, f"Found <b>{len(records)}</b> matching fuels.<br/>Fetching detail pages...")

        # Save session cookies for reuse
        try:
            context.storage_state(path=str(SESSION_FILE))
        except Exception:
            pass

        # Scrape detail HTML for each fuel
        fetched_count = 0
        cached_count = 0
        for i, rec in enumerate(records, start=1):
            fid = rec["fuel_id"]
            url = rec["record_url"]
            out_file = DEBUG_DIR / f"detail_{fid}.html"

            if out_file.exists() and out_file.stat().st_size > 2000:
                content = out_file.read_text(encoding="utf-8", errors="ignore")
                if "Fuel Details for" in content:
                    cached_count += 1
                    if i % 10 == 0 or i == len(records):
                        print(f"[{i}/{len(records)}] Fuel {fid} (cached)")
                    continue

            # Fetch fresh detail page
            update_overlay_status(page, f"Fetching fuel <b>{i} / {len(records)}</b> (ID: {fid})...")
            try:
                # Use page to fetch HTML with current session
                resp = page.request.get(url)
                if resp.status == 200:
                    html = resp.text()
                    out_file.write_text(html, encoding="utf-8", errors="replace")
                    fetched_count += 1
                    print(f"[{i}/{len(records)}] Fuel {fid} -> saved ({len(html)} bytes)")
                else:
                    print(f"[{i}/{len(records)}] Fuel {fid} -> HTTP {resp.status}")
            except Exception as e:
                print(f"[{i}/{len(records)}] Fuel {fid} -> ERROR: {e}")

            time.sleep(0.3)

        print(f"\n[Scraper] Done: {fetched_count} fetched, {cached_count} cached.")
        update_overlay_status(page, f"<b>Scraping Complete!</b><br/>{fetched_count} newly fetched, {cached_count} from cache.<br/>Processing data...")

        browser.close()

    print("\nRunning offline HTML parser on debug/...")
    from subprocess import run
    parse_script = ROOT / "extras" / "fuelblend_data" / "ajftd_parse_html.py"
    if parse_script.exists():
        run([sys.executable, str(parse_script), str(DEBUG_DIR)], check=False)

    print("\nRunning v5 GCxGC composition loader...")
    from fuelblend.data.loaders_ajftd_gcxgc import main as load_gcxgc
    load_gcxgc()

    print("\n" + "=" * 70)
    print("ALL DONE! Check data/real_fuels/v5/ajftd_gcxgc_fuels.csv")
    print("=" * 70)


if __name__ == "__main__":
    scrape_filtered_fuels()
