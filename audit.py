#!/usr/bin/env python
"""Check what worked and what is missing, and produce a small report you can send back.

    python audit.py            # prints a summary, writes AUDIT.md
    python audit.py --bundle   # also writes audit_bundle.zip (AUDIT.md + small diagnostic files, no raw data)

Status codes:  OK | WARN (works but look) | MISSING (expected, not there) | FAIL (there but wrong)
Run from the fuelblend_data folder, after fetch.py / inspect_downloads.py / build_* scripts.
"""
from __future__ import annotations
import hashlib, importlib, json, platform, re, sys, zipfile
from pathlib import Path

RAW, INT, PROC = Path("data/raw"), Path("data/interim"), Path("data/processed")
checks: list[tuple[str, str, str, str]] = []   # status, area, item, detail/fix


def add(status, area, item, detail=""):
    checks.append((status, area, item, detail))


def junk(p: Path) -> bool:
    return "__MACOSX" in p.parts or p.name.startswith("._") or p.name.startswith("~$") or p.name == ".DS_Store"


def find(root: Path, pattern: str) -> Path | None:
    if not root.exists():
        return None
    for p in sorted(root.rglob(pattern)):
        if p.is_file() and not junk(p):
            return p
    return None


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""):
            h.update(c)
    return h.hexdigest()


def nrows(p: Path) -> int:
    n = 0
    with p.open("rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""):
            n += c.count(b"\n")
    return max(n - 1, 0)


# ------------------------------------------------------------------ A. environment
def env():
    add("OK", "env", f"python {platform.python_version()} on {platform.system()}")
    for m in ("pandas", "numpy", "rdkit", "requests", "lxml", "openpyxl"):
        try:
            mod = importlib.import_module(m)
            add("OK", "env", m, getattr(mod, "__version__", ""))
        except Exception as e:  # noqa: BLE001
            add("MISSING" if m != "openpyxl" else "WARN", "env", m, f"pip install {m if m != 'rdkit' else 'rdkit'}  ({e})")
    if len(str(Path.cwd())) > 120:
        add("WARN", "env", "very long folder path", "Windows 260-char limit / OneDrive sync can break deep zip extraction; move project to e.g. C:\\work\\")
    if "onedrive" in str(Path.cwd()).lower():
        add("WARN", "env", "project is inside OneDrive", "sync of ~800 small files + 400 MB slows everything; move to a non-synced folder")


# ------------------------------------------------------------------ B. raw downloads
EXPECT_PDFS = ["NISTIR_6646_RP1.pdf", "CRC_647_world_fuel_sampling.pdf", "CRC_AV-27-18_fuel_properties.pdf",
               "CRC_AV-11-09_freezing_point.pdf", "AFRL_TR-2020-0017_jet_fuel_properties.pdf",
               "CRC_handbook_aviation_fuel_properties.pdf", "Zabetakis_Bulletin_627_flammability.pdf",
               "Li_2020_flash_point_GNN_arXiv2002.11315.pdf", "Leenhouts_2026_mixture_eval_arXiv2605.29698.pdf",
               "Chew_2024_viscosity_SI.pdf"]
TGZ_SHA = "231161b5e443dc1ae0e5da8429d86a88474cb722016e5b790817bb31c58d7ec2"


def raw():
    for g in ("chemixhub", "mixture_evaluation"):
        d = RAW / "git" / g
        if (d / ".git").exists():
            add("OK", "raw", f"git:{g}")
        else:
            add("MISSING", "raw", f"git:{g}", "python fetch.py git   (mixture_evaluation is on gitlab.kuleuven.be; optional)")
    nist = find(RAW / "git" / "chemixhub", "processed_NISTlogV.csv")
    add("OK" if nist else "MISSING", "raw", "CheMixHub nist-logV table", str(nist or ""))
    for label, fn, size in (("leenhouts_benchmark", "datasets.zip", 100e6), ("leenhouts_benchmark", "splits.zip", 100e6),
                            ("bilodeau_nist_viscosity", "pretrained_models.zip", 50e6)):
        p = RAW / "zenodo" / label / fn
        if not p.exists():
            add("MISSING", "raw", f"zenodo:{label}/{fn}", "python fetch.py zenodo  (or browser download)")
        elif p.stat().st_size < size:
            add("FAIL", "raw", f"zenodo:{label}/{fn}", f"only {p.stat().st_size/1e6:.1f} MB - truncated? delete + refetch")
        else:
            add("OK", "raw", f"zenodo:{label}/{fn}", f"{p.stat().st_size/1e6:.0f} MB")
    optional = {"AFRL_TR-2020-0017_jet_fuel_properties.pdf", "CRC_handbook_aviation_fuel_properties.pdf"}   # DTIC; AJFTD covers most of it
    for fn in EXPECT_PDFS:
        p = RAW / "pdfs" / fn
        bad = "WARN" if fn in optional else "FAIL"
        miss = "WARN" if fn in optional else "MISSING"
        if not p.exists():
            add(miss, "raw", f"pdf:{fn}", ("optional - DTIC site; get later in a browser. " if fn in optional else "") +
                "save into data/raw/pdfs/ with this exact name")
        else:
            ok = p.read_bytes()[:5] == b"%PDF-"
            add("OK" if ok and p.stat().st_size > 20_000 else bad, "raw", f"pdf:{fn}",
                f"{p.stat().st_size/1e6:.1f} MB" + ("" if ok else " - not a real PDF (HTML page saved by an old fetch.py): delete it"
                                                    + (" (optional)" if fn in optional else "")))
    t = find(RAW / "thermoml", "*.tgz")
    if not t:
        add("MISSING", "raw", "ThermoML archive", "python fetch.py thermoml")
    else:
        good = sha256(t) == TGZ_SHA
        add("OK" if good else "FAIL", "raw", f"ThermoML {t.name}", f"{t.stat().st_size/1e6:.0f} MB" + ("" if good else " sha256 mismatch - refetch"))


# ------------------------------------------------------------------ C. interim (unzipped benchmark)
EXPECT_INTERIM = {   # filename -> expected rows (from the published inventory)
    "flashpoint_pure_exp.csv": 585, "flashpoint_mix_exp.csv": 1006,
    "viscosity_pure_exp.csv": 3027, "viscosity_mix_exp.csv": 35938,
    "dcn_mix_exp.csv": 484, "mon_mix_exp.csv": None, "NIST_Visc_Data.csv": 239200,
}


def interim():
    if not INT.exists():
        add("MISSING", "interim", "data/interim", "python inspect_downloads.py")
        return
    for fn, exp in EXPECT_INTERIM.items():
        p = find(INT, fn)
        if not p:
            add("MISSING", "interim", fn, "python inspect_downloads.py")
            continue
        n = nrows(p)
        status = "OK" if exp is None or abs(n - exp) <= max(2, 0.001 * exp) else "WARN"
        add(status, "interim", fn, f"{n} rows" + ("" if status == "OK" else f" (expected ~{exp})"))
    j = sum(1 for p in INT.rglob("*") if p.is_file() and junk(p.relative_to(INT)))
    if j:
        add("WARN", "interim", "macOS junk files on disk", f"{j} files; rerun inspect_downloads.py (new version removes them)")
    e = INT / "read_errors.txt"
    if e.exists():
        txt = e.read_text(errors="replace")
        m = re.findall(r"^(\w+): (\d+)", txt, re.M)
        add("OK" if not m else "WARN", "interim", "inspect read errors", ", ".join(f"{k}:{v}" for k, v in m) or "none")


# ------------------------------------------------------------------ D. processed outputs
def processed():
    import pandas as pd
    for fn, label in (("viscosity_pure_hc.csv", "NIST viscosity (pure, hydrocarbon)"),
                      ("viscosity_mix_hc.csv", "NIST viscosity (mixtures, hydrocarbon)"),
                      ("compounds_hc.csv", "hydrocarbon compound table")):
        p = PROC / fn
        add("OK" if p.exists() and nrows(p) > 0 else "MISSING", "processed", label,
            f"{nrows(p)} rows" if p.exists() else "python build_viscosity.py")
    p = PROC / "thermoml_hc.csv"
    if not p.exists():
        add("MISSING", "processed", "thermoml_hc.csv", "python thermoml_extract.py --debug   (needs the archive)")
    else:
        df = pd.read_csv(p)
        if df.empty:
            add("FAIL", "processed", "thermoml_hc.csv is empty", "run with --debug and send me the printed XML structure")
        else:
            if "prop_phase" not in df.columns:
                add("WARN", "processed", "thermoml_hc.csv is from the OLD extractor",
                    "it has excess viscosities, gas-phase densities, high pressures and VLE vapour compositions mixed in and lost "
                    "61% of the temperatures -> rerun: python thermoml_extract.py   (build_unified.py ignores the old file)")
            g = df.groupby(["prop", "n_comp"]).size()
            add("OK", "processed", "thermoml_hc.csv", f"{len(df)} rows; props: {df['prop'].nunique()}; " +
                "; ".join(f"{k[0][:28]}/{k[1]}c={v}" for k, v in g.head(12).items()))
            nT = df.T_K.notna().mean()
            add("OK" if nT > 0.8 else "WARN", "processed", "thermoml: rows with temperature", f"{nT:.0%}")
    a = PROC / "antoine_best.csv"
    if not a.exists():
        add("MISSING", "processed", "antoine_best.csv", "python webbook_antoine.py")
        return
    b = pd.read_csv(a)
    want = pd.read_csv("compounds.csv") if Path("compounds.csv").exists() else pd.DataFrame(columns=["name"])
    miss = sorted(set(want.name) - set(b.name))
    add("OK" if not miss else "WARN", "processed", "Antoine compounds found", f"{len(b)}/{len(want)}" + (f"; missing: {miss}" if miss else ""))
    nc = b[~b.covers_window.astype(bool)]
    add("OK" if nc.empty else "WARN", "processed", "Antoine range covers 40-100 C",
        "all" if nc.empty else f"{len(nc)} compounds do NOT: {list(nc.name)} -> refit from vapour-pressure data")
    norm = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower().replace("n-", "").replace("iso", ""))
    odd = [f"{r['name']}->{r['page_name']}" for _, r in b.iterrows() if norm(r["name"]) not in norm(r["page_name"])
           and norm(r["page_name"]) not in norm(r["name"])]
    add("OK" if not odd else "WARN", "processed", "WebBook page name vs your name (heuristic)",
        "all match" if not odd else f"check by eye: {odd[:12]}")


# ------------------------------------------------------------------ E. what is actually usable (hydrocarbon coverage)
def coverage():
    try:
        import pandas as pd
        from common import is_hydrocarbon, mol_from
    except Exception as e:  # noqa: BLE001
        add("WARN", "coverage", "skipped", str(e)); return
    cache: dict[str, bool] = {}

    def hc(inchi):
        if inchi not in cache:
            cache[inchi] = is_hydrocarbon(mol_from(inchi=inchi))
        return cache[inchi]

    for fn, label in (("flashpoint_pure_exp.csv", "flash point, pure"), ("flashpoint_mix_exp.csv", "flash point, mixtures"),
                      ("viscosity_pure_exp.csv", "viscosity, pure"), ("viscosity_mix_exp.csv", "viscosity, mixtures"),
                      ("dcn_mix_exp.csv", "DCN, mixtures"), ("mon_mix_exp.csv", "MON, mixtures")):
        p = find(INT, fn)
        if not p:
            continue
        df = pd.read_csv(p, low_memory=False)
        cols = [c for c in df.columns if c.startswith("component_inchi_")]
        n_comp = df[cols].notna().sum(axis=1)
        is_hc = df[cols].apply(lambda r: all(hc(x) for x in r.dropna()), axis=1)
        det = f"{int(is_hc.sum())}/{len(df)} rows all-hydrocarbon"
        if is_hc.any():
            det += f"; n_comp among HC rows: {n_comp[is_hc].value_counts().sort_index().to_dict()}"
            if "temperature" in df:
                t = df.loc[is_hc, "temperature"]; det += f"; T {t.min():.0f}-{t.max():.0f} K"
            if "value" in df:
                v = df.loc[is_hc, "value"]; det += f"; value {v.min():.3g} to {v.max():.3g} {df['unit'].iloc[0]}"
        add("OK" if is_hc.sum() > 0 else "WARN", "coverage", label, det)


# ------------------------------------------------------------------ F. manual / gated
def manual():
    for sub, what in (("ajftd", "AJFTD jet-fuel export (needs registration)"), ("rp1", "NIST RP-1 11-blend data"),
                      ("pour_point", "pour/cloud point tables (hand-extracted)"), ("lfl", "LFL values (from Zabetakis)"),
                      ("flashpoint_li", "Li et al. flash-point SI")):
        d = RAW / sub
        has = d.exists() and any(p for p in d.rglob("*") if p.is_file() and "template" not in p.name.lower())
        add("OK" if has else "MISSING", "manual", what, f"data/raw/{sub}/" + ("" if has else "  (empty - expected until you add it)"))


def unified():
    import pandas as pd
    U = Path("data/unified")
    dd, comp = U / "measurements_dedup.csv", U / "compounds.csv"
    if not dd.exists():
        add("MISSING", "unified", "measurements_dedup.csv", "python build_unified.py")
        return
    d = pd.read_csv(dd); c = pd.read_csv(comp)
    hc = d[d.is_hc]
    for prop, g in hc.groupby("property"):
        kinds = g.kind.value_counts().to_dict()
        add("OK", "unified", f"HC records: {prop}", f"{len(g)} ({kinds}); T {g.T_K.min():.0f}-{g.T_K.max():.0f} K; systems={g.system_id.nunique()}")
    m3 = hc[hc.n_comp >= 3]
    add("OK" if len(m3) else "WARN", "unified", "multi-component (>=3) hydrocarbon records", f"{len(m3)}" + ("" if len(m3) else " - need AJFTD/RP-1 for these"))
    add("OK" if not d.conflict.mean() > 0.05 else "WARN", "unified", "conflicting duplicate records", f"{int(d.conflict.sum())}/{len(d)}")
    if "antoine_covers_window" in c.columns:
        nw = int(c.antoine_covers_window.astype(bool).sum())
        add("OK" if nw else "WARN", "unified", "compounds with Antoine valid at 40-100 C", f"{nw}/{len(c)}"
            + ("" if nw else "  -> run thermoml_extract.py, build_unified.py, refit_antoine.py, build_unified.py"))
    c = c.copy()
    for col in ("lfl_vol_pct",):
        if col not in c.columns:
            c[col] = float("nan")
    nl = int(c.lfl_vol_pct.notna().sum())
    if "lfl_source" in c.columns and nl:
        meas = int(c.lfl_source.fillna("").str.startswith("measured").sum())
        add("OK" if meas else "WARN", "unified", "compounds with LFL", f"{nl}/{len(c)} (measured {meas}, estimated {nl - meas}) "
            + ("" if meas else " - all estimated; add measured values to data/raw/lfl/lfl_measured.csv"))
    else:
        add("OK" if nl else "MISSING", "unified", "compounds with LFL", f"{nl}/{len(c)}" + ("" if nl else "  -> python fill_lfl.py, then build_unified.py"))
    ood = U / "measurements_out_of_domain.csv"
    if ood.exists():
        add("OK", "unified", "rows removed by fuel-domain filter", f"{len(pd.read_csv(ood))} (see measurements_out_of_domain.csv)")
    bp = U / "BASELINE_FLASH_REPORT.md"
    add("OK" if bp.exists() else "MISSING", "unified", "flash-point physics baseline report", "" if bp.exists() else "python baselines_flash.py")
    add("OK" if (d.kind == "real_fuel").any() else "MISSING", "unified", "real-fuel records (AJFTD)", f"{int((d.kind == 'real_fuel').sum())}"
        + ("" if (d.kind == "real_fuel").any() else "  -> python ingest_ajftd.py --describe"))


def main():
    env(); raw(); interim()
    try:
        processed()
    except Exception as e:  # noqa: BLE001
        add("FAIL", "processed", "audit crashed", repr(e))
    try:
        coverage()
    except Exception as e:  # noqa: BLE001
        add("FAIL", "coverage", "audit crashed", repr(e))
    try:
        unified()
    except Exception as e:  # noqa: BLE001
        add("FAIL", "unified", "audit crashed", repr(e))
    manual()

    order = {"FAIL": 0, "MISSING": 1, "WARN": 2, "OK": 3}
    cnt = {k: sum(1 for c in checks if c[0] == k) for k in order}
    md = ["# Audit\n", " | ".join(f"{k}: {v}" for k, v in cnt.items()), "\n## Needs attention\n"]
    todo = [c for c in checks if c[0] != "OK"]
    md += [f"- **{s}** [{a}] {i}: {d}" for s, a, i, d in sorted(todo, key=lambda c: (order[c[0]], c[1]))] or ["- nothing"]
    md.append("\n## Everything\n\n| status | area | item | detail |\n|---|---|---|---|")
    md += [f"| {s} | {a} | {i} | {str(d).replace('|', '/')} |" for s, a, i, d in checks]
    Path("AUDIT.md").write_text("\n".join(md), encoding="utf-8")
    print("\n".join(md[:2] + md[2:3]))
    for s, a, i, d in sorted(todo, key=lambda c: (order[c[0]], c[1])):
        print(f"{s:8} [{a}] {i}: {d}")
    print("\nfull report: AUDIT.md")

    if "--bundle" in sys.argv:
        with zipfile.ZipFile("audit_bundle.zip", "w", zipfile.ZIP_DEFLATED) as z:
            for f in ("AUDIT.md", "data/interim/INVENTORY.md", "data/interim/read_errors.txt", "data/raw/MANIFEST.json",
                      "data/raw/fetch_failures.txt", "data/processed/antoine_best.csv", "data/processed/antoine_all_sets.csv",
                      "data/processed/antoine_refit.csv", "data/processed/compounds_hc.csv", "data/unified/UNIFIED_REPORT.md", "data/unified/BASELINE_FLASH_REPORT.md", "data/raw/lfl/lfl.csv", "data/processed/ajftd_fuels.csv",
                      "data/raw/ajftd/ajftd_metrics_inventory.csv",
                      "data/interim/ajftd_description.txt", "data/interim/ajftd_unmapped.txt", "data/raw/lfl/lfl_template.csv"):
                if Path(f).exists():
                    z.write(f)
            p = PROC / "thermoml_hc.csv"
            if p.exists():
                z.writestr("thermoml_hc_head.csv", "".join(p.open(encoding="utf-8").readlines()[:200]))
        print("wrote audit_bundle.zip  (small; contains no raw datasets) - upload it here")


if __name__ == "__main__":
    main()
