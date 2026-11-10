"""Rebuild the 57 CRC 647 compositions with Table G.1 hydrocarbon types (r3).  Descriptors only; no target is read.

Inputs per fuel: ASTM D6379 aromatics/saturates and ASTM D2425 type shares (CRC 647 Table G.1, parsed to
data/raw/pour_point/crc647_hydrocarbon_types_g1.csv), D86 (mean of tests 1 and 2), olefins vol% (D1319, else D1159).
"""
from __future__ import annotations

import importlib.util
import json
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.WARNING)

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.baselines.flash_unifac import _lookup  # noqa: E402
from fuelblend.data.reconstruction_r2 import nominal_sat_split_detail  # noqa: E402
from fuelblend.data.reconstruction_r3 import (ARO_TYPES, SAT_TYPES, SUBKINDS, UNEXPLAINED, SubReps,  # noqa: E402
                                              n_share_prior, reconstruct_r3, subkind)
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider  # noqa: E402

_spec = importlib.util.spec_from_file_location("rebuild_crc_reconstruction_r2",
                                               ROOT / "scripts/v4/rebuild_crc_reconstruction.py")
_r2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_r2)
descriptor_rows, COLS = _r2.descriptor_rows, _r2.COLS

PROPS = ROOT / "data/real_fuels/v4/jet_properties_v4.csv"
JETCLASS = ROOT / "data/real_fuels/v4/jet_class_compositions.csv"
G1_CSV = ROOT / "data/raw/pour_point/crc647_hydrocarbon_types_g1.csv"
R2_CSV = ROOT / "data/real_fuels/v4/compositions_crc647_r2.csv"
OUT_CSV = ROOT / "data/real_fuels/v4/compositions_crc647_r3.csv"
OUT_DIR = ROOT / "results/v4/R3"
OUT_MD = OUT_DIR / "reconstruction_report.md"
OUT_FB = OUT_DIR / "fallback_detail.csv"
LUMP = "LUMP_NO_REALISTIC_REP"
R2_FIXED_RING_SHARE = 0.5260


def g1_code(fuel_id: str) -> str:
    m = re.fullmatch(r"CRC647_(\d{1,3})", fuel_id)
    if m is None:
        raise ValueError(f"cannot map fuel id {fuel_id!r} to a Table G.1 code")
    return m.group(1).zfill(3)


def load_g1(path: Path, fuel_ids) -> dict:
    """{fuel_id: {property: value}}; cells printed with '<' (below detection) are 0.  Stops on a missing/ambiguous map."""
    df = pd.read_csv(path, dtype={"fuel_code": str})
    wide: dict[str, dict] = {}
    for r in df.itertuples(index=False):
        v = 0.0 if str(r.value_raw).startswith("<") else float(r.value_std)
        wide.setdefault(r.fuel_code, {})[r.property] = v
    codes = {fid: g1_code(fid) for fid in fuel_ids}
    dup = {c for c in codes.values() if list(codes.values()).count(c) > 1}
    missing = sorted(fid for fid, c in codes.items() if c not in wide)
    if dup or missing:
        raise SystemExit(f"G.1 mapping problem: ambiguous codes {sorted(dup)}, fuels without a G.1 row {missing}")
    return {fid: wide[c] for fid, c in codes.items()}


def build_compositions(rows_by_fuel, g1_by_fuel, n_share, reps):
    """-> (out_rows, stats[(fid, fallback, n_rows)], kind_mass{subkind: [per-fuel share]}, skipped_fuels, details{fid: dict})."""
    out, stats, kind_mass, skipped_fuels, details = [], [], {k: [] for k in SUBKINDS}, [], {}
    for fid in sorted(rows_by_fuel):
        fr, fb, det = reconstruct_r3(rows_by_fuel[fid], g1_by_fuel[fid], n_share, reps)
        if fr.empty and fb == 1.0:
            skipped_fuels.append(fid)
            continue
        stats.append((fid, fb, len(fr)))
        details[fid] = det
        for k in SUBKINDS:
            kind_mass[k].append(float(fr.loc[fr.kind == k, "fraction"].sum()) if len(fr) else 0.0)
        for r in fr.itertuples():
            out.append(dict(system_id=f"fuel:{fid}", inchikey=r.inchikey, smiles=r.smiles, fraction=r.fraction,
                            basis="mass", is_lumped=False, unexplained_fraction=fb,
                            provenance="reconstructed_r3", verified=False))
        out.append(dict(system_id=f"fuel:{fid}", inchikey=LUMP, smiles="", fraction=fb,
                        basis="mass", is_lumped=True, unexplained_fraction=fb,
                        provenance="reconstructed_r3", verified=False))
    return out, stats, kind_mass, skipped_fuels, details


def _kind_shares(csv_path: Path, fuels) -> dict:
    """Per-kind mean mass share over ``fuels`` from a composition CSV (kinds derived from SMILES with ``subkind``)."""
    df = pd.read_csv(csv_path, keep_default_na=False)
    df = df[df.system_id.isin({f"fuel:{f}" for f in fuels})]
    cache: dict[str, str] = {}

    def kind_of(row):
        if row.inchikey == LUMP:
            return "fallback lump"
        if row.smiles not in cache:
            cache[row.smiles] = subkind(row.smiles) or "other (not a sub-kind)"
        return cache[row.smiles]
    df = df.assign(k=[kind_of(r) for r in df.itertuples()])
    piv = df.pivot_table(index="system_id", columns="k", values="fraction", aggfunc="sum", fill_value=0.0)
    return piv.mean().to_dict(), df.system_id.nunique()


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-share", type=float, default=None, help="override the n-paraffin share of the paraffins (default: GC prior)")
    ap.add_argument("--out-csv", type=Path, default=globals()["OUT_CSV"])
    ap.add_argument("--out-dir", type=Path, default=globals()["OUT_DIR"], help="directory for reconstruction_report.md and fallback_detail.csv")
    args = ap.parse_args(argv)
    OUT_CSV, OUT_DIR = args.out_csv, args.out_dir
    OUT_MD, OUT_FB = OUT_DIR / "reconstruction_report.md", OUT_DIR / "fallback_detail.csv"
    reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    n_share = n_share_prior(JETCLASS) if args.n_share is None else args.n_share
    _, gc_used, _ = nominal_sat_split_detail(JETCLASS)
    rows_by_fuel = descriptor_rows(pd.read_csv(PROPS))
    g1 = load_g1(G1_CSV, sorted(rows_by_fuel))

    out, stats, kind_mass, skipped_fuels, details = build_compositions(rows_by_fuel, g1, n_share, reps)
    df = pd.DataFrame(out, columns=COLS)
    fbs = np.array([s[1] for s in stats])
    print(f"skipped={skipped_fuels}")
    print(f"fuels={len(stats)} mean_fallback={fbs.mean():.4f} max_fallback={fbs.max():.4f} n_share={n_share:.4f}")
    if fbs.mean() > 0.20:
        print("STOP: mean fallback > 0.20; not writing CSV. Bins without a realistic representative:")
        for k in SUBKINDS:
            print(k, "available C#:", reps.ns[k])
        sys.exit(2)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    fuels = [s[0] for s in stats]
    pd.DataFrame([dict(fuel_id=f, bucket=b, mass=v) for f in fuels for b, v in sorted(details[f].items())],
                 columns=["fuel_id", "bucket", "mass"]).to_csv(OUT_FB, index=False)

    # ---- report ---------------------------------------------------------------------------------------------
    buckets = sorted({b for f in fuels for b in details[f]})
    ring = np.array([(g1[f]["hc_type_cycloparaffins_pct"] + g1[f]["hc_type_dicycloparaffins_pct"]
                      + g1[f]["hc_type_tricycloparaffins_pct"]) / sum(g1[f][f"hc_type_{k}_pct"] for k in SAT_TYPES)
                     for f in fuels])
    d_aro = np.array([g1[f]["aromatics_vol_pct"] - rows_by_fuel[f]["aromatics_vol_pct"] for f in fuels])
    r2_share, r2_n = _kind_shares(R2_CSV, fuels) if R2_CSV.exists() else ({}, 0)
    r3_share, r3_n = _kind_shares(OUT_CSV, fuels)
    two_ring = ("cyclo_di", "indan_tetralin", "aromatic_di")

    md = ["# r3 CRC 647 reconstruction report", "",
          "Inputs per fuel: CRC 647 Table G.1 (ASTM D6379 aromatics and saturates in vol%, ASTM D2425 hydrocarbon-type "
          "shares), the D86 curve (mean of tests 1 and 2: IBP, T10, T20, T50, T90, FBP) and olefins vol% (D1319, else "
          "D1159). No target property is read.", "",
          f"CRC fuels written: {len(stats)} (of {len(rows_by_fuel)}); skipped: "
          + (", ".join(skipped_fuels) if skipped_fuels else "none"), "",
          "## n / iso paraffin split (prior)", "",
          f"n-paraffin share of the paraffins = {n_share:.4f} (n / (n + iso)). This is a prior from the same {len(gc_used)} "
          "GC jet fuels r2 used, not a fit to CRC; Table G.1 does not separate normal from branched paraffins.", "",
          "## Mass share per sub-kind (over the fuels written)", "", "| sub-kind | mean | min | max |", "|---|---|---|---|"]
    md += [f"| {k} | {np.mean(v):.4f} | {np.min(v):.4f} | {np.max(v):.4f} |" for k, v in kind_mass.items()]
    md += [f"| fallback (total) | {fbs.mean():.4f} | {fbs.min():.4f} | {fbs.max():.4f} |", "",
           "## Fallback", "",
           f"Fallback mass: mean {fbs.mean():.4f}, max {fbs.max():.4f} (worst fuel {stats[int(fbs.argmax())][0]}).", "",
           "Mean fallback mass by bucket (fuels without that bucket count as 0):", "", "| bucket | mean mass |", "|---|---|"]
    md += [f"| {b} | {np.mean([details[f].get(b, 0.0) for f in fuels]):.4f} |" for b in buckets]
    md += ["", "Per-fuel detail: results/v4/R3/fallback_detail.csv.", "",
           "## Representative per (sub-kind, carbon number) bin", "",
           "Candidates require a usable registry entry, valid pure eta/rho, and an Antoine set plus LFL in the flash "
           "registry. Lowest InChIKey within a bin; `n cand` > 1 marks bins that hold more than one candidate.", "",
           "| sub-kind | C# | n cand | inchikey | smiles |", "|---|---|---|---|---|"]
    for (k, n), v in sorted(reps.by.items(), key=lambda kv: (SUBKINDS.index(kv[0][0]), kv[0][1])):
        md.append(f"| {k} | {n} | {len(v)}{' (>1)' if len(v) > 1 else ''} | {reps.ik[v[0]]} | {reps.smiles[v[0]]} |")
    md += ["", "## Two-ring classes: carbon numbers present in the registry", ""]
    md += [f"- {k}: " + (", ".join(f"C{n}" for n in reps.ns[k]) if reps.ns[k] else "none") for k in two_ring]
    md += ["", "## Ring share of the saturates: measured (G.1) vs r2 fixed split", "",
           "Measured ring share = (cycloparaffins + dicycloparaffins + tricycloparaffins) / sum of the four D2425 "
           "saturate types.", "", "| | mean | min | max |", "|---|---|---|---|",
           f"| G.1 measured (r3 input, {len(fuels)} fuels) | {ring.mean():.4f} | {ring.min():.4f} | {ring.max():.4f} |",
           f"| r2 fixed split (cycloalkane) | {R2_FIXED_RING_SHARE:.4f} | {R2_FIXED_RING_SHARE:.4f} | "
           f"{R2_FIXED_RING_SHARE:.4f} |", "",
           "## Mean mass share per kind: r2 vs r3", "",
           f"Both computed from the composition CSVs, kinds derived from SMILES with `subkind` (r2 fuels found: {r2_n}; "
           f"r3 fuels found: {r3_n}).", "", "| kind | r2 mean | r3 mean |", "|---|---|---|"]
    for k in list(SUBKINDS) + ["other (not a sub-kind)", "fallback lump"]:
        md.append(f"| {k} | {r2_share.get(k, 0.0):.4f} | {r3_share.get(k, 0.0):.4f} |")
    md += ["", f"G.1 total aromatics minus the aromatics r2 used (D1319, else D6379), over the fuels: mean "
           f"{d_aro.mean():+.2f} vol%, min {d_aro.min():+.2f}, max {d_aro.max():+.2f}.", "",
           "## Assumptions and known gaps", "",
           f"- vol% to mass uses the class densities 0.87 (aromatic) / 0.78 (olefin) / 0.79 (saturate) g/cm3, the same as r2.",
           "- D2425 type shares are treated as mass shares inside each HPLC fraction (saturates, aromatics); the table prints "
           "no unit for them.",
           "- Cells below the detection limit (printed with '<') are set to 0.",
           "- The n / iso paraffin split is a prior (above); G.1 gives total paraffins only.",
           "- One D86-based carbon-number distribution is applied to every class.",
           "- cyclo_di, indan_tetralin and aromatic_di use the nearest carbon number available in the registry "
           "(max_gap=99), because the registry holds very few such molecules.",
           "- Tricycloparaffins and the minor aromatic types (indenes, acenaphthenes, acenaphthylenes, tricyclic aromatics) "
           "have no realistic representative class and go to named fallback buckets.",
           "- The G.1 (D6379) aromatics differ from the page-92 D1319/D6379 aromatics that r2 used (see the difference above), "
           "so the aromatic total also changes between r2 and r3.",
           "- r2's report line \"Bicyclic cycloalkane candidates available: 0\" is wrong: it counted rings on a temporary "
           "RDKit Mol, whose RingInfo reports 0 rings. r3 keeps the Mol in a variable."]
    OUT_MD.write_text("\n".join(md) + "\n", encoding="utf8")
    print("wrote", OUT_CSV, OUT_MD, OUT_FB)


if __name__ == "__main__":
    main()
