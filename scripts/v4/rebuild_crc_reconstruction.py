"""Rebuild the 57 CRC 647 compositions with realistic representatives (r2).  Descriptors only; no target is read.

Descriptors per fuel: D86 (mean of tests 1 and 2: IBP,T10,T20,T50,T90,FBP), aromatics vol% (D1319, else D6379),
olefins vol% (D1319, else D1159), diaromatics = naphthalenes vol% (D1840).
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.WARNING)

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.baselines.flash_unifac import _lookup  # noqa: E402
from fuelblend.data.reconstruction_r2 import KINDS, RealisticReps, diaromatic_mass, nominal_sat_split_detail, reconstruct  # noqa: E402
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider  # noqa: E402

PROPS = ROOT / "data/real_fuels/v4/jet_properties_v4.csv"
JETCLASS = ROOT / "data/real_fuels/v4/jet_class_compositions.csv"
OUT_CSV = ROOT / "data/real_fuels/v4/compositions_crc647_r2.csv"
OUT_MD = ROOT / "results/v4/R2/reconstruction_report.md"
D86_MAP = {"distillation_T0_C": "d86_ibp", "distillation_T10_C": "d86_t10", "distillation_T20_C": "d86_t20",
           "distillation_T50_C": "d86_t50", "distillation_T90_C": "d86_t90", "distillation_T100_C": "d86_fbp"}
COLS = ["system_id", "inchikey", "smiles", "fraction", "basis", "is_lumped", "unexplained_fraction", "provenance",
        "verified"]


def descriptor_rows(p: pd.DataFrame) -> dict:
    p = p[p.fuel_id.str.startswith("CRC647") & (p.role == "descriptor")]
    out = {}
    for fid, g in p.groupby("fuel_id"):
        row = {}
        d86 = g[g.method.str.startswith("ASTM D86")]
        for prop, key in D86_MAP.items():
            v = d86.loc[d86.property == prop, "value"]
            if len(v):
                row[key] = float(v.mean())

        def pick(prop, prefs):
            for m in prefs:
                v = g.loc[(g.property == prop) & (g.method == m), "value"]
                if len(v):
                    return float(v.mean())
            return np.nan
        row["aromatics_vol_pct"] = pick("aromatics_vol_pct", ["ASTM D1319", "ASTM D6379"])
        row["olefins_vol_pct"] = pick("olefins_vol_pct", ["ASTM D1319", "ASTM D1159"])
        row["diaromatics_vol_pct"] = pick("naphthalenes_vol_pct", ["ASTM D1840"])
        out[fid] = row
    return out


def build_compositions(rows_by_fuel, split, reps):
    """-> (out_rows, stats[(fid, fallback, n_rows)], kind_mass{kind: [per-fuel share]}, skipped_fuels, diaro)."""
    out, stats, kind_mass, skipped_fuels, diaro = [], [], {k: [] for k in KINDS}, [], []
    for fid in sorted(rows_by_fuel):
        fr, fb = reconstruct(rows_by_fuel[fid], split, reps)
        if fr.empty and fb == 1.0:
            skipped_fuels.append(fid)
            continue
        diaro.append(diaromatic_mass(rows_by_fuel[fid]))
        stats.append((fid, fb, len(fr)))
        for k in KINDS:
            kind_mass[k].append(float(fr.loc[fr.kind == k, "fraction"].sum()) if len(fr) else 0.0)
        for r in fr.itertuples():
            out.append(dict(system_id=f"fuel:{fid}", inchikey=r.inchikey, smiles=r.smiles, fraction=r.fraction,
                            basis="mass", is_lumped=False, unexplained_fraction=fb,
                            provenance="reconstructed_r2", verified=False))
        out.append(dict(system_id=f"fuel:{fid}", inchikey="LUMP_NO_REALISTIC_REP", smiles="", fraction=fb,
                        basis="mass", is_lumped=True, unexplained_fraction=fb,
                        provenance="reconstructed_r2", verified=False))
    return out, stats, kind_mass, skipped_fuels, diaro


def main():
    reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
    reps = RealisticReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    split, used, skipped = nominal_sat_split_detail(JETCLASS)
    rows_by_fuel = descriptor_rows(pd.read_csv(PROPS))

    out, stats, kind_mass, skipped_fuels, diaro = build_compositions(rows_by_fuel, split, reps)
    df = pd.DataFrame(out, columns=COLS)
    fbs = np.array([s[1] for s in stats])
    print(f"skipped={skipped_fuels}")
    print(f"fuels={len(stats)} mean_fallback={fbs.mean():.4f} max_fallback={fbs.max():.4f} split={split}")
    if fbs.mean() > 0.20:
        print("STOP: mean fallback > 0.20; not writing CSV. Bins without a realistic representative:")
        for k in KINDS:
            print(k, "available C#:", reps.ns[k])
        sys.exit(2)
    df.to_csv(OUT_CSV, index=False)

    md = ["# r2 CRC 647 reconstruction report", "",
          "Descriptors only (D86 mean of tests 1 and 2, aromatics, olefins, naphthalenes); no target property is read.", "",
          f"CRC fuels written: {len(stats)} (of {len(rows_by_fuel)}); skipped: "
          + (", ".join(skipped_fuels) if skipped_fuels else "none")
          + " (fewer than 4 D86 points or no aromatics; excluded from the CSV, the fallback stats and the shares below)", "",
          "## Nominal saturates split (prior from GC jet fuels, not a fit to CRC)", "",
          f"n-paraffin {split[0]:.4f}, iso-paraffin {split[1]:.4f}, cycloalkane {split[2]:.4f}", "",
          f"GC fuels used ({len(used)}): " + ", ".join(used), ""]
    if skipped:
        md += ["Skipped: " + "; ".join(f"{a} ({b})" for a, b in skipped), ""]
    md += ["## Mean mass share per kind (over the CRC fuels)", "", "| kind | mean share |", "|---|---|"]
    md += [f"| {k} | {np.mean(v):.4f} |" for k, v in kind_mass.items()]
    md += ["| fallback (LUMP_NO_REALISTIC_REP) | %.4f |" % fbs.mean(), "",
           f"Fallback mass: mean {fbs.mean():.4f}, max {fbs.max():.4f} (worst fuel {stats[int(fbs.argmax())][0]})", "",
           "## Representative per (kind, carbon number) bin", "", "| kind | C# | inchikey | smiles |", "|---|---|---|---|"]
    for (k, n), v in sorted(reps.by.items()):
        md.append(f"| {k} | {n} | {reps.ik[v[0]]} | {reps.smiles[v[0]]} |")
    md += ["", "Candidates also require an Antoine set and LFL in the flash registry and valid pure eta/rho.", ""]
    # known gaps
    from rdkit import Chem
    empty_kinds = [k for k in KINDS if not reps.ns[k]]
    bicyc = [i for (k, n), v in reps.by.items() if k == "cycloalkane" for i in v
             if Chem.MolFromSmiles(reps.smiles[i]).GetRingInfo().NumRings() == 2]
    dmean = float(np.mean(diaro)) if diaro else 0.0
    md += ["## Known gaps", "",
           "- Kinds with no realistic representative at all: " + (", ".join(empty_kinds) if empty_kinds else "none"),
           "- Bicyclic cycloalkane (decalin-type) candidates available: %d" % len(bicyc)
           + (" (none: all cycloalkane mass sits on single-ring compounds)" if not bicyc else ""),
           "- aromatic_di representatives available: %d bin(s) (C# %s)" % (len(reps.ns["aromatic_di"]), reps.ns["aromatic_di"]),
           f"- Diaromatic (naphthalene) mass per fuel, mean {dmean:.4f}, max {max(diaro):.4f}; "
           + ("all of it" if not reps.ns["aromatic_di"] else "the part without a representative")
           + f" goes to the fallback lump, i.e. {dmean / fbs.mean() * 100:.0f}% of the mean fallback mass ({fbs.mean():.4f}).",
           f"- The nominal split is cyclo-heavy (cycloalkane {split[2]:.3f}) and comes from only {len(used)} GC fuels; "
           "it is a prior, not a fit to CRC.",
           "- get() substitutes the nearest-carbon-number representative within max_gap=3 carbons; the table below shows "
           "the requested carbon number, the representative's carbon number and the gap (gap > 0 = substituted).", "",
           "## Carbon-number substitution (kinds other than aromatic_di)", "",
           "| kind | requested C# | representative C# | gap |", "|---|---|---|---|"]
    for k in KINDS:
        if k == "aromatic_di":
            continue
        for n in range(5, 21):
            i = reps.get(k, n)
            if i is None:
                md.append(f"| {k} | {n} | none | - |")
                continue
            rn = int(sum(a.GetSymbol() == "C" for a in Chem.MolFromSmiles(reps.smiles[i]).GetAtoms()))
            if rn != n:
                md.append(f"| {k} | {n} | {rn} | {abs(rn - n)} |")
    md += ["", "(Rows are listed only where the requested bin is empty: substituted or without representative.)"]
    OUT_MD.write_text("\n".join(md) + "\n", encoding="utf8")
    print("wrote", OUT_CSV, OUT_MD)


if __name__ == "__main__":
    main()
