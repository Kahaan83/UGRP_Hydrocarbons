"""v5 L1: re-run the flash evaluation with the LFL overlay (data/processed/v5/lfl_v5.csv) applied in memory.

Readings are declared in results/v5/L1/README.md.  The v4 evaluation code is run unchanged; only the pure-compound
registry returned by fuelblend.baselines.flash_unifac._default_registry is replaced by a copy whose lfl_vol_pct comes
from the overlay.  No file under results/v4, configs/ or data/processed/v2 is written.

Outputs: results/v5/L1/F_r3/* (v4 flash outputs on the overlay), ajftd_g1_flash_overlay.csv, comparison.csv, RESULTS.md
"""
from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import fuelblend.baselines.flash_unifac as fu  # noqa: E402

OUT = ROOT / "results" / "v5" / "L1"
OVERLAY = ROOT / "data" / "processed" / "v5" / "lfl_v5.csv"
R3C = ROOT / "data" / "real_fuels" / "v4" / "compositions_crc647_r3.csv"
BASE_SUMMARY = ROOT / "results" / "v4" / "R3" / "F" / "flash_summary.csv"
G1 = ROOT / "results" / "v5" / "G1" / "predictions.csv"
MAPPED = ROOT / "data" / "real_fuels" / "v5" / "ajftd_gcxgc_mapped.csv"
TABULATED = ("iec60079_20_1_table", "cooptima_compiled")


def overlaid_registry() -> tuple[pd.DataFrame, pd.DataFrame]:
    """-> (registry copy with overlay LFL, table of the changed compounds)."""
    reg = fu._default_registry().copy()
    ov = pd.read_csv(OVERLAY)
    ov = ov[ov.lfl_source_v5.isin(TABULATED)].set_index("inchikey")
    hit = reg.inchikey.isin(ov.index) & reg.lfl_vol_pct.notna()
    ch = reg.loc[hit, ["inchikey", "lfl_vol_pct", "lfl_source"]].rename(columns={"lfl_vol_pct": "lfl_before", "lfl_source": "source_before"})
    ch["lfl_after"] = ch.inchikey.map(ov.lfl_vol_pct)
    ch["source_after"] = ch.inchikey.map(ov.lfl_source_v5)
    still_estimated = ch.source_before.astype(str).str.startswith("estimated")
    ch = ch[still_estimated]                                 # never overwrite a measured value
    reg.loc[ch.index, "lfl_vol_pct"] = ch.lfl_after
    reg.loc[ch.index, "lfl_source"] = "v5 overlay: " + ch.source_after
    return reg, ch.reset_index()


def stats(e) -> dict:
    e = np.asarray(e, float)
    e = e[np.isfinite(e)]
    return dict(N=int(e.size), MAE_C=float(np.abs(e).mean()), bias_C=float(e.mean()), over_gt5C_rate=float((e > 5).mean()))


def ajftd_flash() -> pd.DataFrame:
    g1 = pd.read_csv(G1)
    g1 = g1[g1.property == "flash_point_C"]
    mapped = pd.read_csv(MAPPED)
    rows = []
    for fid, t in g1.groupby("fuel_id"):
        g = mapped[mapped.fuel_id == fid]
        smi, w = list(g.smiles), g.fraction.to_numpy(float)
        x = w / np.array([Descriptors.MolWt(Chem.MolFromSmiles(s)) for s in smi])
        x = x / x.sum()
        for model, use in (("ideal", False), ("unifac", True)):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", fu.UnifacFallbackWarning)
                p = fu.flash_point_unifac_detail(smi, x, use_unifac=use)
            b = t[t.model == model].iloc[0]
            rows.append(dict(fuel_id=fid, fuel_label=b.fuel_label, model=model, measured=b.measured, pred_before=b.pred,
                             pred_after=p["flash_C"], err_before=b.pred - b.measured, err_after=p["flash_C"] - b.measured))
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    reg, changed = overlaid_registry()
    changed.to_csv(OUT / "registry_lfl_changes.csv", index=False)
    fu._default_registry = lambda: reg                       # in-memory only; _lookup resolves the module global at call time

    spec = importlib.util.spec_from_file_location("eval_flash_v4", ROOT / "scripts" / "v4" / "eval_flash_v4.py")
    ef = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ef)
    ef._default_registry = fu._default_registry
    ef.main(["--compositions", str(R3C), "--out", str(OUT / "F_r3")])

    before = pd.read_csv(BASE_SUMMARY)
    after = pd.read_csv(OUT / "F_r3" / "flash_summary.csv")
    cmp = before.merge(after, on=["scope", "slice", "model"], suffixes=("_before", "_after"))
    cmp = cmp[["scope", "slice", "model", "N_before", "MAE_C_before", "MAE_C_after", "bias_C_before", "bias_C_after",
               "over_gt5C_rate_before", "over_gt5C_rate_after"]].rename(columns={"N_before": "N"})
    aj = ajftd_flash()
    aj.to_csv(OUT / "ajftd_g1_flash_overlay.csv", index=False)
    for model, g in aj.groupby("model"):
        b, a = stats(g.err_before), stats(g.err_after)
        cmp.loc[len(cmp)] = ["ajftd_gcxgc_G1", "all eligible", model, b["N"], b["MAE_C"], a["MAE_C"], b["bias_C"], a["bias_C"],
                             b["over_gt5C_rate"], a["over_gt5C_rate"]]
    cmp.to_csv(OUT / "comparison.csv", index=False)

    used_r3 = pd.read_csv(OUT / "F_r3" / "pure_inputs_used.csv")
    n_r3 = int(used_r3.lfl_source.astype(str).str.startswith("v5 overlay").sum())
    mapped_smiles = {Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in pd.read_csv(MAPPED).smiles.unique()}
    n_aj = len(mapped_smiles & set(changed.can))
    key = [("lab_mixtures", "all_126 (n_systems=31)"), ("real_fuels_CRC647", "D93 dev"), ("ajftd_gcxgc_G1", "all eligible")]
    sel = [cmp[(cmp.scope == s) & (cmp.slice == sl) & (cmp.model == "ideal")] for s, sl in key]
    missing = [k for k, d in zip(key, sel) if d.empty]
    L = ["# L1 results (generated by scripts/v5/eval_flash_lfl_overlay_v5.py; readings are in README.md)", "",
         f"Registry compounds whose estimated LFL was replaced: {len(changed)} "
         f"({', '.join(f'{k} {v}' for k, v in changed.source_after.value_counts().items())}). "
         f"Of the pure inputs actually used: {n_r3} of {len(used_r3)} in run A, {n_aj} of {len(mapped_smiles)} in run B.", "",
         "| scope | slice | model | N | MAE before | MAE after | bias before | bias after | over >5 C before | after |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in cmp.itertuples():
        L.append(f"| {r.scope} | {r.slice} | {r.model} | {int(r.N)} | {r.MAE_C_before:.3f} | {r.MAE_C_after:.3f} | {r.bias_C_before:+.3f} | "
                 f"{r.bias_C_after:+.3f} | {r.over_gt5C_rate_before:.3f} | {r.over_gt5C_rate_after:.3f} |")
    L += ["", "## Readings", ""]
    if missing:
        L.append(f"* R2 and R3 not evaluable: slice(s) {missing} not found in the summaries.")
    else:
        d = pd.concat(sel)
        dm = (d.MAE_C_after - d.MAE_C_before).to_numpy()
        do = (d.over_gt5C_rate_after - d.over_gt5C_rate_before).to_numpy()
        r2 = bool((do <= 1e-12).all() and (dm <= 0.2).all())
        r3 = bool((dm < -0.2).sum() >= 2 and (dm <= 1e-12).all())
        L.append(f"* Ideal-model MAE change (lab 126, CRC 647 D93 dev, AJFTD G1): {', '.join(f'{v:+.3f}' for v in dm)} C; "
                 f"change in over-by-5 C fraction: {', '.join(f'{v:+.3f}' for v in do)}.")
        L.append(f"* R2 \"the overlay is safe to adopt\": {'HOLDS' if r2 else 'DOES NOT HOLD'}.")
        L.append(f"* R3 \"the overlay improves accuracy\": {'HOLDS' if r3 else 'DOES NOT HOLD'}.")
    L.append("")
    (OUT / "RESULTS.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
