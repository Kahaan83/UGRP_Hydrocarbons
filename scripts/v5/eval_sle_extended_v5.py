"""v5 data audit: how many ThermoML SLE (liquidus) systems are unused only because fusion constants were never audited?

v4 (scripts/v4/eval_sle.py) predicts a liquidus only when EVERY component has a `verified` Tm and dHfus in
data/processed/v4/pure_props_v4.csv.  That table covers the design compounds, so 51 of 59 systems are never scored and the
freezing margin has 8 groups (9 needed).  This script keeps the v4 model unchanged (ideal Schroeder-van Laar, no solid
solution, no solid-solid transition term) and widens only the constants, in declared tiers:

  v4_verified        : the v4 table, unchanged (must reproduce the v4 system count).
  two_compilations   : + compounds whose melting point has a consensus of tabulated sources and whose dHfus (CRC table in the
                       `chemicals` library) agrees within 2 % with the Broadhurst experimental n-alkane table shipped in
                       extras/fuelblend_data/data/raw/wax/correlation_validation.xlsx.  The two compilations may share a
                       primary measurement, so this is agreement of compilations, not two independent experiments.
  single_compilation : + compounds with a melting-point consensus and a CRC-table dHfus that has no second source.

Rules fixed before scoring
  * Tm: tabulated sources CRC_ORG, WEBBOOK, OPEN_NTBKM, COMMON_CHEMISTRY only (never JOBACK or WIKIDATA).  At least two
    sources; all within 2.0 K, or all but one within 2.0 K (that one is dropped and listed).  Otherwise `disputed`: excluded.
  * dHfus: CRC table only.  An estimated (JOBACK) value is never used.
  * A component whose InChI leaves a stereocentre unassigned is excluded (its constants would belong to one isomer).
  * Nothing is fitted.  The margin is the same group-level split-conformal quantile as configs/uncertainty_v4.yaml
    (residual = predicted - observed, unsafe direction `under`, alpha 0.10, group = system_id).

Outputs: results/v5/SLE2/constants_extended.csv, sle_eval_extended.csv, summary_by_tier.csv, README.md
This script writes nothing to configs/ and changes no v4 output.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.inchi import InchiToInchiKey

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts" / "v4"))
import eval_sle as v4  # noqa: E402
from fuelblend.uq.conformal import conformal_rank, group_scores, one_sided_margin  # noqa: E402

OUT = ROOT / "results" / "v5" / "SLE2"
BROADHURST = ROOT / "extras" / "fuelblend_data" / "data" / "raw" / "wax" / "correlation_validation.xlsx"
TM_SOURCES = ("CRC_ORG", "WEBBOOK", "OPEN_NTBKM", "COMMON_CHEMISTRY")
TM_TOL_K = 2.0
DH_TOL_REL = 0.02
ALPHA = 0.10
TIERS = ("v4_verified", "two_compilations", "single_compilation")


def tm_consensus(vals: dict[str, float], tol: float = TM_TOL_K) -> tuple[float | None, str]:
    """Median of the tabulated melting points if they agree within `tol` (at most one source dropped), else None."""
    v = sorted(vals.items(), key=lambda kv: kv[1])
    if len(v) < 2:
        return None, "single_source"
    x = [b for _, b in v]
    if x[-1] - x[0] <= tol:
        return float(np.median(x)), "agree"
    if len(v) >= 3:
        for drop in (0, len(v) - 1):                      # only an extreme value can be the lone outlier
            y = x[:drop] + x[drop + 1:]
            if y[-1] - y[0] <= tol:
                return float(np.median(y)), f"agree_dropped_{v[drop][0]}={v[drop][1]:.2f}"
    return None, "disputed"


def broadhurst_dhfus() -> dict[int, float]:
    """carbon number -> experimental enthalpy of fusion, kJ/mol (n-alkanes), from the SP-Wax validation workbook."""
    if not BROADHURST.exists():
        return {}
    d = pd.read_excel(BROADHURST, sheet_name="Thermo-phsyical Properties", header=None)
    out = {}
    for _, r in d.iloc[4:36].iterrows():
        name, cn, val = str(r[7]).strip(), r[9], r[10]
        if name.startswith("C") and not name.startswith("iC") and pd.notna(cn) and pd.notna(val):
            out[int(cn)] = float(val)
    return out


def n_alkane_cn(mol) -> int | None:
    if any(a.GetSymbol() != "C" or a.GetDegree() > 2 for a in mol.GetAtoms()) or mol.GetRingInfo().NumRings():
        return None
    if any(b.GetBondType() != Chem.BondType.SINGLE for b in mol.GetBonds()):
        return None
    return mol.GetNumAtoms()


def unassigned_stereo(mol) -> bool:
    return any(c[1] == "?" for c in Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False))


def extended_constants(components: dict[str, tuple[str, str]]) -> pd.DataFrame:
    """components: full InChIKey -> (name, InChI).  One row per component with the tier it qualifies for (or the reason not)."""
    from chemicals import CAS_from_any
    from chemicals.phase_change import Hfus, Hfus_methods, Tm, Tm_methods

    v4c = v4.verified_constants()
    bh = broadhurst_dhfus()
    rows = []
    for ik, (name, inchi) in sorted(components.items(), key=lambda kv: kv[1][0]):
        mol = Chem.MolFromInchi(inchi)
        rec = dict(inchikey=ik, name=name, tier="", tm_K=np.nan, dh_J_mol=np.nan, tm_rule="", tm_sources="", dh_source="",
                   dh_broadhurst_kJ=np.nan, reason="")
        if ik[:14] in v4c:                                # v4 keys on the connectivity block; kept so tier 1 reproduces v4
            rec.update(tier="v4_verified", tm_K=v4c[ik[:14]]["tm_K"], dh_J_mol=v4c[ik[:14]]["dh_J_mol"],
                       tm_rule="v4", dh_source="pure_props_v4 verified")
            rows.append(rec)
            continue
        if unassigned_stereo(mol):
            rec["reason"] = "stereo_unspecified"
            rows.append(rec)
            continue
        try:
            cas = CAS_from_any(inchi)
        except Exception:
            try:
                cas = CAS_from_any(name)
            except Exception:
                rec["reason"] = "no_cas"
                rows.append(rec)
                continue
        tms = {m: float(Tm(cas, method=m)) for m in Tm_methods(cas) if m in TM_SOURCES}
        tm, rule = tm_consensus(tms)
        rec.update(tm_rule=rule, tm_sources=json.dumps({k: round(x, 2) for k, x in tms.items()}))
        if tm is None:
            rec["reason"] = f"tm_{rule}"
            rows.append(rec)
            continue
        if "CRC" not in Hfus_methods(cas):
            rec.update(tm_K=tm, reason="no_tabulated_dHfus")
            rows.append(rec)
            continue
        dh = float(Hfus(cas, method="CRC"))
        rec.update(tm_K=tm, dh_J_mol=dh, dh_source="chemicals CRC table", tier="single_compilation")
        cn = n_alkane_cn(mol)
        if cn in bh:
            rec["dh_broadhurst_kJ"] = bh[cn]
            if abs(dh / 1000.0 - bh[cn]) <= DH_TOL_REL * bh[cn]:
                rec["tier"] = "two_compilations"
        rows.append(rec)
    return pd.DataFrame(rows)


def score(df: pd.DataFrame, const_tab: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for i, tier in enumerate(TIERS):
        ok = const_tab[const_tab.tier.isin(TIERS[: i + 1])]
        full = {r.inchikey: dict(tm_K=r.tm_K, dh_J_mol=r.dh_J_mol) for r in ok.itertuples()}
        pred = []
        for r in df.itertuples():
            x, iks = json.loads(r.x_mole), json.loads(r.inchikeys)
            if any(k not in full for k in iks):
                pred.append(np.nan)
                continue
            ts = []
            for xi, k in zip(x, iks):
                if xi <= 0:
                    continue
                inv = 1.0 / full[k]["tm_K"] - (v4.R_GAS / full[k]["dh_J_mol"]) * math.log(xi)
                ts.append(1.0 / inv if inv > 0 else np.nan)
            pred.append(np.nan if (not ts or any(np.isnan(ts))) else max(ts))
        df[f"ideal_{tier}_K"] = pred
    return df


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    s = df[df.scored]
    rows = []
    for tier in TIERS:
        col = f"ideal_{tier}_K"
        for cls, g in [("ALL", s)] + list(s.groupby("system_class")):
            g = g[g[col].notna()]
            if g.empty:
                continue
            resid = (g[col] - g.T_obs_K).to_numpy()
            grp = g.system_id.to_numpy()
            n = g.system_id.nunique()
            sc = group_scores(resid, grp, "under")
            rows.append(dict(tier=tier, system_class=cls, n_systems=n, n_points=len(g), bias_K=resid.mean(), mae_K=np.abs(resid).mean(),
                             worst_under_K=float(sc.max()), conformal_rank_k=conformal_rank(n, ALPHA) if cls == "ALL" else np.nan,
                             margin_under_K=one_sided_margin(resid, grp, ALPHA, "under") if cls == "ALL" else np.nan))
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = v4.build_rows()
    raw = pd.read_csv(v4.SLE)
    comps: dict[str, tuple[str, str]] = {}
    for r in raw.itertuples():
        for i, n in zip(json.loads(r.inchis), json.loads(r.names)):
            comps[InchiToInchiKey(i)] = (n, i)
    ct = extended_constants(comps)
    ct.to_csv(OUT / "constants_extended.csv", index=False)
    df = score(df, ct)
    s = df[df.scored]
    n_v4 = s[s.ideal_svl_K.notna()].system_id.nunique()
    n_re = s[s.ideal_v4_verified_K.notna()].system_id.nunique()
    same = bool(np.allclose(s.ideal_svl_K.fillna(-1), s.ideal_v4_verified_K.fillna(-1), atol=1e-6))
    df.to_csv(OUT / "sle_eval_extended.csv", index=False)
    summ = summarise(df)
    summ.to_csv(OUT / "summary_by_tier.csv", index=False)

    blocked = s[s[f"ideal_{TIERS[-1]}_K"].isna()]
    reason = dict(zip(ct.inchikey, ct.reason))
    name = dict(zip(ct.inchikey, ct.name))
    why = {}
    for _, g in blocked.groupby("system_id"):
        bad = [f"{name[k]} ({reason[k]})" for k in json.loads(g.inchikeys.iloc[0]) if reason.get(k)]
        why[g.names.iloc[0]] = (len(g), "; ".join(bad))
    a = summ[summ.system_class == "ALL"].set_index("tier")
    L = ["# results/v5/SLE2 - SLE systems recovered by auditing fusion constants (v5 data audit)", "",
         "Generated by `scripts/v5/eval_sle_extended_v5.py`. Model, data file and margin rule are those of v4 Task 13; only the",
         "constants table is wider. Tier rules are in the script docstring and were fixed before scoring. Nothing is fitted and",
         "nothing is written to `configs/`.", "",
         f"Check: the `v4_verified` tier reproduces v4 ({n_re} systems here, {n_v4} in v4; predictions identical: {same}).", "",
         "| tier (cumulative) | systems | points | bias K (pred - obs) | MAE K | worst under-prediction K | rank k | 90 % margin K |",
         "|---|---|---|---|---|---|---|---|"]
    for t in TIERS:
        if t in a.index:
            r = a.loc[t]
            mg = "inf (too few systems)" if math.isinf(r.margin_under_K) else f"{r.margin_under_K:.2f}"
            L.append(f"| {t} | {int(r.n_systems)} | {int(r.n_points)} | {r.bias_K:+.2f} | {r.mae_K:.2f} | {r.worst_under_K:.2f} | {int(r.conformal_rank_k)} | {mg} |")
    L += ["", "By system class, widest tier:", "", "| class | systems | points | bias K | MAE K | worst under K |", "|---|---|---|---|---|---|"]
    for r in summ[(summ.tier == TIERS[-1]) & (summ.system_class != "ALL")].itertuples():
        L.append(f"| {r.system_class} | {r.n_systems} | {r.n_points} | {r.bias_K:+.2f} | {r.mae_K:.2f} | {r.worst_under_K:.2f} |")
    q, nq = ct[ct.tier != ""], ct[ct.tier == ""]
    L += ["", f"Constants: {len(q)} of {len(ct)} components qualify ({', '.join(f'{k} {v}' for k, v in q.tier.value_counts().items())}); "
              f"not qualified: {', '.join(f'{k} {v}' for k, v in nq.reason.value_counts().items())}.", "",
          f"Still blocked: {len(why)} systems, {sum(v[0] for v in why.values())} points:", ""]
    for k, (n, b) in sorted(why.items(), key=lambda kv: -kv[1][0]):
        L.append(f"* {k} ({n} points): {b}")
    L += ["", "## Limits", "",
          "* `two_compilations` is agreement between the CRC table and the Broadhurst table; they may trace to one measurement.",
          "* `single_compilation` dHfus values have no second source. Treat that tier as a diagnostic, not as verified input.",
          "* The model has no solid-solid transition term and no solid-solution term. Long n-alkanes (C20 and above) crystallise",
          "  below their rotator transition and several aromatic pairs form solid solutions, so the class table matters more than",
          "  the pooled row.",
          "* These are lab binaries and ternaries, mostly far from jet-fuel composition; a margin from them does not by itself",
          "  certify a freezing-point design (see results/v5/D4).", ""]
    (OUT / "README.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[6:]))


if __name__ == "__main__":
    main()
