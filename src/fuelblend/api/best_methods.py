"""Best method per property under the v5 to v7 evidence (dev rows, so errors are optimistic).

Two scopes:
  pure   single-compound properties (Track B pool, scaffold 5-fold CV). Method = lowest committed OOF MAE among the
         deployable fixed configurations: RDKit-descriptor random forest, or the fixed long-cap D-MPNN (v1 config 0).
  blend  blend properties. v7 adopted nothing: every v7 candidate was tied, failed its gate or had a CI that includes 0,
         so the designated route stays the v5 route in `predict_route.predict_fuel`.

Each pure entry names the committed CSV and column its MAE comes from; tests/test_v7_best_methods.py reads them back.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Optional

TB = "results/v7/track_b"
RF_SRC = (f"{TB}/descriptor_gbm/cv_summary.csv", "mae_rf")
DM_SRC = (f"{TB}/curves_v1_long/cv_summary.csv", "mae_long")
CONST_SRC = (f"{TB}/curves_v1_long/cv_summary.csv", "mae_const")


@dataclass(frozen=True)
class PureEntry:
    prop: str                  # column in the Track B pool
    unit: str
    method: str                # "rf_rdkit2d" or "dmpnn_long"
    src_csv: str
    src_col: str
    n: int
    significance: str          # reading of the committed paired bootstrap, never stronger
    note: str = ""


PURE = {
    "flash_C": PureEntry("flash_C", "degC", "rf_rdkit2d", *RF_SRC, 837,
                         "significant: beats the long-cap D-MPNN (Holm p 0.0135) and LightGBM"),
    "cetane": PureEntry("cetane", "-", "rf_rdkit2d", *RF_SRC, 259,
                        "tied with the long-cap D-MPNN (8.66); CI of the gain includes 0"),
    "dcn": PureEntry("dcn", "-", "dmpnn_long", *DM_SRC, 89,
                     "nominal: beats the random forest by 1.76 (90% CI 0.68 to 2.66) but Holm p 0.12, not below 0.10"),
    "density": PureEntry("density", "kg/m3", "dmpnn_long", *DM_SRC, 558,
                         "not significant: CI of the gain over LightGBM includes 0"),
    "dhfus_kJ_mol": PureEntry("dhfus_kJ_mol", "kJ/mol", "rf_rdkit2d", *RF_SRC, 81,
                              "significant: beats the long-cap D-MPNN (Holm p 0.005)"),
    "logvisc_mPas_293": PureEntry("logvisc_mPas_293", "log10 mPa.s at 293 K", "rf_rdkit2d", *RF_SRC, 180,
                                  "significant: beats the long-cap D-MPNN (Holm p 0.028)"),
    "mon": PureEntry("mon", "-", "dmpnn_long", *DM_SRC, 363,
                     "nominal: gain over LightGBM 1.27 (CI 0.19 to 1.90) but Holm p 0.36"),
    "ron": PureEntry("ron", "-", "dmpnn_long", *DM_SRC, 378,
                     "nominal: gain over LightGBM 1.50 (CI 0.21 to 2.36) but Holm p 0.35"),
    "tm_K": PureEntry("tm_K", "K", "rf_rdkit2d", *RF_SRC, 547,
                      "tied: lowest MAE of the fixed configurations; the long-cap D-MPNN is 28.01, CI of the gap includes 0"),
    "ysi": PureEntry("ysi", "YSI", "rf_rdkit2d", *RF_SRC, 168,
                     "tied with the long-cap D-MPNN (73.18); simpler model kept. Locked ysi_holdout rows excluded from training"),
}

# Blend properties: designated route stays v5.  `v7_candidate` is the best v7 alternative and why it was not adopted.
BLEND = {
    "viscosity": dict(method="v5 Refutas + optional LightGBM pair residual", status="adopted_v5",
                      v7_candidate="set transformer", why_not="fails the predeclared gate against the pair baseline"),
    "density": dict(method="v5 ideal volume mixing + optional LightGBM residual", status="adopted_v5",
                    v7_candidate="deep sets", why_not="CI of the gain over v5 pair includes 0"),
    "flash_point": dict(method="v5 Le Chatelier + Raoult + Antoine + LFL overlay (designated)", status="adopted_v5",
                        v7_candidate="HANNA gamma, UNIFAC gamma, Antoine from parameter heads",
                        why_not="HANNA no skill against ideal gamma; parameter-head Antoine over-predicts by 16 degC; hydrocarbon binaries: ideal gamma lowest"),
    "flash_point_oxygenates": dict(method="v6 MolPool (exploratory)", status="exploratory",
                                   v7_candidate="v7 MolPool", why_not="v6 arm failed its primary gate on the locked set; v7 no better"),
    "freezing_point": dict(method="v5 designated (Schroeder-van Laar + UNIFAC, residual per decision log)", status="adopted_v5",
                           v7_candidate="joint pure+fuel MolPool LGBM",
                           why_not="tied with fuel-only and v6 residual; exploratory dev-only; real-fuel freezing stays: status: no skill shown"),
    "octane_mixtures": dict(method="v6 MolPool heads (research only)", status="exploratory",
                            v7_candidate="v7 MolPool", why_not="worse than v6 on Leenhouts MON"),
}

PURE_METHOD_LABEL = {"rf_rdkit2d": "random forest on RDKit 2D descriptors (500 trees, median imputation)",
                     "dmpnn_long": "D-MPNN, v1 config 0 (hidden 200, depth 4, dropout 0.1), 1000-epoch cap, patience 100, 5-model ensemble, multi-task"}


def best_method(prop: str, scope: str = "pure") -> dict:
    """Registry entry for a property. scope 'pure' or 'blend'."""
    if scope == "pure":
        if prop not in PURE:
            raise KeyError(f"no pure-compound entry for {prop!r}; known: {sorted(PURE)}")
        d = asdict(PURE[prop])
        d["method_label"] = PURE_METHOD_LABEL[d["method"]]
        return d
    if scope == "blend":
        if prop not in BLEND:
            raise KeyError(f"no blend entry for {prop!r}; known: {sorted(BLEND)}")
        return dict(BLEND[prop])
    raise ValueError("scope must be 'pure' or 'blend'")


def evidence_mae(prop: str, root: Optional[str] = None) -> float:
    """Read the registry MAE of a pure property back from its committed CSV."""
    import pandas as pd
    from pathlib import Path
    e = PURE[prop]
    base = Path(root) if root else Path(__file__).resolve().parents[3]
    df = pd.read_csv(base / e.src_csv).set_index("target")
    return float(df.loc[prop, e.src_col])
