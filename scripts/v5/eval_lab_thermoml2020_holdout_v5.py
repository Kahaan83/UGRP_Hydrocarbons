"""v5 S1b: one-time held-out test of raw Refutas and its lab margins on the 2020 ThermoML systems.

Readings are declared in results/v5/S1b/README.md before scoring.  Nothing is fitted; margins and bias offsets come from S1 dev
systems only.  Writes ONLY to results/v5/S1b/ (S1, configs and splits are read, never written).

Run:  PYTHONPATH=src .venv/Scripts/python.exe scripts/v5/eval_lab_thermoml2020_holdout_v5.py
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

OUT = ROOT / "results" / "v5" / "S1b"
S1_DIR = ROOT / "results" / "v5" / "S1"
S1_SCRIPT = ROOT / "scripts" / "v5" / "eval_lab_known_composition_v5.py"
NEW_TABLE = ROOT / "data" / "real_fuels" / "v5" / "thermoml2020_lab.csv"
V4_YAML = ROOT / "configs" / "uncertainty_v4.yaml"
COMPOUNDS = ROOT / "data" / "unified" / "compounds.csv"
ALPHA = 0.10
TARGET = 0.90
TAIL = 0.05
MIN_SYSTEMS = 9
REGRESSION_TOL = 1e-6
ADAMANTANE_SMARTS = "C1C2CC3CC1CC(C2)C3"
WARM_BIN = "T_ge_283K"
METHODS = ("refutas", "refutas_biascorr")
VERDICT_NA = "NOT EVALUABLE"


# --------------------------------------------------------------------------------------------- pure helpers
def binom_threshold(n: int, p: float = TARGET, tail: float = TAIL) -> int:
    """Smallest k with P(Binomial(n, p) <= k) >= tail; a margin 'holds' if covered >= k."""
    from scipy.stats import binom
    for k in range(int(n) + 1):
        if binom.cdf(k, n, p) >= tail:
            return k
    return int(n)


def coverage_verdict(covered: int, n: int) -> tuple[int | None, str]:
    """-> (k, 'HOLDS' | 'FAILS' | 'NOT EVALUABLE'); not evaluable below MIN_SYSTEMS systems."""
    if n < MIN_SYSTEMS:
        return None, VERDICT_NA
    k = binom_threshold(n)
    return k, ("HOLDS" if covered >= k else "FAILS")


def select_test_systems(ids, split: dict) -> tuple[list[str], list[str]]:
    """-> (test ids, not-independent ids): ids not in / in the S1 split (dev or locked)."""
    known = set(split["dev"]) | set(split["locked"])
    uniq = sorted(set(ids))
    return [s for s in uniq if s not in known], [s for s in uniq if s in known]


def has_adamantane(smiles) -> bool:
    """True iff the molecule contains the adamantane cage (substituted adamantanes included)."""
    if not isinstance(smiles, str):
        return False
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    return bool(mol is not None and mol.HasSubstructMatch(Chem.MolFromSmarts(ADAMANTANE_SMARTS)))


def md_table(df: pd.DataFrame) -> list[str]:
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in df.itertuples(index=False):
        out.append("| " + " | ".join("" if (isinstance(v, float) and math.isnan(v)) else
                                     (f"{v:.4f}" if isinstance(v, float) else str(v)) for v in r) + " |")
    return out


def quiet_margin(resid, groups) -> float:
    from fuelblend.uq.conformal import two_sided_margin
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return two_sided_margin(resid, groups, ALPHA)


def load_script(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def unscorable_reasons(keys, reg, provider, cls, nc, classes) -> list[str]:
    """Reason categories why a row's components cannot be scored (mirrors the checks in S1 score_rows)."""
    out = []
    for k in keys:
        i = reg.ik2i.get(k)
        if i is None:
            out.append("component not in registry")
        elif not provider.valid[i]:
            out.append("no pure-property provider data")
        elif not reg.usable[i]:
            out.append("registry marks component not usable")
        elif cls[i] not in classes:
            out.append(f"class {cls[i]} not in CLASSES")
        elif not np.isfinite(nc[i]):
            out.append("carbon number not finite")
    return sorted(set(out))


# --------------------------------------------------------------------------------------------- margins in use
def margins_in_use(s1, s1_yaml: dict, v4_yaml: dict) -> dict[tuple[str, str], float]:
    """(method, bin) -> margin in use.  Raw warm bin: the v4 config; raw cold bins and biascorr: the S1 yaml.  Absent -> inf."""
    out: dict[tuple[str, str], float] = {}
    for m in METHODS:
        for b in s1.BINS:
            key = s1.method_key(m, b)
            if m == "refutas" and b == WARM_BIN:
                ent = v4_yaml["margins"][key]
            else:
                ent = s1_yaml["margins"].get(key, {})
            out[(m, b)] = float(ent.get("applied_margin", float("inf")))
    return out


def regression_table(s1, dev_res: pd.DataFrame, s1_yaml: dict) -> pd.DataFrame:
    """R-S1b-1: recomputed dev margins vs the S1 yaml (raw cold bins; biascorr all bins)."""
    rows = []
    for m in METHODS:
        for b in s1.BINS:
            if m == "refutas" and b == WARM_BIN:
                continue
            g = dev_res[(dev_res.method == m) & (dev_res.t_bin == b) & (dev_res.split_role == "dev")]
            mine = quiet_margin(g.resid.to_numpy(), g.system_id.to_numpy())
            ref = float(s1_yaml["margins"].get(s1.method_key(m, b), {}).get("applied_margin", float("inf")))
            same = (math.isinf(mine) and math.isinf(ref)) or abs(mine - ref) < REGRESSION_TOL
            rows.append(dict(method=m, t_bin=b, dev_systems=g.system_id.nunique(), recomputed=mine, s1_yaml=ref,
                             abs_diff=(0.0 if math.isinf(mine) and math.isinf(ref) else abs(mine - ref)), match=bool(same)))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------------- coverage etc.
def coverage_rows(test_res: pd.DataFrame, margins: dict, adamantane_ids: set, bins) -> pd.DataFrame:
    from fuelblend.uq.conformal import coverage
    rows = []
    for variant, drop in (("all systems", set()), ("without adamantane systems", adamantane_ids)):
        for m in METHODS:
            for b in bins:
                g = test_res[(test_res.method == m) & (test_res.t_bin == b) & ~test_res.system_id.isin(drop)]
                mar = margins[(m, b)]
                cov, n = coverage(g.resid.to_numpy(), g.system_id.to_numpy(), mar, two_sided=True)
                k, verdict = coverage_verdict(cov, n)
                if math.isinf(mar):
                    verdict = f"{VERDICT_NA} (margin inf)"
                rows.append(dict(reading=("R-S1b-3" if m.endswith("biascorr") else "R-S1b-2") + ("" if not drop else " sensitivity (R-S1b-5)"),
                                 variant=variant, method=m, t_bin=b, margin=mar, n_systems=n, covered=cov, k=k, verdict=verdict))
    return pd.DataFrame(rows)


def accuracy_rows(groups: dict[str, pd.DataFrame], s1_summary: pd.DataFrame, bins) -> pd.DataFrame:
    rows = []
    for name, df in groups.items():
        for m in METHODS:
            for b in bins:
                g = df[(df.method == m) & (df.t_bin == b)]
                if len(g):
                    rows.append(dict(group=name, method=m, t_bin=b, n_systems=g.system_id.nunique(), n_rows=len(g),
                                     bias_ln=g.resid.mean(), rmse_ln=float(np.sqrt((g.resid ** 2).mean()))))
    for role in ("dev", "locked"):
        s = s1_summary[(s1_summary.split_role == role) & s1_summary.method.isin(METHODS) & s1_summary.t_bin.isin(bins)]
        for r in s.itertuples():
            rows.append(dict(group=f"S1 {role}", method=r.method, t_bin=r.t_bin, n_systems=r.n_systems, n_rows=r.n_records,
                             bias_ln=r.bias_ln, rmse_ln=r.rmse_ln))
    return pd.DataFrame(rows)


def per_system_rows(res: pd.DataFrame, group: str, names: dict, adam: set, margins: dict) -> list[dict]:
    rows = []
    raw = res[res.method == "refutas"]
    for (sid, b), g in raw.groupby(["system_id", "t_bin"]):
        mx = float(g.resid.abs().max())
        cor = res[(res.method == "refutas_biascorr") & (res.system_id == sid) & (res.t_bin == b)]
        cmx = float(cor.resid.abs().max()) if len(cor) else float("nan")
        rows.append(dict(group=group, system_id=sid, components=names[sid], adamantane=sid in adam, t_bin=b, n_rows=len(g),
                         mean_ln_resid=float(g.resid.mean()), max_abs_ln_resid=mx, margin_raw=margins[("refutas", b)],
                         covered_raw=("yes" if mx <= margins[("refutas", b)] else "no"),
                         max_abs_ln_resid_biascorr=cmx,
                         covered_biascorr=("" if not len(cor) else ("yes" if cmx <= margins[("refutas_biascorr", b)] else "no")),
                         status="scored", reason=""))
    return rows


# ------------------------------------------------------------------------------------------------- main
def main() -> None:
    import yaml
    OUT.mkdir(parents=True, exist_ok=True)
    s1 = load_script("s1v5", S1_SCRIPT)
    bins = s1.BINS
    split = json.loads((S1_DIR / "lab_systems_splits_v5.json").read_text())
    s1_yaml = yaml.safe_load((S1_DIR / "uncertainty_v5.yaml").read_text(encoding="utf-8"))
    v4_yaml = yaml.safe_load(V4_YAML.read_text(encoding="utf-8"))
    s1_summary = pd.read_csv(S1_DIR / "summary_by_bin_method.csv")
    margins = margins_in_use(s1, s1_yaml, v4_yaml)

    meas4 = pd.read_csv(s1.MEAS, low_memory=False)
    boehm = pd.read_csv(s1.BOEHM, low_memory=False)
    boehm["sources"] = boehm["sources"].fillna("boehm2022")
    new = pd.read_csv(NEW_TABLE, low_memory=False)
    old = pd.concat([meas4, boehm], ignore_index=True)
    old_dev = old[old.system_id.isin(set(split["dev"]))]

    lb = s1._load("build_viscosity_leaderboard_v4", "scripts/v4/build_viscosity_leaderboard_v4.py")
    from fuelblend.mixture.cold_flow import compound_meta
    from fuelblend.mixture.common import CLASSES, CompoundRegistry, PurePropertyProvider
    reg = CompoundRegistry(COMPOUNDS)
    _, cls, nc = compound_meta(COMPOUNDS)
    provider = PurePropertyProvider(reg, ROOT / "data" / "cache" / "pure_grid.npz", verbose=False)
    smiles = pd.read_csv(COMPOUNDS).set_index("inchikey")["smiles"].to_dict()

    # ---- R-S1b-1 first: S1 dev systems only, new rows not scored yet
    dev_res = s1.score_rows(old_dev, reg, provider, cls, nc, lb, CLASSES)
    dev_res["split_role"] = "dev"
    dev_res = s1.add_bias_corrected(dev_res)
    reg_tab = regression_table(s1, dev_res, s1_yaml)
    print("R-S1b-1 regression:\n" + reg_tab.to_string(index=False), flush=True)
    reg_ok = bool(reg_tab.match.all())
    if not reg_ok:
        (OUT / "RESULTS.md").write_text("# v5 S1b results\n\n## R-S1b-1 FAILED: margins recomputed from S1 dev systems do not equal the S1 yaml. "
                                        "Per the README nothing else was read.\n\n" + "\n".join(md_table(reg_tab)) + "\n", encoding="utf-8")
        raise SystemExit("R-S1b-1 failed; stopping")

    # ---- test-system selection (ids only)
    sel = s1.lab_selection(new)
    test_ids, ni_ids = select_test_systems(sel.system_id.unique(), split)
    new_res = s1.score_rows(new[new.system_id.isin(set(test_ids) | set(ni_ids))], reg, provider, cls, nc, lb, CLASSES)
    scored_ids = set(new_res.system_id)
    ns_ids = [s for s in test_ids + ni_ids if s not in scored_ids]
    first = sel.drop_duplicates("system_id").set_index("system_id")
    names = {s: " + ".join(smiles.get(k, k[:14]) for k in json.loads(first.loc[s, "comp_inchikeys"])) for s in first.index}
    adam = {s for s in first.index if any(has_adamantane(smiles.get(k)) for k in json.loads(first.loc[s, "comp_inchikeys"]))}
    (OUT / "test_systems.json").write_text(json.dumps(dict(
        test=test_ids, not_independent=ni_ids, not_scorable=ns_ids, not_scorable_among_test=[s for s in ns_ids if s in test_ids],
        adamantane_systems=sorted(adam & set(test_ids))), indent=1) + "\n")

    # ---- bias-corrected rows: test systems get role 'locked' only for this call
    new_res["split_role"] = new_res.system_id.map(lambda s: "locked" if s in set(test_ids) else "not_independent")
    comb = s1.add_bias_corrected(pd.concat([dev_res[dev_res.method.isin(s1.METHODS)].drop(columns=["resid_raw", "offset"]),
                                            new_res], ignore_index=True))
    test_res = comb[comb.split_role == "locked"].assign(group="test")
    ni_res = comb[comb.split_role == "not_independent"].assign(group="not_independent")
    pd.concat([test_res, ni_res], ignore_index=True).drop(columns=["split_role"]).to_csv(
        OUT / "viscosity_test_residuals.csv", index=False, float_format="%.8g")

    # ---- readings
    cov_df = coverage_rows(test_res, margins, adam, bins)
    cov_df.to_csv(OUT / "coverage.csv", index=False, float_format="%.8g")
    acc_df = accuracy_rows({"test": test_res, "not_independent (raw only)": ni_res[ni_res.method == "refutas"]}, s1_summary, bins)
    acc_df.to_csv(OUT / "accuracy.csv", index=False, float_format="%.6g")

    ps = per_system_rows(test_res, "test", names, adam, margins) + per_system_rows(ni_res, "not_independent", names, adam, margins)
    sel_n = sel.groupby("system_id").size()
    for s in ns_ids:
        ks = json.loads(first.loc[s, "comp_inchikeys"])
        reasons = unscorable_reasons(ks, reg, provider, cls, nc, CLASSES) or ["composition or prediction invalid"]
        ps.append(dict(group="test" if s in test_ids else "not_independent", system_id=s, components=names[s], adamantane=s in adam,
                       t_bin="", n_rows=int(sel_n[s]), status="not scorable", reason="; ".join(reasons)))
    ps_df = pd.DataFrame(ps)
    ps_df.to_csv(OUT / "per_system.csv", index=False, float_format="%.6g")

    write_results(OUT / "RESULTS.md", reg_tab, cov_df, acc_df, ps_df, margins, test_ids, ni_ids, ns_ids, adam, sel, new_res, bins)
    print((OUT / "RESULTS.md").read_text(encoding="utf-8"), flush=True)


def write_results(path, reg_tab, cov_df, acc_df, ps_df, margins, test_ids, ni_ids, ns_ids, adam, sel, new_res, bins) -> None:
    ns_test = [s for s in ns_ids if s in test_ids]
    scored_test = [s for s in test_ids if s not in ns_ids]
    dropped = int(len(sel[sel.system_id.isin(test_ids + ni_ids)]) - new_res[new_res.method == "refutas"].shape[0])
    L = ["# v5 S1b results (generated by scripts/v5/eval_lab_thermoml2020_holdout_v5.py)", "",
         "Readings declared in README.md before scoring. Scored once. These systems now count as read (R-S1b-7).", "",
         f"Test systems (passing lab_selection, not in the S1 split): {len(test_ids)}; scorable {len(scored_test)}; not scorable {len(ns_test)}. "
         f"Not independent (already in the S1 split): {len(ni_ids)} (not scorable among them: {len(ns_ids) - len(ns_test)}). "
         f"Selected rows not scored: {dropped}. Adamantane-containing scored test systems: {len(adam & set(scored_test))}.", "",
         "## R-S1b-1 regression (margins recomputed from S1 dev systems vs uncertainty_v5.yaml)", ""]
    L += md_table(reg_tab.assign(recomputed=reg_tab.recomputed.map(lambda v: f"{v:.9f}"), s1_yaml=reg_tab.s1_yaml.map(lambda v: f"{v:.9f}"),
                                 abs_diff=reg_tab.abs_diff.map(lambda v: f"{v:.2e}")))
    L += ["", f"R-S1b-1: {'HOLDS' if reg_tab.match.all() else 'FAILS'} ({int(reg_tab.match.sum())}/{len(reg_tab)} margins match to < 1e-6).", ""]
    L += ["## Margins in use", ""] + md_table(pd.DataFrame([dict(method=m, t_bin=b, margin=v) for (m, b), v in margins.items()])) + [""]
    for reading, title in (("R-S1b-2", "R-S1b-2 raw Refutas coverage (primary)"), ("R-S1b-3", "R-S1b-3 bias-corrected Refutas coverage"),
                           ("R-S1b-2 sensitivity (R-S1b-5)", "R-S1b-5 raw Refutas coverage without adamantane systems"),
                           ("R-S1b-3 sensitivity (R-S1b-5)", "R-S1b-5 bias-corrected coverage without adamantane systems")):
        g = cov_df[cov_df.reading == reading].drop(columns=["reading", "variant"])
        L += [f"## {title}", ""] + md_table(g) + [""]
        for r in g.itertuples():
            L.append(f"* {reading} {r.method} {r.t_bin}: n={r.n_systems}, covered={r.covered}, k={'n/a' if r.k is None or pd.isna(r.k) else int(r.k)}, "
                     f"margin={r.margin:.4f} -> {r.verdict}")
        L.append("")
    L += ["## R-S1b-4 accuracy (ln(pred/meas); S1 dev and locked from results/v5/S1/summary_by_bin_method.csv)", ""] + md_table(acc_df) + [""]
    L += ["## R-S1b-6 per-system table (test group, then not independent, then not scorable)", ""]
    cols = ["group", "system_id", "components", "adamantane", "t_bin", "n_rows", "mean_ln_resid", "max_abs_ln_resid", "covered_raw",
            "max_abs_ln_resid_biascorr", "covered_biascorr", "status", "reason"]
    L += md_table(ps_df.reindex(columns=cols)) + [""]
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
