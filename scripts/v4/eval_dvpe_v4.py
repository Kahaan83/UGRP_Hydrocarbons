"""Task 15 fix: DVPE diagnostic on Co-Optima "N% <oxygenate> in 2016-NREL BOB-01" rows -> results/v4/O/dvpe_eval_v4.csv.

UNITS-UNCONFIRMED: the raw column Blend_DVPE has an EMPTY units column.  Evidence for the unit of the NREL NB6022
rows (see `unit_evidence()` and the README): other rows of the same column are kPa (Vertifuel 35.4 with source text
"RVP 5.14 psi"; PNNL iso-olefin 17.86 under "ASTM D5191, EPA RVP"), but the NREL BOB rows are 4.2-8.8, impossible
as kPa for a gasoline blend (neat gasolines are ~35-90 kPa) and consistent with psi (E10 = 6.44 -> ethanol +~1 psi).
So the NREL rows are treated as psi (INFERRED, not documented); metrics are reported in the file's native numbers,
and the Raoult baseline converts oxygenate vapour pressures to psi.  Treat every number here as a diagnostic.

Rows excluded from exact targets (all counted in the CSV): not a BOB+oxygenate row; censored ("<x"); DVPE temperature
not 37.8 C (two methylfuran rows at 25 C); named mixtures (ketone mix, furan mix, diisobutylene isomer mix).
No neat-BOB (0 % oxygenate) DVPE row exists: the BOB constant is FITTED in-fold.
Groups = oxygenate CAS; dev/locked use the same CAS groups and folds as the bRON octane split (R19).
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import LeaveOneGroupOut

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.models import octane_blend_v4 as O  # noqa: E402

RAW = ROOT / "extras" / "fuelblend_data" / "data" / "raw" / "cooptima" / "fuel_engine_co_optimization.csv"
SPLITS = ROOT / "splits" / "v4" / "real_fuels_v4_splits.json"
OUT = ROOT / "results" / "v4" / "O" / "dvpe_eval_v4.csv"
PSI_TO_KPA = 6.894757
T_DVPE_K = 310.95                      # 37.8 C
NAME_RE = re.compile(r"^\s*(\d+)%\s+(.+?)\s+in\s+2016-NREL[ -]BOB-01", re.I)
MIX_RE = re.compile(r"mix", re.I)


def fix_cas(s):
    """Normalise typed CAS numbers ('137.32-6', '141-78*-6', '71.36-3') by rebuilding from the digits."""
    d = re.sub(r"\D", "", str(s))
    return f"{d[:-3]}-{d[-3:-1]}-{d[-1]}" if len(d) >= 5 else None


def parse_rows(df):
    """Return (resolvable rows table, exclusion counts)."""
    ex = dict(n_rows_with_dvpe=int(df.Blend_DVPE.notna().sum()))
    x = df[df.Blend_DVPE.notna()]
    rows, c = [], dict(not_nrel_bob_oxygenate=0, censored=0, temp_not_37p8=0, named_mixture=0)
    for r in x.to_dict("records"):
        raw = str(r["Blend_DVPE"]).strip()
        if raw.startswith("<") or raw.startswith(">"):
            c["censored"] += 1
            continue
        m = NAME_RE.match(str(r["Blend_Name"]))
        if not m:
            c["not_nrel_bob_oxygenate"] += 1
            continue
        if MIX_RE.search(str(r["Blend_Name"])) or MIX_RE.search(str(r["Blend_Primary Component"])):
            c["named_mixture"] += 1
            continue
        if not np.isclose(pd.to_numeric(r["Blend_DVPE_temp"], errors="coerce"), 37.8, atol=0.3):
            c["temp_not_37p8"] += 1
            continue
        cas = fix_cas(r["Blend component primary CAS"])
        wt = re.match(r"\s*([\d.]+)\s*wt%", str(r["Blend_Primary Component"]))
        rows.append(dict(name=r["Blend_Name"], cas=cas, group=f"BRON:{cas}", f_vol=int(m.group(1)) / 100.0,
                         wt_pct=float(wt.group(1)) if wt else np.nan, y=float(raw),
                         source=r["Blend_DVPE_source"]))
    ex.update({f"excluded_{k}": v for k, v in c.items()})
    ex["n_resolvable"] = len(rows)
    return pd.DataFrame(rows), ex


def add_features(t):
    from thermo import Chemical
    ds, ps = [], []
    for cas in t.cas:
        d = O.blendstock_descriptors(cas, "")
        ds.append(d)
        ps.append(Chemical(cas, T=T_DVPE_K).Psat / (PSI_TO_KPA * 1000.0))   # Pa -> psi
    t = t.copy()
    t["mw"] = [d["mw"] for d in ds]
    t["rho"] = [d["rho"] for d in ds]
    t["p_ox_psi"] = ps
    nb = t.f_vol * t.rho / t.mw
    na = (1 - t.f_vol) * O.BASE_RHO_ASSUMED / O.BASE_MW_ASSUMED
    t["x_mole"] = nb / (nb + na)
    return t


# ----------------------------------------------------------------------------- models (fit on train only)
class DConst:
    name = "constant_mean"

    def fit(self, tr):
        self.m = tr.y.mean()
        return self

    def predict(self, d):
        return np.full(len(d), self.m)


class DLinearVol:
    """P = P_BOB + b*f_vol with ONE pooled slope b (least squares on training rows)."""
    name = "linear_pooled_vol"

    def fit(self, tr):
        A = np.column_stack([np.ones(len(tr)), tr.f_vol])
        self.c = np.linalg.lstsq(A, tr.y.to_numpy(), rcond=None)[0]
        return self

    def predict(self, d):
        return self.c[0] + self.c[1] * d.f_vol.to_numpy()


class DRaoult:
    """Ideal mixing of vapour pressures by mole fraction: P=(1-x)P_BOB + x P_ox(37.8C); P_BOB fitted in-fold."""
    name = "raoult_ideal_mole"

    def fit(self, tr):
        w = 1 - tr.x_mole.to_numpy()
        self.pb = float((w * (tr.y.to_numpy() - tr.x_mole.to_numpy() * tr.p_ox_psi.to_numpy())).sum() / (w ** 2).sum())
        return self

    def predict(self, d):
        x = d.x_mole.to_numpy()
        return (1 - x) * self.pb + x * d.p_ox_psi.to_numpy()


MODELS = [DConst, DLinearVol, DRaoult]


def evaluate(view, tab, folds):
    P = {m.name: {} for m in (c() for c in MODELS)}
    for tr, te in folds:
        for cls in MODELS:
            m = cls().fit(tab.iloc[tr])        # fresh model per fold
            for i, p in zip(te, m.predict(tab.iloc[te].drop(columns=["y"]))):
                P[m.name][i] = p
    idx = sorted(next(iter(P.values())))
    y, g = tab.y.to_numpy()[idx], tab.group.to_numpy()[idx]
    err = {k: np.abs(np.array([v[i] for i in idx]) - y) for k, v in P.items()}
    out = []
    for k, v in P.items():
        p = np.array([v[i] for i in idx])
        mt = O.metrics(y, p)
        _, lo, hi, ng = O.cluster_boot(g, err[k])
        out.append(dict(record_type="metric", view=view, model=k, ref_model="", n=mt["n"], n_groups=ng, mae=mt["mae"],
                        rmse=mt["rmse"], bias=mt["bias"], ci_lo=lo, ci_hi=hi, excludes_zero="",
                        units="UNITS-UNCONFIRMED (inferred psi)", note=""))
        if k != "constant_mean":
            d, dlo, dhi, ngd = O.cluster_boot(g, err[k], err["constant_mean"])
            out.append(dict(record_type="paired_diff", view=view, model=k, ref_model="constant_mean", n=mt["n"],
                            n_groups=ngd, mae=d, ci_lo=dlo, ci_hi=dhi, excludes_zero=bool(dlo > 0 or dhi < 0),
                            units="UNITS-UNCONFIRMED (inferred psi)",
                            note="mae = MAE(model)-MAE(ref); negative = model better"
                                 + ("; locked has 4 groups, CI NOT reliable" if view.startswith("locked") else "")))
    return out


def main():
    df = pd.read_csv(RAW, low_memory=False, encoding="latin-1")
    t, ex = parse_rows(df)
    t = add_features(t)
    sp = json.loads(SPLITS.read_text())
    locked_groups = set(sp["locked_groups"])
    fold_of_group = {}
    for blend, grp in sp["fuel_group"].items():
        for k, v in sp["folds"].items():
            if blend in v:
                fold_of_group[grp] = k
    t["is_locked"] = t.group.isin(locked_groups)
    t["fold"] = t.group.map(fold_of_group)
    unmapped = sorted(set(t.group[~t.group.isin(set(sp["fuel_group"].values()))]))
    ex["groups_total"] = int(t.group.nunique())
    ex["groups_unmapped_to_bron_split"] = unmapped
    dev = t[~t.is_locked & t.fold.notna()].reset_index(drop=True)
    lock = t[t.is_locked].reset_index(drop=True)
    ex.update(n_dev=len(dev), n_locked=len(lock), groups_dev=int(dev.group.nunique()), groups_locked=int(lock.group.nunique()),
              n_dev_unmapped_dropped=int((~t.is_locked & t.fold.isna()).sum()))
    out = []
    out += evaluate("dev_cv5_stored_folds", dev, [(np.where(dev.fold != k)[0], np.where(dev.fold == k)[0]) for k in sorted(dev.fold.unique())])
    out += evaluate("dev_leave_one_oxygenate_out", dev, list(LeaveOneGroupOut().split(dev, groups=dev.group)))
    both = pd.concat([dev, lock], ignore_index=True)
    out += evaluate("locked_thin_test", both, [(np.arange(len(dev)), np.arange(len(dev), len(both)))])
    res = pd.DataFrame(out)
    # locked rows: evaluate() computes over `both` test indices only (test = locked), so n/groups are locked counts
    OUT.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(OUT, index=False, float_format="%.6g")
    (OUT.parent / "dvpe_prep_report.json").write_text(json.dumps(ex, indent=1, default=str))
    print("wrote", OUT, len(res), ex)


def section(df, ex):
    f = lambda v, d=2: "nan" if pd.isna(v) else f"{v:.{d}f}"
    L = ["", "## DVPE diagnostic (UNITS-UNCONFIRMED)", "",
         "Source: Co-Optima `fuel_engine_co_optimization.csv`, `Blend_DVPE` on rows `N% <oxygenate> in 2016-NREL BOB-01` "
         "(source text `NREL NB6022 pg 016-017`). The units column is empty; DVPE temperature is filled (37.8 C on the rows used). "
         "Evidence: Vertifuel (35.4, source text `RVP 5.14 psi`) shows the column is kPa for some rows, but the NREL rows (4.2-8.8) "
         "are impossible as kPa for gasoline and match psi (E10 = 6.44). Unit for the NREL rows is therefore INFERRED to be psi, "
         "not documented; all numbers below are in the file's native numbers (psi if the inference holds) and are a diagnostic only.",
         "", f"* Rows with a DVPE value: {ex['n_rows_with_dvpe']}; resolvable BOB+oxygenate rows used: {ex['n_resolvable']} "
         f"({ex['groups_total']} oxygenate CAS groups; dev {ex['n_dev']} rows / {ex['groups_dev']} groups, locked {ex['n_locked']} / {ex['groups_locked']}).",
         f"* Excluded and counted: not a BOB+oxygenate row {ex['excluded_not_nrel_bob_oxygenate']}, censored `<x` {ex['excluded_censored']}, "
         f"DVPE temperature not 37.8 C {ex['excluded_temp_not_37p8']} (methylfuran at 25 C), named mixtures {ex['excluded_named_mixture']} "
         "(ketone mix, furan mix, diisobutylene isomer mix).",
         "* No neat-BOB DVPE row exists (the real gap): the BOB constant is fitted in-fold. Nominal volume % from the row name is used "
         "as the composition (wt % given for most rows is not used); the Raoult baseline assumes a BOB of M=100 g/mol, rho=745 kg/m3 (typed) "
         "and 37.8 C oxygenate vapour pressures from `thermo`. Groups = oxygenate CAS; the same CAS dev/locked/fold assignment as the octane split.",
         f"* Unmapped to the bRON split: {ex['groups_unmapped_to_bron_split'] or 'none'}.", ""]
    m = df[df.record_type == "metric"]
    for view in ("dev_cv5_stored_folds", "dev_leave_one_oxygenate_out", "locked_thin_test"):
        s = m[m.view == view]
        L += [f"### {view}" + (" (thin; point estimates only; CIs unreliable)" if view.startswith("locked") else ""), "",
              "| model | N | groups | MAE | 95% cluster CI | RMSE | bias |", "|---|---|---|---|---|---|---|"]
        for r in s.itertuples():
            ci = "n/a" if view.startswith("locked") else f"[{f(r.ci_lo)}, {f(r.ci_hi)}]"
            L.append(f"| {r.model} | {r.n} | {r.n_groups} | {f(r.mae)} | {ci} | {f(r.rmse)} | {f(r.bias)} |")
        L.append("")
    d = df[(df.record_type == "paired_diff") & (~df.view.str.startswith("locked"))]
    L += ["### Paired MAE differences vs constant_mean (negative = model better)", "",
          "| view | model | dMAE | 95% CI | CI excludes 0 |", "|---|---|---|---|---|"]
    for r in d.itertuples():
        L.append(f"| {r.view} | {r.model} | {f(r.mae)} | [{f(r.ci_lo)}, {f(r.ci_hi)}] | {r.excludes_zero} |")
    L += ["", "Only CI-supported statements should be drawn from these tables. The pooled-slope linear model collapses to the "
          "constant mean (oxygenates raise or lower DVPE with opposite signs, so one slope averages to ~0). The ideal-Raoult model, "
          "with only the BOB constant fitted, is CI-separated from the constant mean on dev; because it needs absolute oxygenate "
          "vapour pressures, that is weak circumstantial support for the psi inference but does NOT confirm the unit. "
          "Non-ideality (alcohol/hydrocarbon) is not modelled, and the locked set (4 groups) cannot rank models.", ""]
    return L


if __name__ == "__main__":
    main()
