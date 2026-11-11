"""v5 S2 Question C: flash point of T2 fuels (measured GCxGC composition + measured flash) with ideal and UNIFAC models.
Readings R-C1..R-C4 in results/v5/S2/README.md.  Nothing is tuned; outputs only under results/v5/S2/.

Order of work (the maps are written to disk BEFORE any flash prediction is made):
  1. Parse the r3 representative table (results/v4/R3/reconstruction_report.md) -> (sub-kind, carbon number) -> SMILES.
  2. Map every class x carbon-number bin of AV-33-22 and of the 4 FuelLib fuels with a measured flash point to a representative,
     with the same rule as src/fuelblend/data/reconstruction_r3.SubReps.get: nearest carbon number, max gap 3
     (max gap 99 for cyclo_di, indan_tetralin, aromatic_di).  Write av3322_representative_map.csv / fuellib_representative_map.csv.
  3. Mass fractions of mapped bins -> renormalise -> mole fractions with RDKit MW; predict flash (ideal, original UNIFAC).
  4. Score against coa_flash_point_C (AV-33-22 bulk file) and the FuelLib flash_point_C row.

Primary scoring uses the flash registry exactly as the r3 scores did (mesitylene constants as stored in
data/processed/v2/compounds_anchored_vp.csv).  A second variant substitutes the WebBook mesitylene set chosen in s2_mesitylene.py
(mesitylene_best_input.json); it is reported next to the primary and is NOT what the readings are judged on.
"""
from __future__ import annotations

import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from scipy import stats

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.baselines.flash_unifac import UnifacFallbackWarning, _default_registry, _lookup, flash_point_unifac_detail  # noqa: E402

OUT = ROOT / "results" / "v5" / "S2"
V5 = ROOT / "data" / "real_fuels" / "v5"
REPORT = ROOT / "results/v4/R3/reconstruction_report.md"
FUELLIB_FLASH = ["posf10264", "posf10289", "posf10325", "posf11498"]
MAX_UNMAPPED = 0.10
OVER = 5.0
NEAREST_ANY = ("cyclo_di", "indan_tetralin", "aromatic_di")
MES = "Cc1cc(C)cc(C)c1"

AV_GROUP = {"n-alkane": "n_paraffin", "iso-alkane": "iso_paraffin", "monocycloalkane": "cyclo_mono", "dicycloalkane": "cyclo_di",
            "alkylbenzene": "alkylbenzene", "cycloaromatic": "indan_tetralin", "diaromatic": "aromatic_di", "alkene": "alkene_alkyne",
            "tricycloalkane": "tricycloalkane_no_rep"}
AROMATIC_KINDS = {"alkylbenzene", "indan_tetralin", "aromatic_di"}


def fuellib_kind(gc_bin: str) -> str:
    b = gc_bin
    if re.match(r"n-C\d+$", b):
        return "n_paraffin"
    for key, kind in [("Isoparaffin", "iso_paraffin"), ("Monocycloparaffin", "cyclo_mono"), ("Dicycloparaffin", "cyclo_di"),
                      ("Tricycloparaffin", "tricycloalkane_no_rep"), ("Diaromatic", "aromatic_di"), ("Cycloaromatic", "indan_tetralin"),
                      ("Alkene", "alkene_alkyne")]:
        if key in b:
            return kind
    if b == "Toluene" or re.match(r"C\d+-Benzene$", b):
        return "alkylbenzene"
    raise ValueError(f"unrecognised FuelLib gc_bin {b!r}")


def parse_rep_table(path: Path = REPORT) -> dict:
    """{(sub-kind, carbon number): smiles} from the 'Representative per (sub-kind, C#) bin' table of the r3 report."""
    txt = path.read_text(encoding="utf-8")
    sec = txt.split("## Representative per", 1)[1].split("\n## ", 1)[0]
    out = {}
    for line in sec.splitlines():
        m = re.match(r"\|\s*([a-z_]+)\s*\|\s*(\d+)\s*\|\s*[^|]*\|\s*([A-Z]{14}-[A-Z]{10}-[A-Z])\s*\|\s*([^|\s]+)\s*\|", line)
        if m:
            out[(m.group(1), int(m.group(2)))] = (m.group(4), m.group(3))
    if not out:
        raise RuntimeError("no representative rows parsed")
    return out


def rep_for(table: dict, kind: str, cn: int):
    """Same rule as SubReps.get: nearest cn (ties -> lower cn), gap limit 3 (99 for the few-molecule classes). -> (smiles, ik, gap) or None."""
    ns = sorted(n for (k, n) in table if k == kind)
    if not ns:
        return None
    best = min(ns, key=lambda k: (abs(k - cn), k))
    gap = abs(best - cn)
    if gap > (99 if kind in NEAREST_ANY else 3):
        return None
    s, ik = table[(kind, best)]
    return s, ik, best, gap


def build_maps(table: dict):
    comp = pd.read_csv(V5 / "crc_av3322_composition.csv")
    comp["kind"] = comp.group.map(AV_GROUP)
    assert comp.kind.notna().all(), comp[comp.kind.isna()].group.unique()
    amap = comp.groupby(["group", "kind", "carbon_number"]).agg(mean_mass_pct=("mass_pct", "mean"), n_fuels=("fuel_id", "nunique"),
                                                               any_plus_bin=("plus", "max")).reset_index()
    fl = pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv")
    fl = fl[fl.fuel_id.isin([f"FUELLIB_{s}" for s in FUELLIB_FLASH])].copy()
    fl["kind"] = fl.gc_bin.map(fuellib_kind)
    fmap = fl.groupby(["gc_bin", "kind", "n_carbon"]).agg(mean_mass_fraction=("fraction", "mean"), n_fuels=("fuel_id", "nunique")).reset_index()

    reg = _default_registry()

    def resolve(kind, cn):
        if kind == "tricycloalkane_no_rep":
            return dict(rep_smiles="", rep_inchikey="", rep_cn=np.nan, cn_gap=np.nan, mapped=False, reason="no realistic representative (tricycloparaffin)")
        r = rep_for(table, kind, int(cn))
        if r is None:
            return dict(rep_smiles="", rep_inchikey="", rep_cn=np.nan, cn_gap=np.nan, mapped=False,
                        reason="no representative in sub-kind" if not any(k == kind for k, _ in table) else "nearest carbon number > gap limit")
        s, ik, rcn, gap = r
        can = Chem.MolToSmiles(Chem.MolFromSmiles(s))
        if _lookup(can, reg) is None:
            return dict(rep_smiles=can, rep_inchikey=ik, rep_cn=rcn, cn_gap=gap, mapped=False, reason="representative lacks Antoine or LFL in flash registry")
        return dict(rep_smiles=can, rep_inchikey=ik, rep_cn=rcn, cn_gap=gap, mapped=True, reason="")

    amap = pd.concat([amap, pd.DataFrame([resolve(k, c) for k, c in zip(amap.kind, amap.carbon_number)])], axis=1)
    fmap = pd.concat([fmap, pd.DataFrame([resolve(k, c) for k, c in zip(fmap.kind, fmap.n_carbon)])], axis=1)
    return comp, amap, fl, fmap


def mass_to_mole(smiles, w):
    mw = np.array([Descriptors.MolWt(Chem.MolFromSmiles(s)) for s in smiles])
    n = np.asarray(w, float) / mw
    return n / n.sum()


def fuel_mixture(df: pd.DataFrame, key_cols, mass_col, scale, mp: pd.DataFrame):
    """-> (smiles list, mass fractions renormalised over mapped mass, coverage dict).  df rows: bins of ONE fuel."""
    m = df.merge(mp[key_cols + ["rep_smiles", "mapped", "reason"]], on=key_cols, how="left")
    tot = m[mass_col].sum() * scale
    mapped = m[m.mapped.astype(bool)]
    mass_by = mapped.groupby("rep_smiles")[mass_col].sum() * scale
    unm = m[~m.mapped.astype(bool)].groupby("reason")[mass_col].sum() * scale
    cov = dict(total_mass=tot, mapped_mass=float(mass_by.sum()), unmapped_mass=float(tot - mass_by.sum()),
               mapped_fraction=float(mass_by.sum() / tot), **{f"unmapped[{k}]": float(v) for k, v in unm.items()})
    return list(mass_by.index), (mass_by / mass_by.sum()).to_numpy(), cov


def predict(smi, x, reg):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UnifacFallbackWarning)
        i = flash_point_unifac_detail(smi, x, registry=reg, use_unifac=False)
        u = flash_point_unifac_detail(smi, x, registry=reg, use_unifac=True)
    return i, u


def metrics(err):
    e = np.asarray(err, float)
    e = e[np.isfinite(e)]
    if e.size == 0:
        return dict(N=0)
    return dict(N=int(e.size), mean_bias_C=float(e.mean()), MAE_C=float(np.abs(e).mean()), frac_over_gt5C=float((e > OVER).mean()),
                n_over_gt5C=int((e > OVER).sum()), worst_over_C=float(e.max()), worst_under_C=float(e.min()))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    table = parse_rep_table()
    comp, amap, fl, fmap = build_maps(table)
    # ---- step 2: the maps are on disk before any prediction
    amap.to_csv(OUT / "av3322_representative_map.csv", index=False, float_format="%.6g")
    fmap.to_csv(OUT / "fuellib_representative_map.csv", index=False, float_format="%.6g")

    bulk = pd.read_csv(V5 / "crc_av3322_bulk.csv")
    props = pd.read_csv(V5 / "crc_av3322_properties.csv")
    flp = pd.read_csv(V5 / "fuellib_properties.csv")
    t10_av = props[props.property == "distillation_T10_C"].set_index("fuel_id").value
    t10_fl = flp[flp.property == "distillation_T10_C"].set_index("fuel_id").value
    flash_method = props[props.property == "flash_point_C"].set_index("fuel_id").method

    reg0 = _default_registry()
    best = json.load(open(OUT / "mesitylene_best_input.json")) if (OUT / "mesitylene_best_input.json").exists() else None
    regs = {"registry_as_stored": reg0}
    if best:
        reg1 = reg0.copy()
        can = Chem.MolToSmiles(Chem.MolFromSmiles(MES))
        reg1.loc[can, ["antoine_A", "antoine_B", "antoine_C"]] = best["abc"]
        regs["mesitylene_corrected"] = reg1

    items = []   # (fuel_id, source, df, key_cols, mass_col, scale, mp, measured, method, T10, arom_kinds)
    amp = amap[["group", "kind", "carbon_number", "rep_smiles", "mapped", "reason"]]
    for fid, g in comp.groupby("fuel_id"):
        b = bulk[bulk.fuel_id == fid].iloc[0]
        items.append(dict(fuel_id=fid, source="AV-33-22", df=g[["group", "kind", "carbon_number", "mass_pct"]], keys=["group", "kind", "carbon_number"],
                          mass="mass_pct", scale=0.01, mp=amp, measured=float(b.coa_flash_point_C), method=str(flash_method.get(fid, "")),
                          t10=float(t10_av.get(fid, np.nan)), fuel_type=b.fuel_type))
    fmp = fmap[["gc_bin", "kind", "n_carbon", "rep_smiles", "mapped", "reason"]]
    for s in FUELLIB_FLASH:
        fid = f"FUELLIB_{s}"
        g = fl[fl.fuel_id == fid][["gc_bin", "kind", "n_carbon", "fraction"]]
        meas = flp[(flp.fuel_id == fid) & (flp.property == "flash_point_C")].iloc[0]
        items.append(dict(fuel_id=fid, source="FuelLib", df=g, keys=["gc_bin", "kind", "n_carbon"], mass="fraction", scale=1.0, mp=fmp,
                          measured=float(meas.value), method=str(meas.method), t10=float(t10_fl.get(fid, np.nan)), fuel_type=meas.fuel_name))

    cov_rows, pred_rows = [], []
    for it in items:
        smi, w, cov = fuel_mixture(it["df"], it["keys"], it["mass"], it["scale"], it["mp"])
        arom = float(it["df"][it["df"].kind.isin(AROMATIC_KINDS)][it["mass"]].sum() * it["scale"] / (it["df"][it["mass"]].sum() * it["scale"]))
        included = (1 - cov["mapped_fraction"]) <= MAX_UNMAPPED
        cov_rows.append(dict(fuel_id=it["fuel_id"], source=it["source"], fuel_type=it["fuel_type"], aromatic_mass_fraction=arom,
                             included=included, n_reps=len(smi), **cov))
        if not included:
            continue
        x = mass_to_mole(smi, w)
        assert abs(x.sum() - 1) < 1e-9
        rec = dict(fuel_id=it["fuel_id"], source=it["source"], fuel_type=it["fuel_type"], flash_method=it["method"], measured_C=it["measured"],
                   aromatic_mass_fraction=arom, d86_T10_C=it["t10"], mapped_fraction=cov["mapped_fraction"], n_components=len(smi),
                   mesitylene_mole_fraction=float(sum(xi for s_, xi in zip(smi, x) if s_ == Chem.MolToSmiles(Chem.MolFromSmiles(MES)))))
        for tag, reg in regs.items():
            i, u = predict(smi, x, reg)
            rec[f"{tag}_ideal_C"], rec[f"{tag}_unifac_C"] = i["flash_C"], u["flash_C"]
            rec[f"{tag}_ideal_status"], rec[f"{tag}_unifac_status"] = i["status"], u["status"]
            rec[f"{tag}_ideal_err_C"] = i["flash_C"] - it["measured"]
            rec[f"{tag}_unifac_err_C"] = u["flash_C"] - it["measured"]
        pred_rows.append(rec)
    cov = pd.DataFrame(cov_rows)
    cov.to_csv(OUT / "t2_fuel_composition_coverage.csv", index=False, float_format="%.6g")
    pred = pd.DataFrame(pred_rows)
    pred.to_csv(OUT / "t2_flash_predictions.csv", index=False, float_format="%.6g")

    # ---- R-C1 metrics
    mrows = []
    for tag in regs:
        ok = pred[(pred[f"{tag}_ideal_status"] == "ok") & (pred[f"{tag}_unifac_status"] == "ok")]
        groups = {"all_eligible": ok, "FuelLib-4": ok[ok.source == "FuelLib"], "AV-33-22": ok[ok.source == "AV-33-22"]}
        for meth, d in ok[ok.source == "AV-33-22"].groupby("flash_method"):
            groups[f"AV-33-22 CoA method={meth}"] = d
        for gname, d in groups.items():
            for model in ("ideal", "unifac"):
                mrows.append(dict(registry=tag, group=gname, model=model, **metrics(d[f"{tag}_{model}_err_C"])))
    met = pd.DataFrame(mrows)
    met.to_csv(OUT / "t2_flash_metrics.csv", index=False, float_format="%.6g")

    # ---- readings (primary = registry_as_stored)
    P = "registry_as_stored"
    ok = pred[(pred[f"{P}_ideal_status"] == "ok") & (pred[f"{P}_unifac_status"] == "ok")].copy()
    av = ok[ok.source == "AV-33-22"]
    rd = []
    m_av_i = metrics(av[f"{P}_ideal_err_C"])
    rd.append(dict(reading="R-C2", item="AV-33-22 ideal mean bias (C)", value=m_av_i["mean_bias_C"], threshold="<= +2.0", met=bool(m_av_i["mean_bias_C"] <= 2.0)))
    rd.append(dict(reading="R-C2", item="AV-33-22 ideal fraction over-predicted >5 C", value=m_av_i["frac_over_gt5C"], threshold="<= 0.20", met=bool(m_av_i["frac_over_gt5C"] <= 0.20)))
    rd.append(dict(reading="R-C2", item="overall (both conditions)", value=np.nan, threshold="both", met=bool(rd[-1]["met"] and rd[-2]["met"])))
    rd.append(dict(reading="R-C2 ref", item="r3 dev reference bias / frac>5C", value="+0.44 / 0.13", threshold="", met=""))
    rd.append(dict(reading="R-C2 ref", item="r2 reference bias / frac>5C", value="+4.22 / 0.53", threshold="", met=""))
    # pairing with CRC 647 (v4_overlap)
    fl_comp = pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv")
    fl_ov = fl_comp[fl_comp.fuel_id.isin([f"FUELLIB_{s}" for s in FUELLIB_FLASH])].groupby("fuel_id").v4_overlap.agg(lambda s: ";".join(sorted(set(s.dropna()))))
    pair = pd.DataFrame({"fuel_id": fl_ov.index, "v4_overlap_value": fl_ov.values})
    pair["overlaps_CRC647"] = pair.v4_overlap_value.str.startswith("CRC647")
    pair = pd.concat([pair, pd.DataFrame([dict(fuel_id="CRC_AV3322_01..38", v4_overlap_value="(no v4_overlap field in AV-33-22 files)", overlaps_CRC647=False)])])
    pair.to_csv(OUT / "t2_pairing.csv", index=False)
    n_pair = int(pair.overlaps_CRC647.sum())
    rd.append(dict(reading="R-C2 pairing", item="T2 fuels also in CRC 647", value=n_pair, threshold=">= 3 to pair", met=bool(n_pair >= 3)))
    v4p = pd.read_csv(ROOT / "results/v4/R3/F/flash_real_fuels.csv")
    v4p = v4p[v4p.fuel_id == "POSF-10325"].iloc[0]
    p10325 = pred[pred.fuel_id == "FUELLIB_posf10325"].iloc[0]
    rd.append(dict(reading="info", item="POSF-10325: v4 r3 class x cn (D2887-derived) ideal err / measured used (C)", value=f"{v4p.ideal_err_C:.2f} / {v4p.measured_C}", threshold="", met=""))
    rd.append(dict(reading="info", item="POSF-10325: FuelLib GCxGC ideal err / measured used (C)", value=f"{p10325[f'{P}_ideal_err_C']:.2f} / {p10325.measured_C}", threshold="", met=""))

    # R-C3
    med = ok.aromatic_mass_fraction.median()
    up = ok[ok.aromatic_mass_fraction > med]
    ei, eu = up[f"{P}_ideal_err_C"].mean(), up[f"{P}_unifac_err_C"].mean()
    rd.append(dict(reading="R-C3", item="median aromatic mass fraction (eligible fuels)", value=med, threshold="", met=""))
    rd.append(dict(reading="R-C3", item=f"upper half n={len(up)}: mean err ideal (C)", value=ei, threshold="", met=""))
    rd.append(dict(reading="R-C3", item=f"upper half n={len(up)}: mean err UNIFAC (C)", value=eu, threshold="opposite signs", met=bool(np.sign(ei) != np.sign(eu))))
    lo = ok[ok.aromatic_mass_fraction <= med]
    rd.append(dict(reading="R-C3 info", item=f"lower half n={len(lo)}: mean err ideal / UNIFAC (C)", value=f"{lo[f'{P}_ideal_err_C'].mean():.3f} / {lo[f'{P}_unifac_err_C'].mean():.3f}", threshold="", met=""))
    # R-C4
    def pear(d, col):
        d = d.dropna(subset=["d86_T10_C", col])
        if len(d) < 3:
            return len(d), np.nan, np.nan
        r, p = stats.pearsonr(d.d86_T10_C, d[col])
        return len(d), float(r), float(p)
    for src, d in [("FuelLib-4 (declared)", ok[ok.source == "FuelLib"]), ("AV-33-22 (not declared; file has D86 T10)", av), ("all eligible (info)", ok)]:
        for model in ("ideal", "unifac"):
            n, r, p = pear(d, f"{P}_{model}_err_C")
            rd.append(dict(reading="R-C4", item=f"{src} {model}: Pearson r(err, D86 T10)", value=r, threshold=f"n={n}; p={p:.3g}" if np.isfinite(p) else f"n={n}", met=""))
    rdf = pd.DataFrame(rd)
    rdf.to_csv(OUT / "t2_readings.csv", index=False, float_format="%.6g")

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    pd.set_option("display.max_colwidth", 80)
    print(amap[~amap.mapped].to_string())
    print(fmap[~fmap.mapped].to_string())
    print(cov[["fuel_id", "source", "mapped_fraction", "included", "aromatic_mass_fraction", "n_reps"]].to_string())
    print(met.to_string())
    print(rdf.to_string())
    print(pred[["fuel_id", "flash_method", "measured_C", f"{P}_ideal_C", f"{P}_unifac_C", "mesitylene_corrected_ideal_C" if best else P + "_ideal_C", "mesitylene_mole_fraction"]].to_string())


if __name__ == "__main__":
    main()
