"""v6 flash arm step 1: stage the Leenhouts flash mixtures (all chemistries), attach measured pure flash points, drop rows with
pure coverage < 0.9 (logged, never filled), assign the outcome-blind dev/locked split per chemistry class.
Writes results/v6/flash_oxy/rows.csv, lock.json, build_log.json and the staged canonical table. No model is fitted here."""
import sys, io, json, hashlib, zipfile
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
R = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(R / "src"))
from fuelblend.data import loaders_leenhouts_flash as L   # noqa: E402

ZIP = R / "extras/fuelblend_data/data/raw/zenodo/leenhouts_benchmark/datasets.zip"
PURE_MEMBER = "datasets/flashpoint/processed_data/flashpoint_pure_exp.csv"
OUT = R / "results/v6/flash_oxy"
STG = R / "data/interim/external_staging"
LOCK_ROW_MAX, COV_MIN = 0.25, 0.9
log = {}


def pure_flash():
    with zipfile.ZipFile(ZIP) as z:
        p = pd.read_csv(io.BytesIO(z.read(PURE_MEMBER)))
    assert set(p["unit"].unique()) <= {"Celcius"}
    p["ik"] = [Chem.MolToInchiKey(m) if (m := Chem.MolFromInchi(i)) is not None else None for i in p.component_inchi_0]
    p = p.dropna(subset=["ik", "value"])
    log["pure"] = dict(rows=int(len(p)), compounds=int(p.ik.nunique()))
    return p.groupby("ik").value.median(), {ik: Chem.MolToSmiles(Chem.MolFromInchi(i)) for ik, i in zip(p.ik, p.component_inchi_0)}


def smiles_for(d, pure_smi):
    """InChIKey -> SMILES for every mixture component, read from the benchmark file itself."""
    raw = L._read_raw(ZIP)
    m = {}
    for k in range(L.MAX_COMP):
        for inchi in raw[f"component_inchi_{k}"].dropna().unique():
            mol = Chem.MolFromInchi(inchi.strip())
            if mol is not None:
                m[Chem.MolToInchiKey(mol)] = Chem.MolToSmiles(mol)
    return m


def assign_roles(df):
    lock = {}
    for cls, g in df.groupby("chem_class"):
        rows_per = g.groupby("system_id").size()
        groups = sorted(rows_per.index, key=lambda s: hashlib.md5(f"v6f:{s}".encode()).hexdigest())
        locked, n = [], 0
        for gp in groups:
            if n + rows_per[gp] <= LOCK_ROW_MAX * len(g):
                locked.append(gp)
                n += rows_per[gp]
        lock[cls] = dict(n_groups=len(groups), locked=locked)
    df["role"] = [("locked" if s in lock[c]["locked"] else "dev") for c, s in zip(df.chem_class, df.system_id)]
    return lock


if __name__ == "__main__":
    mix = L.build(ZIP, log.setdefault("loader", {}))
    mix.to_csv(STG / "canonical_leenhouts_flash_mix_measurements.csv", index=False)
    pf, pure_smi = pure_flash()
    smi = smiles_for(mix, pure_smi)
    out, low = [], []
    for r in mix.itertuples():
        ks, x = json.loads(r.comp_inchikeys), np.array(json.loads(r.x_mole), float)
        have = np.array([k in pf.index for k in ks])
        cov = float(x[have].sum())
        if cov < COV_MIN:
            low.append(dict(fuel_id=r.fuel_id, system_id=r.system_id, chem_class=r.chem_class, cov=round(cov, 3)))
            continue
        base = float((x[have] * np.array([pf[k] for k, h in zip(ks, have) if h])).sum() / cov)
        out.append(dict(fuel_id=r.fuel_id, chem_class=r.chem_class, system_id=r.system_id, n_comp=r.n_comp, y=r.value, base_linear=base,
                        cov=cov, comps=r.comp_inchikeys, smiles=json.dumps([smi[k] for k in ks]), x=r.x_mole, w=r.w_mass,
                        pure_flash=json.dumps([float(pf[k]) if k in pf.index else None for k in ks])))
    df = pd.DataFrame(out)
    log["coverage"] = dict(kept=int(len(df)), dropped_low_coverage=len(low), cov_min=COV_MIN, dropped=low,
                           pure_missing_components=sorted({k for r in mix.itertuples() for k in json.loads(r.comp_inchikeys) if k not in pf.index}))
    lock = assign_roles(df)
    df = df.sort_values(["chem_class", "system_id", "fuel_id"], kind="stable").reset_index(drop=True)
    df.to_csv(OUT / "rows.csv", index=False)
    summ = df.groupby(["chem_class", "role"]).agg(rows=("y", "size"), systems=("system_id", "nunique")).reset_index()
    strat = df[df.n_comp >= 3].groupby(["chem_class", "role"]).system_id.nunique().reset_index(name="systems_n3")
    log["counts"] = summ.to_dict("records"); log["counts_n3_systems"] = strat.to_dict("records")
    json.dump(dict(row_max=LOCK_ROW_MAX, rule="md5('v6f:'+system_id) ascending; per chem_class; lock a system while locked rows <= row_max of the class's rows",
                   classes=lock, rows_csv_sha256=hashlib.sha256((OUT / "rows.csv").read_bytes()).hexdigest()), open(OUT / "lock.json", "w"), indent=1)
    json.dump(log, open(OUT / "build_log.json", "w"), indent=1, default=str)
    print(summ.to_string()); print(strat.to_string())
    print({k: v for k, v in log["coverage"].items() if k not in ("dropped",)})
