"""Track C step 1: build the deduplicated mixture rows (KAUST Table 4, Leenhouts/Kuzhagaliyeva MON mixtures, Boehm 2022 cold
viscosity/density), attach the physics baselines, and assign the dev/locked group split. Writes results/v6/mixture/rows.csv,
lock.json and build_log.json. Does not touch data/unified/, splits/ or v5 BASE.

Octane heads deliberately ignore the is_hc filter (oxygenates, e.g. ethanol, are in scope). Staged files are raw rows, so rows are
deduplicated here on (head, prop, group, composition rounded to 3 decimals, T) with the median value (AGENTS.md dedup rule)."""
import sys, json, math, hashlib
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
R = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(R / "scripts/external_ingest")); sys.path.insert(0, str(R / "src"))
from _common import mol_info                                                # noqa: E402
from kaust_registry import COMPONENTS, densities, load_tables, pure_props   # noqa: E402
from fuelblend.data import loaders_boehm2022 as B                           # noqa: E402

OUT = R / "results/v6/mixture"
STG = R / "data/interim/external_staging"
LOCK_ROW_MAX = 0.25
COV_MIN = 0.9
log = {}


def canon_ik(smi):
    return Chem.MolToInchiKey(Chem.MolFromSmiles(smi))


def kaust_rows():
    t3, t4, _ = load_tables()
    pp, dens = pure_props(t3, t4), densities()
    mw = {k: mol_info(v["smiles"])[1] for k, v in COMPONENTS.items()}
    ik = {k: canon_ik(v["smiles"]) for k, v in COMPONENTS.items()}
    out = []
    for r in t4.itertuples():
        comp = json.loads(r.composition_vol_pct)
        if len(comp) < 2:
            continue  # single-component anchor rows define the pure values
        tot = sum(comp.values())
        v = {k: x / tot for k, x in comp.items()}
        n = {k: v[k] * dens[k][0] / mw[k] for k in v}
        nt = sum(n.values())
        x = {k: n[k] / nt for k in n}
        grp = "+".join(sorted(c for c in comp if c != "ethanol"))
        names = sorted(x)
        for prop in ("ron", "mon"):
            out.append(dict(head="kaust", prop=prop, group=grp, system_id=f"kaust:{grp}", n_comp=len(comp), T_K=np.nan,
                            y=float(getattr(r, prop)), base=sum(x[k] * pp[k][prop] for k in x), cov=1.0,
                            comps=json.dumps([ik[k] for k in names]), smiles=json.dumps([COMPONENTS[k]["smiles"] for k in names]),
                            x=json.dumps([x[k] for k in names]), w="[]"))
    return pd.DataFrame(out)


def leenhouts_rows(pool_mon):
    d = pd.read_csv(STG / "canonical_leenhouts_mon_mix_measurements.csv")
    cmp_ = pd.read_csv(STG / "canonical_leenhouts_mon_compounds.csv")
    smi = dict(zip(cmp_.inchikey, cmp_.smiles))
    recs = []
    for r in d.itertuples():
        ks = json.loads(r.comp_inchikeys)
        x = np.array(json.loads(r.x_mole), float)
        x = x / x.sum()
        recs.append((r.system_id, tuple(ks), tuple(np.round(x, 3)), float(r.value)))
    df = pd.DataFrame(recs, columns=["system_id", "ks", "xr", "y"])
    n_raw = len(df)
    df = df.groupby(["system_id", "ks", "xr"], as_index=False).y.median()
    log["leenhouts"] = dict(raw_rows=n_raw, dedup_rows=len(df))
    out, low = [], 0
    for r in df.itertuples():
        ks, x = list(r.ks), np.array(r.xr, float)
        x = x / x.sum()
        have = np.array([k in pool_mon.index for k in ks])
        cov = float(x[have].sum())
        if cov < COV_MIN:
            low += 1
            continue
        base = float((x[have] * np.array([pool_mon[k] for k, h in zip(ks, have) if h])).sum() / cov)
        out.append(dict(head="leenhouts", prop="mon", group=r.system_id, system_id=f"leenhouts:{r.system_id}", n_comp=len(ks), T_K=np.nan,
                        y=r.y, base=base, cov=cov, comps=json.dumps(ks), smiles=json.dumps([smi[k] for k in ks]),
                        x=json.dumps(list(x)), w="[]"))
    log["leenhouts"].update(dropped_cov_below_min=low, kept=len(out), cov_min=COV_MIN)
    return pd.DataFrame(out)


def boehm_rows():
    b = pd.read_csv(R / "data/real_fuels/v5/boehm2022_viscosity.csv")
    ik2smi = {canon_ik(s): s for s in B.SMILES.values()}
    b = b[b.property.isin(["viscosity_kin_cSt", "density_kg_m3"])]
    pure = b[b.n_comp == 1].copy()
    pure["ik"] = pure.comp_inchikeys.map(lambda s: json.loads(s)[0])
    pv = {(p, ik, round(T, 2)): v for p, ik, T, v in zip(pure.property, pure.ik, pure.T_K, pure.value)}
    out, skipped = [], 0
    for r in b[b.n_comp >= 2].itertuples():
        ks = json.loads(r.comp_inchikeys)
        w = np.array(json.loads(r.w_mass), float)
        x = np.array(json.loads(r.x_mole), float)
        if any(k not in ik2smi for k in ks):
            skipped += 1
            continue
        pu = [pv.get((r.property, k, round(r.T_K, 2))) for k in ks]
        if any(p is None for p in pu):
            skipped += 1
            continue
        pu = np.array(pu, float)
        if r.property == "viscosity_kin_cSt":
            vbn = 14.534 * np.log(np.log(pu + 0.8)) + 10.975
            vb = float((w * vbn).sum())
            nu = math.exp(math.exp((vb - 10.975) / 14.534)) - 0.8
            y, base, prop = math.log(r.value), math.log(nu), "ln_visc"
        else:
            y, base, prop = float(r.value), float(1.0 / (w / pu).sum()), "density"
        out.append(dict(head="boehm", prop=prop, group=r.system_id, system_id=f"boehm:{r.system_id}", n_comp=len(ks), T_K=round(r.T_K, 2),
                        y=y, base=base, cov=1.0, comps=json.dumps(ks), smiles=json.dumps([ik2smi[k] for k in ks]),
                        x=json.dumps(list(x)), w=json.dumps(list(w))))
    log["boehm"] = dict(mixture_rows_in=int((b.n_comp >= 2).sum()), kept=len(out), skipped_no_pure_or_smiles=skipped)
    return pd.DataFrame(out)


def assign_roles(df):
    """Outcome-blind lock: groups in md5('v6c:'+group) order are locked while the locked share of rows stays <= LOCK_ROW_MAX
    (a group that would exceed it is skipped), with at least one group locked."""
    lock = {}
    for head, g in df.groupby("head"):
        rows_per = g.groupby("group").size()
        groups = sorted(rows_per.index, key=lambda s: hashlib.md5(f"v6c:{s}".encode()).hexdigest())
        locked, n = [], 0
        for gp in groups:
            if n + rows_per[gp] <= LOCK_ROW_MAX * len(g):
                locked.append(gp)
                n += rows_per[gp]
        lock[head] = dict(n_groups=len(groups), locked=locked)
    df["role"] = [("locked" if gp in lock[h]["locked"] else "dev") for h, gp in zip(df["head"], df["group"])]
    return lock


if __name__ == "__main__":
    pool = pd.read_csv(R / "results/v6/gnn/pool.csv")
    pool_mon = pool[pool.target == "mon"].set_index("inchikey").value
    df = pd.concat([kaust_rows(), leenhouts_rows(pool_mon), boehm_rows()], ignore_index=True)
    df["xr"] = df.x.map(lambda s: str(tuple(np.round(json.loads(s), 3))))
    n0 = len(df)
    key = ["head", "prop", "group", "comps", "xr", "T_K"]
    agg = {c: "first" for c in df.columns if c not in key + ["y"]}
    agg["y"] = "median"
    df = df.groupby(key, dropna=False, as_index=False).agg(agg).drop(columns=["xr"])
    log["dedup_total"] = dict(before=n0, after=len(df))
    lock = assign_roles(df)
    df = df.sort_values(["head", "prop", "group", "T_K"], kind="stable").reset_index(drop=True)
    df.to_csv(OUT / "rows.csv", index=False)
    summ = df.groupby(["head", "prop", "role"]).agg(n=("y", "size"), groups=("group", "nunique")).reset_index()
    log["counts"] = summ.to_dict("records")
    json.dump(dict(row_max=LOCK_ROW_MAX, rule="md5('v6c:'+group) ascending; lock a group while locked rows <= row_max of the head's rows; per head", heads=lock,
                   rows_csv_sha256=hashlib.sha256((OUT / "rows.csv").read_bytes()).hexdigest()), open(OUT / "lock.json", "w"), indent=1)
    json.dump(log, open(OUT / "build_log.json", "w"), indent=1, default=str)
    print(summ.to_string())
    print(json.dumps({k_: v for k_, v in log.items() if k_ != "counts"}, default=str))
