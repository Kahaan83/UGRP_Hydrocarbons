"""v4 Task 7: group splits over the v4-clean real-fuel corpora.

Output: splits/v4/real_fuels_v4_splits.json, results/v4/P2/coverage_matrix.csv.

Universe (ruling R1): fuel_ids of inventory_v4.csv  UNION  included jet-GC fuels (``included_jet_fuel_ids``, which
drops the partial JP5/JP8 Potter-Simmons rows)  UNION  gasoline blend_ids.  Disjointness of the inventory and
jet-GC ids is asserted, never assumed.

GROUPS: see GROUP_DEFINITION.

LOCKED SET
  20 % of GROUPS (seed 20261030), chosen FIRST, stratified by source family (largest remainder; families with < 5
  groups contribute none), then 5-fold GroupKFold (sample = fuel_id, shuffled with the same seed) over the
  remaining groups.  Nothing tuned may touch the locked set.
  The headline set (headline_ok) is n=1: ``headline_statistics_allowed`` is False - never report statistics over it.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.data.loaders_jet import included_jet_fuel_ids  # noqa: E402

V4 = ROOT / "data" / "real_fuels" / "v4"
SEED = 20261030
N_FOLDS = 5
LOCKED_FRAC = 0.20
OUT_SPLITS = ROOT / "splits" / "v4" / "real_fuels_v4_splits.json"
OUT_COV = ROOT / "results" / "v4" / "P2" / "coverage_matrix.csv"
GROUP_DEFINITION = (
    "The leakage unit; a group is never split across locked/folds. CRC647_<n> and jet-GC fuels: one group per "
    "fuel_id. NIST_RP1 + its 3 RP-1 surrogates (Binary_2comp, Edwards_3comp, NIST_14comp): ONE group 'RP1' "
    "(surrogates model the same physical fuel). FACE gasolines: group = base_fuel (ethanol blends stay with their "
    "base). Co-Optima bRON blends (blend_id BRON:<CAS>:<frac>): group = 'BRON:<CAS>' of the blendstock, across both "
    "base-gasoline lumps (R19): all fractions of one blendstock stay together, so locked/fold bRON blendstocks are "
    "unseen at train time. The base lump's RON/MON stays a model INPUT (R17); this is a stricter unseen-blendstock "
    "test, and a blendstock tested on both bases is held out for both."
)
RP1_GROUP = ["NIST_RP1", "RP1_Binary_2comp", "RP1_Edwards_3comp", "RP1_NIST_14comp"]
DIESEL_RE = re.compile(r"diesel|lfp|ulsd|en590", re.I)
DIESEL_RAW_DIRS = [ROOT / "data" / "raw" / "pour_point", ROOT / "data" / "raw" / "gc_composition",
                   ROOT / "extras" / "fuelblend_data" / "data" / "raw" / "pour_point",
                   ROOT / "extras" / "fuelblend_data" / "data" / "raw" / "wax"]


DIESEL_SKIPPED_FILES: list[str] = []


def diesel_qualifying_fuels() -> list[str]:
    """Diesel fuels with a raw-file CLASS/COMPONENT composition (``class_or_component`` + ``pct`` columns)."""
    found: set[str] = set()
    DIESEL_SKIPPED_FILES.clear()
    for d in DIESEL_RAW_DIRS:
        for f in sorted(d.glob("*.csv")) if d.exists() else []:
            try:
                df = pd.read_csv(f, dtype=str)
            except (pd.errors.ParserError, UnicodeDecodeError, OSError) as exc:
                DIESEL_SKIPPED_FILES.append(f"{f.name}: {type(exc).__name__}")
                print(f"[diesel scan] skipped unreadable {f}: {exc}", file=sys.stderr)
                continue
            if not {"fuel_code", "class_or_component", "pct"} <= set(df.columns):
                continue
            names = df["fuel_code"].fillna("")
            if "fuel_name" in df.columns:
                names = names + " " + df["fuel_name"].fillna("")
            found |= set(df.loc[names.str.contains(DIESEL_RE), "fuel_code"])
    return sorted(found)


def _family(src: str) -> str:
    if src.startswith("cooptima_bron"):
        return "cooptima_bron"
    if src.startswith("face:face_gasoline_set_II"):
        return "face_set_II"
    if src.startswith("face:"):
        return "face_set_I"
    return src


def _load():
    inv = pd.read_csv(V4 / "inventory_v4.csv")
    jet = pd.read_csv(V4 / "jet_class_compositions.csv")
    gas = pd.read_csv(V4 / "gasoline_octane_blends.csv")
    included = included_jet_fuel_ids(jet)
    assert not set(inv.fuel_id) & set(included), "inventory and jet-GC ids overlap (double count)"
    assert not set(gas.blend_id) & (set(inv.fuel_id) | set(included)), "gasoline ids overlap other corpora"
    return inv, jet, gas, included


def _units(inv, gas, included) -> pd.DataFrame:
    rows = []
    for fid in inv.fuel_id:
        in_rp1 = fid in RP1_GROUP
        rows.append((fid, "RP1" if in_rp1 else fid, "nist_rp1" if in_rp1 else "crc647", "inventory_v4"))
    for fid in included:
        rows.append((fid, fid, "jet_gc", "jet_class_compositions"))
    for bid, base, src in zip(gas.blend_id, gas.base_fuel, gas.source):
        grp = ":".join(bid.split(":")[:2]) if bid.startswith("BRON:") else base
        rows.append((bid, grp, _family(src), "gasoline_octane_blends"))
    u = pd.DataFrame(rows, columns=["fuel_id", "group", "family", "corpus"])
    assert u.fuel_id.is_unique, "duplicate fuel_id in the universe"
    norm = u.fuel_id.str.replace(r"[-_\s]", "", regex=True).str.upper()
    assert norm.is_unique, f"same physical fuel under two spellings: {u.fuel_id[norm.duplicated(keep=False)].tolist()}"
    assert not u.fuel_id.str.startswith("POSF_").any(), "v3-style POSF_ id present"
    gf = u.groupby("group").family.nunique()
    assert (gf == 1).all(), f"group spans families: {gf[gf > 1].index.tolist()}"
    return u


def _choose_locked(groups: pd.DataFrame, seed: int = SEED) -> list[str]:
    n_lock = int(round(LOCKED_FRAC * len(groups)))
    fam_n = groups.groupby("family").size()
    elig = fam_n[fam_n >= 5]
    quota = elig / elig.sum() * n_lock
    alloc = quota.astype(int)
    rem = (quota - alloc).sort_values(ascending=False, kind="stable")
    for fam in rem.index[: n_lock - int(alloc.sum())]:
        alloc[fam] += 1
    rng = np.random.default_rng(seed)
    out: list[str] = []
    for fam in sorted(alloc.index):
        gl = sorted(groups.loc[groups.family == fam, "group"])
        pick = rng.choice(len(gl), size=int(alloc[fam]), replace=False)
        out += [gl[i] for i in sorted(pick)]
    return sorted(out)


def _coverage(u, jet, gas) -> pd.DataFrame:
    jp = pd.read_csv(V4 / "jet_properties_v4.csv")
    jflags = jet.set_index("fuel_id")
    meta = u.set_index("fuel_id")
    gas_i = gas.set_index("blend_id")
    rows = []
    for fid in u.fuel_id:
        if fid in jflags.index:
            r = jflags.loc[fid]
            hl, ci, cc = bool(r.headline_ok), bool(r.cn_split_class_independent), bool(r.composition_complete)
        else:
            hl, ci, cc = False, False, None
        base = {"fuel_id": fid, "group": meta.at[fid, "group"], "source_family": meta.at[fid, "family"],
                "headline_ok": hl, "cn_split_class_independent": ci, "composition_complete": cc}
        rows_f = []
        for (prop, role), g in jp[jp.fuel_id == fid].groupby(["property", "role"]):
            rows_f.append({**base, "property": prop, "present": True, "role": role, "n_values": len(g),
                           "source": "+".join(sorted(g.source_file.unique())),
                           "provenance": "+".join(sorted(g.provenance.unique()))})
        if fid in gas_i.index:
            r = gas_i.loc[fid]
            for prop, ok in (("ron", r.target_ok_ron), ("mon", r.target_ok_mon)):
                if pd.notna(r[prop]):
                    rows_f.append({**base, "property": prop, "present": True,
                                   "role": "target" if ok else "descriptor", "n_values": 1,
                                   "source": r.source, "provenance": r.provenance})
        if not rows_f:
            rows_f.append({**base, "property": "<none>", "present": False, "role": "", "n_values": 0,
                           "source": "", "provenance": ""})
        rows += rows_f
    return pd.DataFrame(rows).sort_values(["fuel_id", "property"], kind="stable").reset_index(drop=True)


def build(seed: int = SEED) -> dict:
    inv, jet, gas, included = _load()
    u = _units(inv, gas, included)
    groups = u.drop_duplicates("group")[["group", "family"]].reset_index(drop=True)
    locked_groups = _choose_locked(groups, seed)
    dev = u[~u.group.isin(locked_groups)].sort_values("fuel_id").reset_index(drop=True)
    role = {f: "locked" for f in u.fuel_id[u.group.isin(locked_groups)]}
    gkf = GroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
    for k, (_, te) in enumerate(gkf.split(dev, groups=dev.group)):
        for f in dev.fuel_id.iloc[te]:
            role[f] = f"fold{k}"
    folds = {f"fold{k}": sorted(f for f, r in role.items() if r == f"fold{k}") for k in range(N_FOLDS)}
    members = {g: sorted(s.fuel_id) for g, s in u.groupby("group")}
    hl = sorted(jet.loc[jet.headline_ok & jet.fuel_id.isin(included), "fuel_id"])
    meta = u.set_index("fuel_id")
    cov = _coverage(u, jet, gas)
    cov["split_role"] = cov.fuel_id.map(role)
    tgt = cov[(cov.role == "target") & cov.present]
    tcount = {sp: {k: int(v) for k, v in tgt[tgt.split_role == sp].groupby("property").fuel_id.nunique().items()}
              for sp in ["locked"] + [f"fold{k}" for k in range(N_FOLDS)]}
    fsz = sorted(len(v) for v in folds.values())
    med = float(np.median(fsz))
    locked_set = [g for g in hl if role[g] == "locked"]
    gc_locked = sorted(f for f in included if role[f] == "locked")
    splits = {
        "seed": seed, "n_folds": N_FOLDS, "locked_fraction_of_groups": LOCKED_FRAC,
        "group_definition": GROUP_DEFINITION,
        "procedure": "locked groups chosen first (family-stratified, largest remainder, families < 5 groups "
                     "contribute none); GroupKFold(shuffle=True, random_state=seed) over the remaining fuel_ids "
                     "grouped by `group`",
        "universe": {"inventory_v4": int((u.corpus == "inventory_v4").sum()),
                     "jet_class_compositions_included": int((u.corpus == "jet_class_compositions").sum()),
                     "jet_class_compositions_excluded": sorted(set(jet.fuel_id) - set(included)),
                     "gasoline_octane_blends": int((u.corpus == "gasoline_octane_blends").sum())},
        "counts": {"fuels": len(u), "groups": len(members),
                   "locked_fuels": sum(r == "locked" for r in role.values()), "locked_groups": len(locked_groups),
                   "locked_fraction_of_fuels": round(sum(r == "locked" for r in role.values()) / len(u), 4),
                   "locked_fraction_of_groups_actual": round(len(locked_groups) / len(members), 4),
                   "fold_size_max_over_median": round(max(fsz) / med, 3),
                   "target_fuels_by_split_and_property": tcount,
                   "diesel_scan_skipped_files": list(DIESEL_SKIPPED_FILES),
                   "fold_fuels": {k: len(v) for k, v in folds.items()},
                   "groups_by_family": groups.family.value_counts().sort_index().to_dict(),
                   "fuels_by_family": u.family.value_counts().sort_index().to_dict()},
        "locked": sorted(f for f, r in role.items() if r == "locked"),
        "locked_groups": locked_groups,
        "folds": folds,
        "fuel_role": dict(sorted(role.items())),
        "fuel_group": dict(sorted(meta["group"].items())),
        "fuel_source_family": dict(sorted(meta["family"].items())),
        "groups": members,
        "headline_set": hl,
        "headline_statistics_allowed": False,
        "headline_note": f"headline_ok has n={len(hl)} (cn_split_class_independent=True); report no statistics over it. "
                         f"Headline fuel(s) in locked set: {locked_set or 'none'} (headline fuel POSF-10325 is in dev "
                         f"fold {role.get('POSF-10325')}, NOT held out); jet-GC fuels locked: {gc_locked}.",
        "diesel": {"qualifying_fuels": diesel_qualifying_fuels(), "decision": "SD-02 (docs/scope_decisions.md)"},
    }
    return {"splits": splits, "coverage": cov.drop(columns="split_role")}


def main() -> None:
    r = build()
    OUT_SPLITS.parent.mkdir(parents=True, exist_ok=True)
    OUT_COV.parent.mkdir(parents=True, exist_ok=True)
    OUT_SPLITS.write_text(json.dumps(r["splits"], indent=1, sort_keys=True) + "\n", encoding="utf-8")
    r["coverage"].to_csv(OUT_COV, index=False)
    c = r["splits"]["counts"]
    md = ["# Coverage matrix summary (v4 splits)", "", "Locked-set / per-fold TARGET fuels by property "
          "(fuels with role=target in coverage_matrix.csv):", ""]
    for sp, d in c["target_fuels_by_split_and_property"].items():
        md.append(f"- **{sp}**: " + (", ".join(f"{k}={v}" for k, v in sorted(d.items())) or "none"))
    md += ["", f"Fold sizes (fuels): {c['fold_fuels']}; max/median = {c['fold_size_max_over_median']}.",
           f"Locked: {c['locked_groups']} groups, {c['locked_fuels']} fuels "
           f"({c['locked_fraction_of_fuels']:.1%} of fuels, {c['locked_fraction_of_groups_actual']:.1%} of groups).", "",
           r["splits"]["headline_note"], ""]
    (OUT_COV.parent / "coverage_matrix_summary.md").write_text(chr(10).join(md), encoding="utf-8")
    print(json.dumps(c, indent=1))


if __name__ == "__main__":
    main()
