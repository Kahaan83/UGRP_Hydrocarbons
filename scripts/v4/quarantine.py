#!/usr/bin/env python
"""v4 Task 3: quarantine non-raw rows of the v2/v3 real-fuel tables and build the clean v4 tables.

Inputs: results/v4/P1/provenance_ledger.csv (Task 1), results/v4/P1/identity_mismatches.csv (Task 2).

Keep rule (per ledger row):
  raw_file                      -> provenance "raw_file",      verified=True
  cited_typed                   -> provenance "cited_typed",   verified=False
  formula_derived whose evidence starts with "reconstructed:" (CRC 647 reconstructed
  compositions, ruling R10)     -> provenance "reconstructed", verified=False
Everything else (formula_derived, target_derived_composition, placeholder_constant, unknown) is
quarantined. A composition row whose InChIKey does not belong to its SMILES (identity_mismatches.csv,
joined on (table, row_id)) is quarantined too, whatever its provenance class. Fuels left with no
clean non-lump composition row are dropped from v4 (all their rows are quarantined with reason
``fuel_dropped:no_clean_composition``) and listed in data/MISSING.md and
data/quarantine/dropped_fuels.csv.

Outputs
  data/quarantine/inventory_quarantined.csv     long: ledger row + reason, every audited inventory table
  data/quarantine/compositions_quarantined.csv  full composition row + table,row_id,provenance_class,reason
  data/quarantine/dropped_fuels.csv
  data/real_fuels/v4/inventory_v4.csv           wide, one row per kept fuel; quarantined cells are NaN;
                                                columns provenance (JSON field->label), unverified_fields,
                                                verified, n_composition_rows_quarantined
  data/real_fuels/v4/compositions_v4.csv        kept composition rows + provenance, verified
v4 is built from the v3 pair (inventory_v3.csv / compositions_v3.csv), a superset of the v2 tables.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

V3_INV = "data/real_fuels/v3/inventory_v3.csv"
V3_COMP = "data/real_fuels/v3/compositions_v3.csv"
PAIRS = [
    ("data/real_fuels/inventory.csv", "data/real_fuels/compositions.csv"),
    ("data/real_fuels/inventory_measured.csv", "data/real_fuels/compositions_measured.csv"),
    ("data/real_fuels/inventory_reconstructed.csv", "data/real_fuels/compositions_reconstructed.csv"),
    (V3_INV, V3_COMP),
]
LEDGER = "results/v4/P1/provenance_ledger.csv"
MISMATCH = "results/v4/P1/identity_mismatches.csv"
QDIR = "data/quarantine"
VDIR = "data/real_fuels/v4"
MARK_BEGIN = "<!-- v4-task3-dropped:begin -->"
MARK_END = "<!-- v4-task3-dropped:end -->"
NO_CLEAN = "fuel_dropped:no_clean_composition"


def keep_label(cls: str, evidence: str) -> str | None:
    """Provenance label of a kept ledger row, or None when the row is quarantined."""
    if cls == "raw_file":
        return "raw_file"
    if cls == "cited_typed":
        return "cited_typed"
    if cls == "formula_derived" and str(evidence).startswith("reconstructed:"):
        return "reconstructed"
    return None


def comp_ledger_keys(comp: pd.DataFrame) -> pd.Series:
    """(fuel_id, field) key of each composition row, as the Task 1 ledger names it."""
    ordn = comp.groupby(["system_id", "inchikey"]).cumcount()
    fid = comp["system_id"].astype(str).str.replace(r"^fuel:", "", regex=True)
    return fid + "|fraction:" + comp["inchikey"].astype(str) + "#" + ordn.astype(str)


def split_composition(comp: pd.DataFrame, comp_ledger: pd.DataFrame, mismatches: pd.DataFrame,
                      table: str) -> pd.DataFrame:
    """Annotate every composition row with provenance_class, evidence, label, reason (``""`` = clean)."""
    led = comp_ledger.assign(_k=comp_ledger["fuel_id"].astype(str) + "|" + comp_ledger["field"])
    assert not led["_k"].duplicated().any()
    idx = led.set_index("_k")
    keys = comp_ledger_keys(comp)
    missing = set(keys) - set(idx.index)
    assert not missing, f"{table}: composition rows without a ledger row: {sorted(missing)[:3]}"
    assert len(set(keys)) == len(comp) == len(idx), f"{table}: ledger/composition row mismatch"
    out = comp.copy()
    out["row_id"] = range(len(comp)) if comp.index.equals(pd.RangeIndex(len(comp))) else comp.index
    out["provenance_class"] = idx.loc[keys, "provenance_class"].to_numpy()
    out["evidence"] = idx.loc[keys, "evidence"].to_numpy()
    mm = mismatches[mismatches["table"] == table]
    kind = dict(zip(mm["row_id"].astype(int), mm["kind"]))
    labels, reasons = [], []
    for rid, cls, ev in zip(out["row_id"], out["provenance_class"], out["evidence"]):
        lab = keep_label(cls, ev)
        why = []
        if lab is None:
            why.append(f"provenance:{cls}")
        if int(rid) in kind:
            why.append(f"identity:{kind[int(rid)]}")
            lab = None
        labels.append(lab)
        reasons.append(";".join(why))
    out["label"] = labels
    out["reason"] = reasons
    return out


def clean_systems(ann: pd.DataFrame) -> set[str]:
    """Fuel ids with at least one clean, non-lumped composition row."""
    lump = ann["is_lumped"].astype(str).str.lower() == "true"
    ok = ann[ann["label"].notna() & ~lump]
    return set(ok["system_id"].astype(str).str.replace(r"^fuel:", "", regex=True))


def build_v4(inv: pd.DataFrame, comp_ann: pd.DataFrame, inv_ledger: pd.DataFrame):
    """Return (inv_v4, comp_v4, inv_quarantined, comp_quarantined, dropped_ids) for one table pair."""
    inv_ids = set(inv["fuel_id"].astype(str))
    keep_ids = clean_systems(comp_ann) & inv_ids
    dropped = sorted(inv_ids - keep_ids)
    fid_of_row = comp_ann["system_id"].astype(str).str.replace(r"^fuel:", "", regex=True)
    in_keep = fid_of_row.isin(keep_ids)

    comp_keep = comp_ann[comp_ann["label"].notna() & in_keep].copy()
    comp_keep["provenance"] = comp_keep["label"]
    comp_keep["verified"] = comp_keep["label"] == "raw_file"
    comp_keep = comp_keep.drop(columns=["label", "reason", "provenance_class", "evidence", "row_id"])

    comp_q = comp_ann[~(comp_ann["label"].notna() & in_keep)].copy()
    comp_q["reason"] = [r if r else NO_CLEAN for r in comp_q["reason"]]
    comp_q = comp_q.drop(columns=["label"])

    led = inv_ledger.copy()
    labs = [keep_label(c, e) for c, e in zip(led["provenance_class"], led["evidence"])]
    led["label"] = labs
    led["_keep"] = [l is not None and f in keep_ids for l, f in zip(labs, led["fuel_id"].astype(str))]
    inv_q = led[~led["_keep"]].copy()
    inv_q["reason"] = [f"provenance:{c}" if l is None else NO_CLEAN
                       for c, l in zip(inv_q["provenance_class"], inv_q["label"])]
    inv_q = inv_q.drop(columns=["label", "_keep"])

    # wide v4 inventory
    kept_led = led[led["_keep"]]
    by_fuel = {f: g for f, g in kept_led.groupby(kept_led["fuel_id"].astype(str))}
    n_q = comp_ann[in_keep & comp_ann["label"].isna()].groupby(fid_of_row[in_keep & comp_ann["label"].isna()]).size()
    ledger_value_fields = set(led.loc[~led["field"].str.startswith("claim:"), "field"])
    rows = []
    for _, r in inv[inv["fuel_id"].astype(str).isin(keep_ids)].iterrows():
        fid = str(r["fuel_id"])
        g = by_fuel.get(fid)
        kept_fields = {} if g is None else {f: l for f, l in zip(g["field"], g["label"])}
        out = r.copy()
        for f in ledger_value_fields:
            if f in out.index and f not in kept_fields:
                out[f] = float("nan")
        claims = [f[len("claim:"):] for f in kept_fields if f.startswith("claim:")]
        if "available_properties" in out.index:
            out["available_properties"] = ";".join(claims)
        vals = {f: l for f, l in kept_fields.items() if not f.startswith("claim:")}
        out["provenance"] = json.dumps(vals, sort_keys=True)
        unv = sorted(f for f, l in vals.items() if l != "raw_file")
        out["unverified_fields"] = ";".join(unv)
        out["verified"] = not unv
        out["n_composition_rows_quarantined"] = int(n_q.get(fid, 0))
        rows.append(out)
    inv_v4 = pd.DataFrame(rows).reset_index(drop=True)
    return inv_v4, comp_keep.reset_index(drop=True), inv_q, comp_q, dropped


def run(root: Path = ROOT, write: bool = True) -> dict:
    ledger = pd.read_csv(root / LEDGER, dtype={"value": str})
    mism = pd.read_csv(root / MISMATCH)
    inv_qs, comp_qs, summary = [], [], []
    v4 = None
    for inv_p, comp_p in PAIRS:
        inv = pd.read_csv(root / inv_p)
        comp = pd.read_csv(root / comp_p)
        is_frac = ledger["field"].str.startswith("fraction:")
        comp_led = ledger[(ledger["table"] == comp_p) & is_frac]
        inv_led = ledger[(ledger["table"] == inv_p) & ~is_frac]
        ann = split_composition(comp, comp_led, mism, comp_p)
        res = build_v4(inv, ann, inv_led)
        inv_v4, comp_v4, inv_q, comp_q, dropped = res
        comp_q.insert(0, "table", comp_p)
        inv_qs.append(inv_q)
        comp_qs.append(comp_q)
        summary.append({"inventory": inv_p, "compositions": comp_p, "comp_rows": len(comp),
                        "comp_kept": len(comp_v4), "comp_quarantined": len(comp_q),
                        "inv_ledger_rows": len(inv_led), "inv_kept": int(len(inv_led) - len(inv_q)),
                        "inv_quarantined": len(inv_q), "fuels": len(inv), "fuels_kept": len(inv_v4),
                        "fuels_dropped": len(dropped)})
        if inv_p == V3_INV:
            v4 = (inv_v4, comp_v4, dropped, inv, comp_led, inv_led)
    inv_v4, comp_v4, dropped, inv3, comp_led3, inv_led3 = v4
    inv_quar = pd.concat(inv_qs, ignore_index=True)
    comp_quar = pd.concat(comp_qs, ignore_index=True)
    result = {"inv_v4": inv_v4, "comp_v4": comp_v4, "inv_quar": inv_quar, "comp_quar": comp_quar,
              "dropped": dropped, "summary": pd.DataFrame(summary), "inv_v3": inv3}
    if write:
        _write(root, result)
    return result


def _write(root: Path, res: dict) -> None:
    (root / QDIR).mkdir(parents=True, exist_ok=True)
    (root / VDIR).mkdir(parents=True, exist_ok=True)
    res["inv_quar"].to_csv(root / QDIR / "inventory_quarantined.csv", index=False)
    res["comp_quar"].to_csv(root / QDIR / "compositions_quarantined.csv", index=False)
    inv3 = res["inv_v3"]
    fam = dict(zip(inv3["fuel_id"].astype(str), inv3["fuel_family"]))
    pd.DataFrame({"fuel_id": res["dropped"], "fuel_family": [fam[f] for f in res["dropped"]],
                  "reason": NO_CLEAN}).to_csv(root / QDIR / "dropped_fuels.csv", index=False)
    res["inv_v4"].to_csv(root / VDIR / "inventory_v4.csv", index=False)
    res["comp_v4"].to_csv(root / VDIR / "compositions_v4.csv", index=False)
    _update_missing(root, res["dropped"], fam)


def _update_missing(root: Path, dropped: list[str], fam: dict) -> None:
    p = root / "data" / "MISSING.md"
    text = p.read_text(encoding="utf-8")
    by_fam: dict[str, list[str]] = {}
    for f in dropped:
        by_fam.setdefault(str(fam[f]), []).append(f)
    lst = "; ".join(f"{k} ({len(v)}): {', '.join(v)}" for k, v in sorted(by_fam.items()))
    block = (f"{MARK_BEGIN}\n"
             f"| 2026-11-10 | v4 Task 3: fuels with no clean composition ({len(dropped)} fuels dropped from "
             f"`data/real_fuels/v4/`) | A raw-file, cited or reconstructed (R10) molecule list for each fuel; "
             f"the v2/v3 compositions of these fuels are target-derived, placeholder or identity-mismatched "
             f"(see `results/v4/P1/AUDIT.md`) | Rows moved to `data/quarantine/`; no replacement data. "
             f"Dropped, by family: {lst}. List also in `data/quarantine/dropped_fuels.csv`. |\n"
             f"{MARK_END}")
    if MARK_BEGIN in text:
        a = text.index(MARK_BEGIN)
        b = text.index(MARK_END) + len(MARK_END)
        text = text[:a] + block + text[b:]
    else:
        text = text.rstrip("\n") + "\n" + block + "\n"
    p.write_text(text, encoding="utf-8")


def main() -> int:
    res = run()
    print(res["summary"].to_string(index=False))
    for lab in ("raw_file", "cited_typed", "reconstructed"):
        print(lab, "kept comp rows:", int((res["comp_v4"]["provenance"] == lab).sum()))
    print("dropped fuels:", len(res["dropped"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
