#!/usr/bin/env python
"""v4 Task 1: provenance audit of the v2/v3 real-fuel tables.

Runs ``fuelblend.data.provenance.classify_fuel_rows`` on every
inventory/composition table pair under ``data/real_fuels`` (including ``v3/``)
with the two ingest scripts as static evidence, and writes

* ``results/v4/P1/provenance_ledger.csv`` -- one row per value
  (``table, fuel_id, field, value, provenance_class, evidence``), sorted.
* ``results/v4/P1/AUDIT.md`` -- counts per class per family, the plan claims
  confirmed / not confirmed from code+data, extra findings, and the v2/v3
  deliverables that read non-``raw_file`` rows.

Read-only with respect to v2/v3 data, scripts and results.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.data.provenance import (  # noqa: E402
    CLASSES, DEFAULT_RAW_ROOTS, INGEST_SCRIPTS, LEDGER_COLUMNS, PROPERTY_VAR_NAMES,
    classify_fuel_rows, load_source_code,
)

RF = "data/real_fuels"
PAIRS = [  # (inventory, compositions) -- every inventory*/compositions* file and data/real_fuels/v3/*
    (f"{RF}/inventory.csv", f"{RF}/compositions.csv"),
    (f"{RF}/inventory_measured.csv", f"{RF}/compositions_measured.csv"),
    (f"{RF}/inventory_reconstructed.csv", f"{RF}/compositions_reconstructed.csv"),
    (f"{RF}/v3/inventory_v3.csv", f"{RF}/v3/compositions_v3.csv"),
]
OUT_DIR = ROOT / "results" / "v4" / "P1"
KEEP = ("raw_file", "cited_typed")


def build_ledger() -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    src = load_source_code(ROOT)
    assert set(src) == set(INGEST_SCRIPTS), f"missing ingest scripts: {set(INGEST_SCRIPTS) - set(src)}"
    parts, invs = [], {}
    for inv_p, comp_p in PAIRS:
        inv = pd.read_csv(ROOT / inv_p)
        comp = pd.read_csv(ROOT / comp_p)
        invs[inv_p] = inv
        led = classify_fuel_rows(inv, comp, src)
        is_frac = led["field"].str.startswith("fraction:")
        led.insert(0, "table", [comp_p if f else inv_p for f in is_frac])
        # coverage: every inventory row and every composition row produced ledger rows
        missing = set(inv["fuel_id"].astype(str)) - set(led.loc[~is_frac, "fuel_id"])
        assert not missing, f"{inv_p}: fuels without a ledger row: {sorted(missing)[:5]}"
        assert int(is_frac.sum()) == len(comp), f"{comp_p}: {int(is_frac.sum())} ledger rows vs {len(comp)}"
        parts.append(led)
        print(f"  {inv_p}: {len(inv)} fuels, {comp_p}: {len(comp)} rows -> {len(led)} ledger rows")
    ledger = pd.concat(parts, ignore_index=True)
    assert not ledger.duplicated(["table", "fuel_id", "field"]).any()
    assert set(ledger["provenance_class"]) <= set(CLASSES)
    ledger = ledger.sort_values(["table", "fuel_id", "field"], kind="mergesort").reset_index(drop=True)
    return ledger[["table"] + LEDGER_COLUMNS], invs


def family_map(invs: dict[str, pd.DataFrame]) -> dict[str, str]:
    fam = {}
    for inv in invs.values():
        for fid, f in zip(inv["fuel_id"].astype(str), inv["fuel_family"].astype(str)):
            fam.setdefault(fid, f)
    return fam


def inchikey_mismatches(paths: list[str]) -> pd.DataFrame:
    """Rows whose stored InChIKey is not the InChIKey of the stored SMILES (RDKit)."""
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:  # pragma: no cover
        return pd.DataFrame()
    rows, cache = [], {}
    for p in paths:
        c = pd.read_csv(ROOT / p)
        for sid, ik, smi in zip(c["system_id"], c["inchikey"].astype(str), c["smiles"]):
            if not isinstance(smi, str) or not smi.strip() or not re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", ik):
                continue
            if smi not in cache:
                m = Chem.MolFromSmiles(smi)
                cache[smi] = Chem.MolToInchiKey(m) if m is not None else None
            if cache[smi] and cache[smi] != ik:
                rows.append({"table": p, "system_id": sid, "smiles": smi, "stored_inchikey": ik,
                             "inchikey_of_smiles": cache[smi]})
    return pd.DataFrame(rows)


def consumers() -> dict[str, list[str]]:
    """Scripts / modules that read each table (grep for the file name), excluding the writers."""
    names = {Path(p).name for pair in PAIRS for p in pair}
    out: dict[str, list[str]] = {n: [] for n in sorted(names)}
    for base in ("scripts", "src"):
        for py in sorted((ROOT / base).rglob("*.py")):
            rel = py.relative_to(ROOT).as_posix()
            if rel.startswith("scripts/v4/") or rel.endswith("data/provenance.py"):
                continue
            text = py.read_text(encoding="utf-8", errors="replace")
            for n in names:
                if re.search(r"[\"'/]" + re.escape(n) + r"[\"']", text) or (
                        n.endswith("_v3.csv") and n in text):
                    tag = rel + (" (writer)" if rel in INGEST_SCRIPTS and _writes(text, n) else "")
                    out[n].append(tag)
    return out


def _writes(text: str, name: str) -> bool:
    return bool(re.search(re.escape(name) + r"[\"']\s*\)?,?\s*index\s*=\s*False", text)) or \
        bool(re.search(r"to_csv\([^)]*" + re.escape(name), text))


def md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def claim_checks(L: pd.DataFrame, invs: dict[str, pd.DataFrame]) -> list[tuple[str, str, str]]:
    """Plan claims, each re-derived from the ledger (code+data evidence)."""
    v3 = L[L["table"].str.contains("/v3/")]
    v2 = L[L["table"].isin([f"{RF}/inventory.csv", f"{RF}/compositions.csv"])]
    out = []

    def add(claim, ok, detail):
        out.append((claim, "CONFIRMED" if ok else "NOT CONFIRMED", detail))

    for name, led in (("v3", v3), ("v2", v2)):
        f = led[led["fuel_id"].str.fullmatch(r"CRC_FACE_[A-J]") & led["field"].str.startswith("fraction:")]
        n_t = f.loc[f["provenance_class"] == "target_derived_composition", "fuel_id"].nunique()
        tol = f[f["field"].str.contains("toluene|UAEPNZWRGJTJPN|YXFVVABEGXRONW", regex=True)]
        add(f"FACE gasoline A-J compositions computed from each fuel's own RON/MON/density ({name})",
            n_t == 10 and tol["evidence"].str.contains("dens - 680").all(),
            f"{n_t}/10 systems target_derived_composition; toluene rows cite `(dens - 680.0) / 200.0`")
    labels = sorted({str(t) for inv in invs.values()
                     for fid, t in zip(inv["fuel_id"].astype(str), inv["composition_type"])
                     if re.fullmatch(r"CRC_FACE_[A-J]", fid)})
    add("...while labelled composition_type=measured_dha_gcxgc", labels == ["measured_dha_gcxgc"],
        f"composition_type of FACE A-J rows in all inventories: {labels} (label not used by the classifier)")
    e = v3[v3["fuel_id"].str.fullmatch(r"CRC_FACE_[A-J]_E\d+") & (v3["field"] == "dvpe_kpa")]
    add("Ethanol-blend DVPE = base + 7.0 (<=E20) else +3.0",
        len(e) > 0 and (e["provenance_class"] == "formula_derived").all() and e["evidence"].str.contains(
            r"\+ 7.0|\+ 3.0").all(), f"{(e['provenance_class'] == 'formula_derived').sum()}/{len(e)} rows formula_derived")
    h = v3[v3["fuel_id"].str.fullmatch(r"CRC_FACE_[A-J]") & (v3["field"] == "heat_of_combustion_MJ_kg")]
    add("Every FACE gasoline heat of combustion is the constant 43.5",
        len(h) == 10 and (h["value"] == "43.5").all() and (h["provenance_class"] == "placeholder_constant").all(),
        f"{len(h)} rows, values {sorted(h['value'].unique())}, classes {sorted(h['provenance_class'].unique())}; "
        "E-blends use round(43.5*w_base + 26.8*w_eth, 1) (formula_derived)")
    d = v3[v3["fuel_id"].str.startswith("CRC_FACE_Diesel_") & (v3["field"] == "flash_point_C")]
    add("FACE diesel flash = 52 + (rho-800)*0.25",
        len(d) == 9 and (d["provenance_class"] == "formula_derived").all(), f"{len(d)} rows formula_derived")
    for name, led in (("v3", v3), ("v2", v2)):
        dc = led[led["fuel_id"].str.startswith("CRC_FACE_Diesel_") & led["field"].str.startswith("fraction:")]
        n_t = dc.loc[dc["provenance_class"] == "target_derived_composition", "fuel_id"].nunique()
        resolved = dc[~dc["field"].str.contains("LUMP_")]
        add(f"FACE diesel compositions back-solved from DCN and density ({name})",
            n_t == 9 and (dc["provenance_class"] == "target_derived_composition").all()
            and resolved["evidence"].str.contains("dcn").all() and resolved["evidence"].str.contains("dens").all(),
            f"{n_t}/9 systems target_derived_composition; every resolved row's evidence cites dcn and dens")
    fr = v3[(v3["field"] == "freezing_point_C") & v3["evidence"].str.contains("cloud_point_C - 3")]
    add("freezing = cloud - 3", len(fr) > 0, f"{len(fr)} fuels: {', '.join(sorted(fr['fuel_id']))}")
    po = v3[(v3["field"] == "pour_point_C") & v3["evidence"].str.contains("freezing_point_C - 4")]
    add("pour = freeze - 4", len(po) > 0, f"{len(po)} fuels: {', '.join(sorted(po['fuel_id']))}")
    na = v3[v3["fuel_id"].str.fullmatch(r"NASA_LFP[35]") & (v3["field"] == "flash_point_C")]
    add("NASA flash is the constant 58.0", len(na) > 0 and (na["value"] == "58").all()
        and (na["provenance_class"] == "placeholder_constant").all(),
        f"{len(na)} fuels ({', '.join(sorted(na['fuel_id']))}), value 58.0, placeholder_constant")
    dn = v3[(v3["field"] == "density_kg_m3") & v3["evidence"].str.contains('804.0 if "RP1"', regex=False)]
    add('RP-1/POSF densities are 804.0 if "RP1" else 795.0',
        len(dn) > 0 and (dn["provenance_class"] == "placeholder_constant").all(),
        f"{len(dn)} carried-over jet fuels: {(dn['value'] == '804').sum()} x 804, {(dn['value'] == '795').sum()} x 795")
    po5 = v3[v3["fuel_id"].str.fullmatch(r"POSF_(10264|10289|11498|12344|12345)") & v3["field"].str.startswith(
        "fraction:")]
    add("All five POSF additional fuels get the same 48/32/20 split",
        po5["fuel_id"].nunique() == 5 and (po5["provenance_class"] == "placeholder_constant").all(),
        f"{po5['fuel_id'].nunique()} systems, identical vector, placeholder_constant")
    c = L[L["table"].str.endswith("inventory_reconstructed.csv") & (L["field"] == "claim:freezing_point_C")]
    add("CRC 647 measured freezing points exist in the raw PDF table (crc647_fuel_properties.csv)",
        len(c) == 57 and c["evidence"].str.contains("raw exists: .*crc647_fuel_properties.csv").all(),
        f"{len(c)}/57 `claim:freezing_point_C` rows name the raw value (`raw exists:`); they stay `unknown` because "
        "inventory_reconstructed.csv stores no value for it")
    cv = v3[v3["fuel_id"].str.startswith("CRC647_") & (v3["field"] == "freezing_point_C")]
    add("...but inventory_v3.csv stores one constant (-48.0) for all 57 CRC 647 fuels instead",
        len(cv) == 57 and (cv["provenance_class"] == "placeholder_constant").all(),
        f"{len(cv)} rows placeholder_constant; evidence lists the contradicting raw value per fuel")
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Building provenance ledger ...")
    L, invs = build_ledger()
    L.to_csv(OUT_DIR / "provenance_ledger.csv", index=False, lineterminator="\n")
    print(f"Wrote {OUT_DIR / 'provenance_ledger.csv'} ({len(L)} rows)")

    fam = family_map(invs)
    L2 = L.assign(family=L["fuel_id"].map(fam).fillna("?"),
                  kind=L["field"].map(lambda f: "composition" if f.startswith("fraction:")
                                      else ("claim" if f.startswith("claim:") else "value")))
    ct = (L2.groupby(["table", "family", "provenance_class"]).size().unstack(fill_value=0)
          .reindex(columns=list(CLASSES), fill_value=0).reset_index())
    ct_kind = (L2.groupby(["table", "kind", "provenance_class"]).size().unstack(fill_value=0)
               .reindex(columns=list(CLASSES), fill_value=0).reset_index())

    # systems (whole composition vectors) per class, per family, v3
    comp_v3 = L2[L2["table"].str.endswith("compositions_v3.csv")]
    sys_cls = comp_v3.groupby(["family", "fuel_id"])["provenance_class"].agg(
        lambda s: s.value_counts().index[0]).reset_index()
    sys_ct = sys_cls.groupby(["family", "provenance_class"]).size().unstack(fill_value=0).reindex(
        columns=list(CLASSES), fill_value=0).reset_index()

    checks = claim_checks(L, invs)
    mism = inchikey_mismatches([p for _, p in PAIRS])
    cons = consumers()

    v3 = L2[L2["table"].str.contains("/v3/")]
    conflicts = v3[v3["evidence"].str.contains("raw CONFLICT") & (v3["kind"] == "value")]
    ron_conf = conflicts[conflicts["field"].isin(["ron", "mon"])]
    crc_conf = conflicts[conflicts["fuel_id"].str.startswith("CRC647_")]
    coinc = v3[v3["fuel_id"].str.startswith("CRC647_") & v3["evidence"].str.contains("coincides with the raw value")]
    basis = L2[L2["evidence"].str.contains("BASIS MISMATCH")]
    nres = L2[L2["evidence"].str.contains("resolved components but the composition table")]
    identical = L2[L2["evidence"].str.contains("identical fraction vector") & L2["table"].str.contains("/v3/")]

    non_kept = L2[~L2["provenance_class"].isin(KEEP)]
    per_table_nonkeep = non_kept.groupby("table").size().to_dict()

    md = []
    md.append("# v4 P1 — Provenance audit of the v2/v3 real-fuel tables\n")
    md.append("Generated by `scripts/v4/audit_provenance.py` (classifier: `src/fuelblend/data/provenance.py`). "
              "Ledger: `results/v4/P1/provenance_ledger.csv` "
              f"({len(L)} rows, sorted by `table, fuel_id, field`; deterministic).\n")
    md.append("## What was audited\n")
    md.append("Every row of: " + ", ".join(f"`{p}`" for pair in PAIRS for p in pair) + ".\n")
    md.append("Ledger rows per table: one per numeric inventory value (`ron, mon, dvpe_kpa, density_kg_m3, "
              "viscosity_kin_cSt, flash_point_C, freezing_point_C, cloud_point_C, pour_point_C, dcn, "
              "heat_of_combustion_MJ_kg, number_of_resolved_components, unresolved_fraction, ethanol_vol_pct, "
              "ethanol_mass_pct, oxygenate_pct`; NaN cells are *no value* and are skipped), one `claim:<prop>` row "
              "per property listed in `available_properties` (the v2 inventories carry no property values, only "
              "these claims), and one `fraction:<inchikey>#<n>` row per composition row (`n` = ordinal among rows "
              "of the same system with the same InChIKey). The script asserts every inventory fuel and every "
              "composition row is covered.\n")
    md.append("## Method (labels are never trusted)\n")
    md.append("* **Static (AST):** both ingest scripts (`" + "`, `".join(INGEST_SCRIPTS) + "`) are parsed and the "
              "code that appends each row dict is *replayed* in a whitelisted evaluator (literal fuel tables, "
              "arithmetic, `min/max/round`; no I/O). For each value the ledger records the expression "
              "(`script:Lnnn`), whether it is a literal shared by many fuels, a per-fuel literal tuple entry, or "
              "arithmetic, and which per-fuel inputs it depends on. Per-fuel loop variables named "
              f"{', '.join('`' + k + '`' for k in sorted(PROPERTY_VAR_NAMES))} are treated as *property* inputs.\n"
              "* **Data:** exact relations are tested on the stored values (`freeze == cloud-3`, `pour == freeze-4`, "
              "`flash == round(52+(rho-800)*0.25,1)`, E-blend `dvpe == base + 7.0/3.0`, E-blend volume-weighted "
              "density, the FACE gasoline/diesel allocation replayed from RON/MON/density/DCN with the scripts' "
              "clip + normalise + round, identical fraction vectors across systems, constants shared by fuels "
              "built by the same dict literal, claimed vs actual component counts).\n"
              "* **Raw files:** `raw_file` requires a per-fuel, per-property match in "
              + ", ".join(f"`{r}`" for r in DEFAULT_RAW_ROOTS) + ": CRC 647 (`pour_point/crc647_fuel_properties.csv`, "
              "fuel_code = numeric suffix), PRF+ethanol (`octane/face_gasoline_set_I.csv`), FACE+ethanol "
              "(`octane/face_gasoline_set_II.csv`), TIGAS / Co-Optima 4-comp / NREL BOB-01 "
              "(`cooptima/fuel_engine_co_optimization.csv`), NIST RP-1 molecule list (`rp1/rp1_blends_composition.csv`), "
              "surrogate recipes (`surrogates/fuel_surrogates.csv`, SMILES identity + common renormalisation factor "
              "+ basis), POSF GC classes (`gc_composition/fuel_gc_compositions.csv`). There is no global "
              "\"value appears somewhere\" search.\n")
    md.append("## Class definitions as applied\n")
    md.append("* `raw_file` — matches the raw value for this fuel and property (or the raw molecule list with the same "
              "SMILES and basis, up to one common renormalisation factor). A `claim:` row (value `listed`) is "
              "`raw_file` **only** when the same table stores a value for that property and that stored value is "
              "itself `raw_file`; when a raw value exists but the table stores nothing (v2 inventories) or stores a "
              "placeholder/formula value, the claim is `unknown` and its evidence reads `raw exists: <path> -> "
              "<value>`. Claims carry no data; a keep-filter should use the value rows.\n"
              "* Two independent detection routes: the code replay, and pure data rules that need no code (the "
              "FACE gasoline/diesel allocation replayed from the stored RON/MON/density/DCN, E-blend fractions "
              "proportional to a target-derived base, `freeze == cloud-3`, `pour == freeze-4`, the diesel flash "
              "formula, E-blend DVPE/density, and one property value on >= 3 fuels of a family -> "
              "`placeholder_constant`). `classify_fuel_rows(inv, comp, {})` still flags the FACE compositions, the "
              "43.5 HoC and the CRC 647 constants (tested).\n"
              "* Property inputs of the replay are found three ways: the variable is written to exactly one property "
              "column of a row dict in its loop, its name (list above), or a subscripted read of a property key "
              "(`face_map[b][\"dens\"]`).\n"
              "* `cited_typed` — narrow: a per-fuel literal typed in the script (or copied from a `data/processed` "
              "table that carries a per-row literature reference) with **no** formula, no sharing with other fuels "
              "and no conflict with a mapped raw file. Not verified.\n"
              "* `formula_derived` — reproduced by a formula / replayed expression, including vol→mass conversions with "
              "assumed densities, counts derived from the composition table, and model reconstructions.\n"
              "* `target_derived_composition` — a fraction that depends (through the replayed code) on the same fuel's "
              "RON/MON/density/DCN/...; a composition vector is one object, so all rows of such a system "
              "(including its lump/normaliser) carry the class.\n"
              "* `placeholder_constant` — one literal assigned to >=2 fuels by the same dict literal (or one fraction "
              "vector shared by >=2 systems), even where it coincides with a raw value.\n"
              "* `unknown` — not shown to be any of the above; includes typed values that **conflict** with the raw "
              "file the script names as its source, mole fractions stored as `basis=mass`, and structural zeros "
              "(`0.0` ethanol / unresolved fraction literals).\n")
    md.append("### CRC 647 reconstructed compositions\n")
    md.append("They are mapped to **`formula_derived`** with evidence starting `reconstructed:` "
              "(`scripts/v2/ingest_real_fuels_v2.py` builds them with `bulk_items`/`to_vector` from fuels.csv bulk "
              "aromatics/olefins/diaromatics vol% and D86 points, with a fixed saturate split 0.45/0.35/0.20). They "
              "are not measured molecule lists, but they are not computed from the fuels' target properties "
              "either (note: D86 is also listed as an available property, so a D86 model on them is partly circular). "
              "**Consequence for Task 3:** a rule that keeps only `raw_file`/`cited_typed` would quarantine them; if "
              "they should stay usable with a `reconstructed` label, Task 3 must whitelist "
              "`formula_derived` rows whose evidence starts with `reconstructed:`. Left for the controller to decide.\n")

    md.append("## Counts per class per family (all ledger rows)\n")
    md.append(md_table(ct) + "\n")
    md.append("## Counts per class by row kind\n")
    md.append(md_table(ct_kind) + "\n")
    md.append("## v3 composition systems by majority class\n")
    md.append(md_table(sys_ct) + "\n")

    md.append("## Plan claims checked against code + data\n")
    md.append(md_table(pd.DataFrame(checks, columns=["claim", "status", "ledger evidence"])) + "\n")
    n_unconf = sum(1 for _, s, _ in checks if s != "CONFIRMED")
    md.append(f"{len(checks) - n_unconf}/{len(checks)} claims confirmed." + (
        " Unconfirmed claims are listed with status NOT CONFIRMED above." if n_unconf else "") + "\n")

    md.append("## Additional findings (not in the plan)\n")
    md.append(f"1. **FACE / PRF ethanol-blend RON/MON vs raw.** {len(ron_conf)} v3 RON/MON values conflict with "
              "the raw octane file the v3 `source` column names (e.g. FACE A+10 is 92/88 in "
              "`face_gasoline_set_II.csv` but 91.7/86.3 in v3; FACE I+10 RON is 79.9 raw vs 100.8 v3). PRF+ethanol "
              "RON/MON do match `face_gasoline_set_I.csv` (`raw_file`). Conflicting values are `unknown`.\n"
              f"2. **CRC 647 values in v3 are constants.** {len(crc_conf)} CRC 647 property values in inventory_v3.csv "
              "conflict with the raw PDF table (density 798.0, viscosity 1.55, flash 45.0, freezing -48.0, pour -52.0, "
              f"HoC 43.2 for all 57 fuels); {len(coinc)} of these constants happen to equal the raw value for "
              f"that fuel ({', '.join(sorted(coinc['fuel_id'] + '.' + coinc['field']))}) and are still "
              "placeholders.\n"
              "3. **POSF-10325/13511/13516 GC compositions** are hydrocarbon-class percentages (including sub-totals and "
              "`*_total` rows, e.g. POSF-10325 sums to 152.28 %) rescaled and mapped to representative molecules by a "
              "substring rule (every `*paraffin*` class, including cycloparaffins, becomes n-dodecane) — "
              "`formula_derived`, not a measured molecule list.\n"
              f"4. **Surrogate basis mislabel.** {len(basis)} composition rows store raw *mole* fractions from "
              "`fuel_surrogates.csv` with `basis=mass` (Dooley, HyChem, Kim/Violi, Violi JP-8) — `unknown`. Seven "
              "HyChem rows have no SMILES in the raw file and are stored under the key `SURR_nan`.\n"
              f"5. **Component counts.** {nres['fuel_id'].nunique()} fuels claim far more resolved components than "
              "their composition rows (FACE gasolines claim 78/79 with 4-5 rows; FACE diesels 55 with 4; POSF "
              "additional 24 with 3; NASA 8 with 8 n-alkanes only).\n"
              f"6. **Shared recipes.** In v3, {identical['fuel_id'].nunique()} systems share an identical fraction "
              "vector with another system (POSF additional x5, NASA LFP3/LFP5, Co-Optima 4-comp = TIGAS = NREL BOB-01).\n"
              f"7. **InChIKey does not match SMILES** in {len(mism)} composition rows (RDKit check; "
              f"{mism['stored_inchikey'].nunique() if len(mism) else 0} distinct wrong keys). "
              "Examples: the NASA n-alkane rows (v2 and v3) carry the isooctane key on n-C13, the n-heptane key on "
              "n-C16 and the n-hexadecane key on n-C10; v3 toluene uses `UAEPNZWRGJTJPN` (v2 uses the correct "
              "`YXFVVABEGXRONW`). Any join on InChIKey (e.g. to `data/unified/compounds.csv`) picks up the wrong "
              "molecule for these rows.\n"
              "8. **TRF-60..90** contain only isooctane + n-heptane (no toluene) and set RON = MON = isooctane vol%.\n")
    if len(mism):
        top = mism.groupby(["stored_inchikey", "smiles", "inchikey_of_smiles"]).size().reset_index(name="rows")
        md.append("InChIKey/SMILES mismatches (all tables):\n\n" + md_table(top.sort_values("rows", ascending=False))
                  + "\n")

    md.append("## v2/v3 deliverables that consume non-`raw_file` rows\n")
    md.append("Rows that are neither `raw_file` nor `cited_typed` per table: " + ", ".join(
        f"`{k}`: {v}" for k, v in sorted(per_table_nonkeep.items())) + ". `raw_file` alone: " + ", ".join(
        f"`{k}`: {int((L2[L2['table'] == k]['provenance_class'] == 'raw_file').sum())}/{int((L2['table'] == k).sum())}"
        for k in sorted(L2["table"].unique())) + ".\n")
    md.append("Every table below contains non-`raw_file` rows, so every reader consumes them unless it filters by "
              "fuel. Readers found by grepping `scripts/` and `src/` for the file name:\n")
    for n, users in cons.items():
        md.append(f"* `{n}`: " + (", ".join(f"`{u}`" for u in users) if users else "(no reader found)"))
    md.append("")
    md.append("## Limitations\n")
    md.append("* `cited_typed` values are typed literals that were **not** verified; only the raw files listed above "
              "are mapped. FACE base gasoline RON/MON/density/DVPE, FACE diesel DCN/density/viscosity/cloud/pour, "
              "POSF additional properties and the diesel surrogate properties have no raw file in the repo.\n"
              "* A property bound under a new name that is never written to a property column nor read by a "
              "property key would not taint the replay; the pure data rules still catch the FACE/diesel cases on "
              "the v3 tables, but not on the v2 tables, which store no property values.\n"
              "* The v2 inventories hold no property values; for them only `claim:` rows and metadata are audited.\n"
              "* The replay covers the dict literals in the two ingest scripts' `main()`; rows produced elsewhere "
              "would show as `unknown`.\n")
    (OUT_DIR / "AUDIT.md").write_text("\n".join(md) + "\n", encoding="utf-8", newline="\n")
    print(f"Wrote {OUT_DIR / 'AUDIT.md'}")
    for c, s, d in checks:
        print(f"  [{s}] {c}")


if __name__ == "__main__":
    main()
