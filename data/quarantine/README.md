# data/quarantine

Written by `scripts/v4/quarantine.py` (v4 Task 3) from `results/v4/P1/provenance_ledger.csv` (Task 1) and `results/v4/P1/identity_mismatches.csv` (Task 2). Do not hand-edit; rerun the script.

| file | content |
|---|---|
| `inventory_quarantined.csv` | one row per ledger value/claim row removed from v4, for every audited inventory table (`table, fuel_id, field, value, provenance_class, evidence, reason`) |
| `compositions_quarantined.csv` | full composition rows removed from v4 plus `table, row_id, provenance_class, evidence, reason, also_in`; rows duplicated across compositions.csv / compositions_measured.csv / compositions_v3.csv appear once (v3 wins; `also_in` lists every table) |
| `dropped_fuels.csv` | v3 fuels not in v4 (92 of 153: 84 with no clean non-lump composition row, 8 whose kept fractions sum to < 0.95), with the reason and kept fraction; also listed in `data/MISSING.md` |

`reason` is `provenance:<class>` (class not `raw_file`/`cited_typed`, nor a CRC 647 `reconstructed:` row), `identity:<kind>` (stored InChIKey does not belong to the SMILES), or `fuel_dropped:no_clean_composition` / `fuel_dropped:incomplete_composition` (the row is otherwise clean, but its fuel has no clean composition or, after quarantine, keeps less than `MIN_KEPT_FRACTION` = 0.95 of its composition, lumped share included; R15).

Kept in `data/real_fuels/v4/` with a `provenance` label: `raw_file` (verified=True), `cited_typed` (typed from a named paper, no local raw file; verified=False; excluded from headline metrics until a raw file confirms it) and `reconstructed` (CRC 647 compositions built from measured bulk aromatics/olefins and D86; verified=False).

Quarantined data must not be used for training, validation or any reported number. Retracted claims are in `docs/claims_register.md` (CLM-20 to CLM-29); pinned tests are in `tests/quarantine/` (excluded by `pytest.ini`).

v4 inventory: 61 fuels, all jet (57 CRC 647 reconstructed, NIST_RP1 and 3 others). `unresolved_fraction` = 1 - kept composition fraction (never NaN); `verified` is True only when every kept value and composition row is `raw_file`, no composition row was quarantined and the kept sum is >= 0.99.
