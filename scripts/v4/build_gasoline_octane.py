"""v4 Task 5: build data/real_fuels/v4/gasoline_octane_blends.csv + results/v4/P2/gasoline_corpus.md."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.data.loaders_octane import load_cooptima_bron, load_face_sets  # noqa: E402

RAW = ROOT / "extras" / "fuelblend_data" / "data" / "raw"
OUT = ROOT / "data" / "real_fuels" / "v4" / "gasoline_octane_blends.csv"
REPORT = ROOT / "results" / "v4" / "P2" / "gasoline_corpus.md"


def main() -> None:
    face = load_face_sets([RAW / "octane" / "face_gasoline_set_I.csv",
                           RAW / "octane" / "face_gasoline_set_II.csv"])
    bron = load_cooptima_bron(RAW / "cooptima" / "bron_data_with_cas.xls")
    df = pd.concat([face, bron], ignore_index=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)

    def ethanol_known(s):
        return any(c["name"] == "Ethanol" and c["frac"] is not None for c in json.loads(s))

    def ethanol_unknown(s):
        return any(c["name"] == "Ethanol" and c["frac"] is None for c in json.loads(s))

    lines = ["# Gasoline octane corpus (v4 Task 5)", "",
             f"Output: `{OUT.relative_to(ROOT).as_posix()}` ({len(df)} rows).  Built by "
             "`scripts/v4/build_gasoline_octane.py`; raw inputs under `extras/fuelblend_data/data/raw/` "
             "(gitignored).", "",
             "## Rows per source", "", "| source | rows | with RON | with MON | composition_resolved | target_ok | target_ok_ron | target_ok_mon |",
             "|---|---|---|---|---|---|---|---|"]
    for src, g in df.groupby(df.source.str.split(":").str[0]):
        lines.append(f"| {src} | {len(g)} | {g.ron.notna().sum()} | {g.mon.notna().sum()} | "
                     f"{int(g.composition_resolved.sum())} | {int(g.target_ok.sum())} | "
                     f"{int(g.target_ok_ron.sum())} | {int(g.target_ok_mon.sum())} |")
    lines.append(f"| **all** | {len(df)} | {df.ron.notna().sum()} | {df.mon.notna().sum()} | "
                 f"{int(df.composition_resolved.sum())} | {int(df.target_ok.sum())} | "
                 f"{int(df.target_ok_ron.sum())} | {int(df.target_ok_mon.sum())} |")
    f_eth = face[face.components.map(ethanol_known)]
    lines += ["", "## Composition coverage", "",
              f"* FACE rows: {len(face)}; with an ethanol volume fraction parsed from the Name suffix: "
              f"{len(f_eth)}; all other FACE rows are an unresolved lump (fraction 1.0) with no ethanol, except the "
              f"{int(face.components.map(ethanol_unknown).sum())} "
              "`TPRF 4+N` rows whose Ethanol cell is present but share unknown: Ethanol kept as a component "
              "(frac=None, basis=unresolved) and in unresolved_members.",
              f"* FACE rows whose hydrocarbon fractions are known: 0 (the CSVs contain component NAMES only, "
              "no fractions; see the loader docstring).",
              f"* FACE rows on a FACE A-J base lump: {int(face.base_fuel.str.startswith('FACE_').sum())}.",
              f"* bRON rows: {len(bron)} from {bron.blend_id.str.split(':').str[1].nunique()} blendstocks; "
              f"the blendstock fraction is a measured volume fraction in every row, the base gasoline "
              f"(distinct base RON values: {sorted(bron.base_ron.dropna().unique().tolist())}) is a lump.",
              "* FACE rows: target_ok False (lump, no fractions).  bRON rows (ruling R17): blend RON is a target "
              "only where the base RON is supplied as an input (target_ok_ron), MON likewise (target_ok_mon); "
              "`target_ok` is their OR.  Ethanol share basis for FACE is `volume_inferred`.",
              "", "## Provenance", "",
              "* `provenance` (ron/mon values): raw_file for all rows.",
              "* `composition_provenance`: raw_file for bRON (numeric vol fractions in file); "
              "unknown for FACE (ethanol % inferred from the Name convention; hydrocarbon ratios absent).",
              "", "## Skipped / not ingested", "",
              "* Nothing skipped inside the two files (every data row yields a row).",
              "* `TPRF 4+N` (N=1..15) not interpreted as ethanol percentages (series index): lump, "
              "plus an unresolved Ethanol component where the Ethanol cell is filled.",
              "* Co-Optima `fuel_engine_co_optimization.csv` Blend_RON1..3/percent columns and "
              "`bob_data.xlsx` are NOT ingested here (outside Task 5 scope).",
              "* v3's typed FACE A-J base RON/MON and compositions (`ingest_phase1_expansion.py`) are NOT "
              "used: 47 of them conflict with the raw file (P1 audit) and no raw GC composition exists.",
              "* Diisobutylene (CAS 107-39-1) carries the InChIKey of the `chemicals` database entry; "
              "the commercial material is an isomer mixture.  'Ketone Mixture' and 'Furan Mixture' "
              "(no CAS) are blendstock lumps.", ""]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
