# FACE gasoline composition recovery report (computed)

Produced by `scripts/external_ingest/resolve_face_gasoline_dha.py`. All counts below are computed from the run.

| Metric | Value |
|---|---|
| Legacy rows read | 144 |
| Rows resolved | 144 |
| Rows unresolved | 0 |
| Surrogate blends resolved to named compounds (Table 4) | 120 |
| FACE + ethanol blends resolved on a volume basis (Table 5) | 24 |

Surrogate rows by number of components: {2: 21, 3: 64, 4: 35}.
FACE base fuels present: FACE A, FACE C, FACE F, FACE G, FACE I, FACE J.

Notes:
- InChIKeys and average molecular weights come from RDKit; densities from `data/unified` (median near 293 K) or the
  CRC handbook values in AI4fuel (see `octane_blending_report.md`).
- `system_id` uses the production helper; FACE rows use `fuel:kaust_jameel2018:<FACE X>+ethanol`.
- FACE rows carry volume fractions only; mass and mole fractions are left empty because the FACE lump has no staged
  molecular weight or density.
- Output: `data/interim/external_staging/canonical_resolved_face_gasoline.csv`. Staging only; production is not modified.
