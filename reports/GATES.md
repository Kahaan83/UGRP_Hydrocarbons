# GATES.md (updated under Rule R1, Phase A)

Thresholds per handoff 17.4:
- **G1**: `n_binary_ternary >= 1000` AND `n_multi >= 100` AND `n_real >= 20`
- **G2**: `(100 <= n_binary_ternary <= 1000)` OR `(5 <= n_real <= 19)`
- **G3**: below G2 thresholds

Under **Rule R1** (`data/processed/model_table.csv`):
- L1 (class totals) and L0 (bulk properties / D86 without GC×GC) fuels are **excluded** (moved to `data/processed/excluded_by_R1.csv`).
- Consequently, uncharacterized `real_fuel` rows count as 0 for model training gates until L2 GC×GC data (e.g. Vozka 2019, Caltech FM99-4) are ingested in Phase B.
- System counts below represent unique chemical systems (`system_id`):

| property | n_multi_systems | n_binary_systems | n_bt_systems | n_real_L2 | gate | Note |
|---|---|---|---|---|---|---|
| `viscosity_dyn_mPas` | 19 | 360 | 379 | 0 | **G2** | Strong binary + ternary coverage; eligible for residual learning. |
| `density_kg_m3` | 26 | 209 | 235 | 0 | **G2** | Strong binary + multi coverage; eligible for residual learning. |
| `flash_point_C` | 1 | 30 | 31 | 0 | **G3** | L3 binaries/ternaries currently limited (31 systems); baseline + physics priority. Needs Phase B USNA/Li data. |
| `viscosity_kin_cSt` | 0 | 1 | 1 | 0 | **G3** | Pure-heavy, sparse mixture rows without pure density. |
| `dcn` | 19 | 19 | 38 | 0 | **G3** | Sparse multi systems; baseline evaluation only. |
| `mon` | 15 | 7 | 22 | 0 | **G3** | Sparse multi systems; baseline evaluation only. |
| `vapor_pressure_kPa` | 4 | 61 | 65 | 0 | **G3** | Used primarily for Antoine refit. |
| `freezing_point_C` | 0 | 0 | 0 | 0 | **G3** | No open mixture freezing data yet; awaiting Phase B (Daridon/Coutinho). |
| `pour_point_C` | 0 | 0 | 0 | 0 | **G3** | No L3/L2 pour point mixture data yet; awaiting Phase B. |
| `smoke_point_mm` | 0 | 0 | 0 | 0 | **G3** | Previous real fuels lacked L2 composition; awaiting Mensch/Boehm TSI in Phase B. |
| `net_heat_content_btu_lb`| 0 | 0 | 0 | 0 | **G3** | Bulk fuels excluded under R1. |
