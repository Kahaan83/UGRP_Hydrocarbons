# Master Claims Register (UGRP Hydrocarbons v2)

> Policy (Ground Rule 2): Any earlier claim that fails a v2 check is marked `RETRACTED` or `REVISED` in this register with the reason, replacing result, and reference to verification files. Old reports are not silently edited.

| Claim ID | Original Claim | Source File / Location | v2 Status | Replacement Finding / Reference | Reason & Evidence |
|---|---|---|---|---|---|
| **CLM-01** | Feasible Jet A 2-component sparse blend (dodecane + tetralin proxy) passes all constraints with 3.4 °C safety margin | `FULL_REPORT.md` line 222 | **RETRACTED** | `results/v2/A1/design_audit.csv`, `docs/progress/v2/A1.md` | Dodecane melts at -9.6 °C; blend crystallizes at -15.2 °C (fails $\le -40$ °C by 24.8 °C). Completely solid at -20 °C (fails kinematic viscosity $\le 8.0$ cSt). Tetralin is an aromatic (ASTM D1319), so 25.5 wt% violates the 25.0 wt% limit. |
| **CLM-02** | Feasible 3-family inverse designs in `results/t6/inverse_designs_three_families.csv` (Jet, Diesel, Gasoline) pass all ASTM specifications | `FULL_REPORT.md` line 221, `results/t6/inverse_designs_three_families.csv` | **RETRACTED** | `results/v2/A1/design_audit.csv`, `docs/progress/v2/A1.md` | Optimizer script used `antoine_A` as proxy for flash point and GCV, omitted WAT and viscosity constraints, and hardcoded `"all_constraints_passed": True`. Jet CSV design freezes at +20.1 °C; Diesel freezes at +26.5 °C; Gasoline density is 1009 kg/m³ with solids melting at +116 °C. |
| **CLM-03** | DeepSets viscosity -20 °C log-RMSE is 0.1633 | Progress report Nov 3 | **REVISED** | Actual bundle evaluation: 0.2720 (`docs/progress/T0.1_repro_note.md`) | Earlier headline claim was sensitive to checkpoint selection and did not reproduce from frozen tag bundles. |
