# Scope Decisions Log (`docs/scope_decisions.md`)

This document records formal scoping decisions governed by empirical data availability and scientific requirements defined in the task specification.

---

## 1. Decision SD-01: Scoping Out Pour Point Modeling for Hydrocarbon Mixtures (Task T2.12)

- **Date:** 2026-11-04
- **Task ID:** T2.9 – T2.12 (Low-Temperature Behavior)
- **Governing Rule:**
  > "T2.12: Search for pour point data with known composition. If fewer than 30 mixtures with composition are found, record the decision to scope pour point out and document it in `docs/scope_decisions.md`."

### Data Audit Findings:
1. **Speciated Chemical Mixtures with Measured Pour Point:**
   - NIST ThermoData Engine (TDE): **0 mixtures** (NIST focuses exclusively on thermodynamic equilibrium properties; pour point ASTM D97 is a kinetic/rheological flow-cessation test).
   - NIST TRC ThermoML Archive: **0 mixtures**.
   - Leenhouts et al. Chemical Mixture Benchmark: **0 mixtures**.
   - CheMixHub Mixture Repository: **0 mixtures**.
   - **Total speciated chemical mixtures found with measured pour point:** **0** (threshold requirement: $\ge 30$).
2. **Real Fuels with Bulk / Semi-Speciated Data:**
   - CRC Report No. 647: 57 fuels with ASTM D5949 pour point, but only bulk D86 distillation curves, not speciated molecular formulations.
   - NASA NTRS 19860018789: 8 fuels with n-alkane mass fractions and static equilibrium pour points.
   - Kirgina et al. (2019/2021): 10 diesel fuels with pour points and bulk boiling cuts.

### Formal Decision:
**Pour point is formally scoped out of the mixture machine learning target list.**

### Scientific Justification:
- Pour point (ASTM D97 / D5949) is an empirical, non-equilibrium standard test measuring the temperature at which a gelled wax crystal network halts gravity-driven fluid movement. Unlike thermodynamic phase transitions (e.g., freezing point, cloud point, or wax appearance temperature), pour point depends heavily on cooling rate, shear history, and trace wax crystal modifiers/pour point depressants (PPDs) rather than pure equilibrium solution thermodynamics.
- For low-temperature aviation and diesel fuel modeling, **freezing point (ASTM D5972 / D2386)**, **cloud point (ASTM D2500)**, and **wax appearance temperature (SLE Schroeder-van Laar)** are retained as the authoritative, thermodynamically grounded cold-flow targets.

---

## 2. Decision SD-02: Diesel real-fuel evaluation waits on Task 8 documents

- **Date:** 2026-11-10
- **Task ID:** v4 Task 7 (diesel corpus decision and splits)
- **Governing rule:** "Include only fuels with raw-file composition. If fewer than 5 qualify, record that the diesel real-fuel evaluation waits on Task 8 documents."

### Data audit findings:
- Searched `data/raw/pour_point`, `data/raw/gc_composition`, `extras/fuelblend_data/data/raw/pour_point` and `.../wax` for diesel fuels with a raw-file class/component composition (`diesel|LFP|ULSD|EN590` in `fuel_code`/`fuel_name`). **0 qualify** (threshold: 5).
- The only diesel fuels in the repo (NASA LFP-3, LFP-5; the Kirgina diesels; FACE diesels; surrogate diesels) have their compositions and properties typed into v3 scripts or in a non-raw table. Task 3 dropped them (`fuel_dropped:no_clean_composition`) and the CRC 647 / AV-27-18 raw property files hold jet fuels only. The wax folder holds only SLE/Coutinho model outputs, not fuel compositions.

### Formal Decision:
**The diesel real-fuel corpus is empty in v4.** `splits/v4/real_fuels_v4_splits.json` contains no diesel fuel, and no diesel real-fuel evaluation is reported. It resumes once Task 8 documents (raw diesel composition tables, e.g. Kirgina et al. and the NASA LFP report) are obtained and transcribed into a raw file.
