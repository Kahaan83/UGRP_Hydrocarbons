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
