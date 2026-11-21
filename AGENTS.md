# UGRP Hydrocarbons

Predict fuel blend properties (Jet A/A-1, diesel, RP-1) from composition, then inverse-design blends. Authors: Kahaan Shah, Namyaa Parmar, Arijit Ganguli, Jayendra Bhalodiya, A. S. K. Sinha (Ahmedabad University).
Docs (read only when needed): `PROJECT_HANDOFF.md`, `docs/GUIDE.md`, `docs/final_numbers_v5.md` (source of truth for numbers).

## Model (v5)
- Viscosity: Refutas blending + optional LightGBM pair residual.
- Flash point: Le Chatelier + Raoult + Antoine + LFL overlay.
- Density: ideal volume mixing + optional LightGBM residual.
- Freezing point: Schroeder-van Laar + UNIFAC.
- Uncertainty: one-sided 90% conformal margins from held-out lab systems.
- Public API: `src/fuelblend/api/predict_route.py` (`predict_fuel()`).

## Commands (Windows, PowerShell, Python 3.12 in `.venv`)
- Tests: `.venv\Scripts\python.exe -m pytest -q` (run one file/test, full suite only at the end)
- v5 baseline check: `python scripts/v5/base/reproduce_base.py`
- v4 pipeline: `python scripts/v4/reproduce_all.py` (do not rerun unless asked)

## Data
- Canonical units: cSt, mPa.s, kg/m3, degC, kPa, MJ/kg, mm. Temperatures in K inside equations.
- Bases: `w` mass (viscosity, density inputs), `x` mole (flash, VLE), `v` volume (pour point).
- Antoine: log10(P/kPa) = A - B/(T_K + C). NIST bar to kPa: A_kPa = A_bar + 2.
- Always `GroupKFold` on `system_id`. Dedup on (property, system_id, composition, T, P): median, keep spread and flags.
- Production data: `data/unified/`. `data/quarantine/` is formula-derived, never train on it.

## Traps
1. Viscosity log units: NIST ln(Pa.s) vs Leenhouts ln(mPa.s); offset ln(1000).
2. ThermoML: constraint T is in `Constraint > ConstraintID > ConstraintType`; liquid phase only; skip excess properties.
3. Never evaluate un-refit Antoine constants for C10+ alkanes in the 40-100 degC flash window.
4. Pandas 3: no `.loc[k, [cols]] = [...]` with mixed types; assign columns one by one.
5. Skip `__MACOSX/._*` and `~$*` files in parsers.

## Rules
1. Never invent or alter numbers; every figure must trace to `docs/final_numbers_v5.md` or committed results.
2. State "status: no skill shown" where a model does not beat a constant (real-fuel freezing point, cold viscosity).
3. Design certification is model-level only; no lab measurements exist.
4. Keep `data/unified/`, `scripts/v5/`, `src/fuelblend/` unchanged unless the task requires it.
5. Academic, direct tone. No em dashes in manuscripts or decks. No AI attribution in commits or PRs.
6. Keep personal names, emails, and project details out of external services.

## Usage limits
- Never bulk-read `data/`, `results/` CSV/JSON, `.venv/`, `.docx`, `.pdf`. Use Grep, `head`, `Read` with `limit`, or short Python summaries.
- NO subagents, `/code-review max|ultra`, or web fetching unless asked.
- Cap command output (`Select-Object -First 40`); send long logs to `scratch/`.
- Sonnet  and haiku for edits, Opus only for hard decisions, never `fable`.
- Terse replies, no diff recaps, one recommendation instead of option surveys. New session per task.
