# Audit of branch v5 (S1-S7, G1, L1, SLE2, P1/D4), 2026-11-11

Read-only audit. No repo script or test was run (they rewrite `results/`, and the S4 locked stage is run-once). One probe imported `fuelblend.mixture.common` and read `data/cache/pure_grid.npz` to look up pure viscosities; the cache file's modification time (2026-11-11 03:20) is unchanged. Every number below was recomputed from the committed CSVs with throwaway scripts outside the repo, or is quoted from a README and marked as such. The only file written is this one.

## 0. Do these before any paper claim is written

| # | Problem | Evidence | Effect |
|---|---|---|---|
| B1 | S3 silently lost the 40 C viscosity of all 38 AV-33-22 fuels. `pred_M` is NaN with `status_M = ok` for 38/38 rows. | `results/v5/S3/s3_predictions.csv`; cause: the "n-alkane C3 and lower" bin (0.000 mass) maps to propane, whose nu at 313 K is 0.193 cSt, so `ln(ln(nu+0.8))` is NaN and `0 * NaN` poisons the VBN sum (`scripts/v4/build_viscosity_leaderboard_v4.py:151-160`; zero-mass reps are not dropped at `s3_composition_vs_physics.py:~301`). | The S3 README sentence "AV-33-22 has -40, -20 C viscosity only" is false. The 40 C row (n = 11, AJFTD only) could be about 46 groups. AV-33-22 has measured viscosity at 9 temperatures (-40 to 100 C); S3 used 2. |
| B2 | S3 freezing headline includes 5 censored ("<") AJFTD values as if measured (fuels 275, 281, 293, 348, 24567). S5 excludes them; S3 and G1 do not (`targets()` in `eval_ajftd_gcxgc_v5.py` has no qualifier filter). One of them (AJFTD:293) has an M error of -103 C. | Recomputed without them: n = 48, M bias +8.2 / MAE 8.2, R bias +1.6 / MAE 4.7, CI of R-M (-4.7, -2.0). | Verdict direction unchanged, but the published row (n = 53, M MAE 10.1, R MAE 5.6) is wrong. G1's freezing "MAE 15.2 vs 7.3" has the same defect (S5 already says so). |
| B3 | S4 heat of combustion is the only T3 "beats constant" verdict, and most of the labels are probably calculated from the model's own features. | Method per row in `ajftd_fuels_long.csv`: D3338 773 rows (D3338 is an estimate from density, aromatics and distillation), codes 400C/400D/400E 1,139 rows (sd 0.118 MJ/kg; meaning of the code not verified), 231A/B 303, D4809/D240 (bomb calorimeter) 365. Dev LOGO ridge MAE on D3338-only fuels is 0.053 against 0.130 for the constant. The locked set is 162 of 215 fuels from one cell (JA 2011), all code 400. | Circular labelling risk. On D4809-only fuels the gain survives in point estimate (dev n = 122: hgb 0.218 vs const 0.378; locked n = 26: 0.249 vs 0.520), but the constant there is a mean over mixed methods (D4809 mean is 0.34 MJ/kg above D3338), so part of that is a method offset. Needs a declared re-run on measured-only rows with a within-method constant. |
| B4 | The two flash correlations in S4 have coefficients "quoted from memory, UNVERIFIED" (README). | `results/v5/S4/README.md` section 1 | "Published correlations lose to a constant" cannot be claimed until the forms are checked against the source. |
| B5 | No v5 claims exist in `docs/claims_register.md` (last entry CLM-41); there is no `scripts/v5/reproduce_all.py` and no `SHA256SUMS_v5.txt`. | file listing | Plan Phase 3 items; nothing on v5 is drift-protected. |
| B6 | "Declared before scoring" is not checkable from git for any item: each README (S1, S2, S4 included) has a single commit that also carries the results. S4 alone has a mechanical guard (README hash in `dev_meta.json` equals the committed file: verified, `dce7250a...`), which proves only that the README did not change between the dev stage and the locked run. S1's README states that one sentence of its declaration was edited after scoring. | `git log` | Say in the methods that readings were self-declared; do not call it pre-registration. |

## 1. Scientific audit

Spot recomputations (all from CSVs): S3 flash row (bias -0.36 / +11.6, MAE 2.58 / 11.67, over-by-5 0.073 / 0.891, CI 7.17-10.79) matched. S3 density, both cold viscosity rows and the 40 C row matched. S5 headline (C 6.41, I 7.99, F0 7.95, F1 3.46, P 4.80; F1-vs-C CI 0.44-4.52) matched. S7 factors (1.569 / 3.508 / 2.434), GP coverage before (0.770 / 0.694) and after (35/36), ridge 30/36 matched. S4 locked heat (const 0.183, hgb 0.094, ridge 0.111) matched. S1 margins (0.4154 = 19th of 20 system maxima; 0.6084 = 16th of 16) matched. S6 eligible-4 rows (density +5.8 / 15.0, viscosity -0.41, flash -12.8) matched. G1 flash (n 37, -0.43, 3.54, 0.108) matched. Not matched: S3 freezing (B2).

### Per item

**S1 (cold viscosity margin).** Numbers follow. Limits the text understates:
* T_263_283K: k = 16 of 16, so the margin (0.608) is the single worst system; the next two maxima are 0.376 and 0.172. One system sets it.
* Both cold margins are essentially one source (Boehm 13-14 of 20 and 14 of 16 systems); components are shared across the nested pre-blends, so "independent systems" is by InChIKey set only.
* The margin is two-sided while the error is one-sided: dev bias -0.11 ln below 263 K, and on real fuels -0.27 ln at -20 C. Under-prediction is the unsafe side of a maximum-viscosity limit. Check against real fuels: 4 of 55 S3 headline fuels have M error below -0.415 (5 of 55 outside +/-0.415), so the 0.415 margin happens to cover about 91 % on real fuels, but only because it is wide, not because the model is centred.
* The Arrhenius rule definition is the author's assumption (disclosed).

**D3 / D4 / P1.** "0 certified" follows. D4's v5pool was chosen after D3 failed (disclosed). Freezing clears the limit by 0.5-0.7 C against a best-case SLE margin of 12.3 K, and the Jet A recipe has 0.3 wt% aromatics, so category B should not be described as "meets every limit" without those two qualifiers. trans-decalin dHfus "verified" rests on two values that probably trace to one measurement (disclosed).

**S2 (flash).** A: holds, with the caveat already written (the "best" mesitylene set is still a 31-42 K extrapolation checked at one point at 437.8 K). I checked whether mesitylene is a one-off: of the 79 pure inputs in the flash run only 22 have a reference Psat ratio at 37.8 C (`results/v5/L1/F_r3/pure_inputs_used.csv`); among those only mesitylene is off (0.52), the rest are 1.00-1.20. The other 57 have no reference in that file. A separate boiling-point test I ran (Psat at Tb vs 101.3 kPa) resolved 59 of the 79 used sets, all of the `anchored_thermo_ambient` kind (it did not cover the refit or WebBook sets, so not mesitylene, tetralin or indane). It flags 12, of which 11 are C13-C18 species far above the 273-373 K anchoring window (expected extrapolation); the one near-window outlier is trans-1,4-dimethylcyclohexane (0.36 at 392 K), worth a look. C: R-C2 "MET" is correct by the rule, but the error is an under-prediction of about 1 C, CoA flash is D56 for 27 fuels and IP 170 for 10 against a model calibrated to neither, and the corrected-mesitylene variant was not propagated to CRC 647. The undeclared r = -0.51 with D86 T10 is labelled as an extra (good).

**S3 (composition vs physics).** Central flash and density rows follow from the CSV. Problems:
* B1 and B2 above.
* Commit message and R-S3-5 say "cold viscosity is physics-limited"; R-S3-1's own declared label is not met at headline level (README says so). The label was effectively redefined in the prose. Also "physics" is not isolated: the lab-mixture bias below 263 K is -0.11 ln, the real-fuel bias is -0.27 / -0.64, so the mixing rule as seen on lab mixtures does not obviously explain it. The cause is not isolated; candidates are the representative choice, pure-component viscosity fits extrapolated below their data, and the mixing rule on 40-component mixtures. One relevant check: a Walther line fitted to each AV-33-22 fuel's own measured viscosity at 15-100 C extrapolates to -20 C with +0.075 ln (sd 0.028) and to -40 C with +0.14 (sd 0.065), so real fuels follow D341 closely and the temperature dependence of the fuel itself is ordinary. Note also that "unbiased at 40 C" rests on 11 AJFTD fuels, and the FACE diesels were -0.41 ln at 40 C.
* The verdict labels "R better than M" were added after scoring (disclosed).
* M and R use different representative tables (SubReps vs RealisticReps) and R is r2, not r3, so "R minus M" is composition plus representative choice (disclosed), and the +11.6 C R flash bias cannot be compared with the +4.2 C quoted for r2 on CRC 647. The light-end explanation is untested.
* Bulk-match grouping tolerances are stated as fixed in advance; not verifiable. Strict-group sensitivity is given, which is adequate.
* Unreported result sitting in the CSV: on the 34 uncensored AV-33-22 fuels the unfitted R route (ideal SvL on a D86-reconstructed composition) has freezing MAE 3.5 C against 5.9 C for a leave-one-out constant (my bootstrap CI of the difference 0.8-4.1); pooled n = 48: 4.7 vs 6.3 (0.2-3.0). That contradicts "freezing models show no skill over a constant" for this source and should be examined, not ignored. It was not a declared reading, so it is a lead only.

**S4 (AJFTD bulk).** Protocol is the strongest in the repo (hash guard, run-once). Limits: B3, B4; locked set is 69 % one cell and the bootstrap runs over 7-8 cells, so "no skill shown" on locked for flash, freezing and viscosity is low power, not a null (hgb locked viscosity MAE 0.323 vs 0.654, CI -0.97 to +0.03). Flash mixes D56, D93/IP170 and 924 "220x" coded rows in one target; the constant's MAE differs fourfold by family (1.5 vs 5.3), so a method-stratified constant is the fair baseline. JFK freezing default values are disclosed. Only 3 flash and 1 aromatics rows are censored; 64 freezing rows are ("<"), and S4's handling of them was not checked.

**S5 (freezing).** Conclusion ("null not shown turned") follows and is stated carefully. Remaining issues: F1's win depends on one fuel (disclosed); freezing methods are mixed (D5972, D2386, D7153, IP 529, CoA); P's alpha hit the grid value 10 in every fold. n-alkanes C8-C10, C17, C19 have no audited constants and are treated as inert, which is a model handicap, not a physics result; the F0x/F1x sensitivity shows it matters (F1x MAE 3.07).

**S6 (diesel).** Gate result follows (9 of 9 by the declared rule). The 10 % unmapped cut-off is not in the declared readings; it is the G1/S2/S3 value and is disclosed, but it reduces R-S6-2 to 4 light fuels that are not representative. With n = 4 the commit-message phrase "registry-based predictions poor" is descriptive only. The -12.8 C flash bias is most likely the PIONA C3-C6 trace bins dominating the Le Chatelier sum (not tested). The transcription check is a self-reported `match` column; I did not re-open the PDFs. The 0.186 margin is hard-coded in `s6_diesel.py:397`.

**S7 (viscosity loose ends).** Arithmetic reproduces. Validity problems: (i) the dev threshold 0.88 cannot fail (42/45 by construction); (ii) the locked threshold 0.80 was set after v4 had published locked coverage 0.69 on the same 12 fuels, and those 12 CRC 647 locked fuels have now been read in v4 V, R2, R3, L1 ("shown") and S7; (iii) 36 locked rows are 12 fuels x 3 temperatures. "Adopted" should read "passes a rule that had little power to fail". The constant-width reading for ridge is an after-the-fact interpretation (disclosed). The locked 20 C pairwise rows are correctly flagged as unreliable.

**G1.** Flash and viscosity follow. Freezing is contaminated by censored values (B2). The 11 S4-locked fuels were read here and again in L1, S3 (apart) and S5 (apart).

**L1.** "Safe to adopt, no accuracy gain" follows. The CRC 647 dev over-prediction count moving 6 to 2 of 45 is 4 fuels; do not present it as a safety result.

**SLE2.** Follows. The 12.3 K margin is k = 11 of 11 (worst single system). The "two compilations" tier may be one measurement (disclosed).

### Units and bases
Seen in code and consistent: flash paths convert mass to mole fractions with RDKit molecular weights before the Le Chatelier sum (`predict_comp` in S3, `mass_to_mole` in S2); Refutas blends VBN on mass fractions of kinematic viscosity (dynamic / density at the same T); density is ideal-volume (`1 / sum(w/rho)`); SvL uses mole fractions and K; registry Antoine is labelled `log10(P/kPa)`, and the S2 comparison reproduces 101.2 kPa at the mesitylene boiling point with the WebBook set read as kPa. Known approximation: Boehm volume fractions are converted with 15 C densities for measurements at 233-288 K (disclosed in S1). Not checked: the AJFTD viscosity temperature parsing, the API-gravity density reference temperature (288.71 K rows are scored as 15 C within a 0.6 K tolerance), the FACE diesel unit conversions, and the UNIFAC path.

### Leakage and double counting
* No untouched T2 fuel remains: AV-33-22 was read in S2, S3, S5; AJFTD GCxGC in G1, L1, S3, S5; FuelLib in S2, S3, S5.
* AV-33-22 vs AJFTD GCxGC overlap was checked by composition distance (min 0.072); AV-33-22 vs the AJFTD bulk table (2,672 fuels, T3) was not. Needed before any T2 + T3 pooling.
* AV-33-22 -20 C viscosity is a median of UDRI-measured and supplier-CoA values (75 rows for 38 fuels).

## 2. Code audit

| Sev | Finding | Where |
|---|---|---|
| High | Zero-mass representative gives NaN Refutas with status "ok" (B1). `refutas_nu` has no guard for nu + 0.8 <= 1; the caller does not drop zero-weight components. | `scripts/v4/build_viscosity_leaderboard_v4.py:151`, `scripts/v5/s3_composition_vs_physics.py` (AV-33-22 block). S5 keeps zero-mass rows the same way but its freezing models are not affected. |
| High | Censored values not filtered in `targets()`. | `scripts/v5/eval_ajftd_gcxgc_v5.py:86-96` (used by G1 and S3) |
| Med | Two mapping implementations. `map_fuel()` calls `SubReps.get`; S2 re-implements the rule (`rep_for`) on a table regex-parsed out of a markdown report (`results/v4/R3/reconstruction_report.md`), and S3/S5 consume S2's CSV output. A report edit changes the science. S2 also marks a bin unmapped when the rep lacks Antoine/LFL, then S3 reuses that flash-specific map for viscosity and density. | `eval_ajftd_gcxgc_v5.py:35-66`, `s2_t2_flash.py:64-115` |
| Med | Scripts import each other by file path (`importlib` of `scripts/v4/...` and `scripts/v5/eval_ajftd_gcxgc_v5.py`); shared functions (`map_fuel`, `refutas_nu`, `ideal_sle_C`, bootstrap, conformal factor) belong in `src/fuelblend`. Four separate paired-bootstrap implementations (S3, S4, S5, S7). | S3, S5, S6, S7 scripts |
| Med | Loaders and two scripts read from gitignored `extras/` (`loaders_ajftd.py:36`, `loaders_ajftd_gcxgc.py:24`, `build_lfl_overlay_v5.py:20`, `eval_sle_extended_v5.py:45`), so a clean checkout cannot rebuild the v5 tables. | as listed |
| Med | Coarse class mapping: every cycloaromatic C11-C15 goes to tetralin (C10) and every dicycloalkane C8-C18 to trans-decalin (gap up to 8); that is about 4.7 % and 10 % of mean AV-33-22 mass with a carbon number that is too low, which should lower predicted viscosity and raise volatility; the size of the effect was not computed. | `results/v5/S2/av3322_representative_map.csv` |
| Low | Bare `except Exception` in `audit_antoine_range_v5.py:113`, `build_lfl_overlay_v5.py:36`, `eval_sle_extended_v5.py:117,120`. Hard-coded margin in `s6_diesel.py:397`. No absolute user paths found. | as listed |
| Low | FuelLib HEFA-blend rows carry `T_C = NaN` (`loaders_fuellib.py:192-204`; the loader comment says the source file does not state the temperature), so 48 property rows are unusable until the temperature is confirmed from the paper. | loader |

Tests: mostly behavioural and adequate (S5 censoring, LOGO isolation, S7 order statistic). Weak ones: `test_s3_composition.py:54-59,89` and `test_s6_diesel.py:71-75` assert on source text; `test_s6_diesel.py` is skipped wholesale if a CSV is absent and `test_v5_data_audit.py` uses `importorskip`, so a broken import passes silently; the transcription test checks a self-written `match` column; S6 unit tests check value ranges only. Missing: a test that a prediction with status "ok" is finite; a test that S3/G1 targets exclude qualified values.

## 3. Better use of the data

Inventory and verdict:

| Data | Used for | To potential? |
|---|---|---|
| AJFTD bulk, 2,672 fuels | S4 only (7 features, 4 targets) | Partly. Method codes not exploited; 578 assumed-T viscosity fuels set aside. |
| AJFTD GCxGC, 113 fuels | G1, S3, S5 physics baselines | No. 75 fuels by my count from the `has_*` flags (the handoff says 74) have composition and no property row; nothing trained. |
| AV-33-22, 38 fuels | flash, density 15 C, viscosity at 2 of 9 T, freezing | No. Viscosity and density at 9 T, vapour pressure (266 rows), heat capacity, surface tension, smoke, heat of combustion, simdis are not read by any `scripts/v5` file (grep). |
| FuelLib | 4 POSF fuels | No. 18 Jet-A/HEFA blends (plus neat Jet-A and 3 neat HEFAs) with density and viscosity and computable compositions unused (temperature not stated in the source). |
| Boehm 2022 | S1 margin | Partly. Never used to fit or test a cold pure-component correction. |
| FACE diesel (9) + NIST surrogates (4) | S6 | Limited by representatives (up to 38 % mass dropped). |
| LFL overlay, SLE2 constants | side tables | Not adopted (user decision pending). |
| NISTIR 6659, NISTIR 6646, wax appearance T | none | Unused. |

Proposals, ranked by value for the paper over effort. All are code-only unless marked.

1. **Fix B1 and B2 and re-issue S3 with all 9 AV-33-22 temperatures** (hours). Gives the paper a viscosity-error-versus-temperature curve on 38 measured-composition fuels (-40 to 100 C) and density likewise. This directly separates "representative/mixing error at warm T" from "cold extrapolation error". Add as a dated note with new readings; the fuels are already read, so label it exploratory.
2. **Cold-viscosity correction from AV-33-22** (1 day). Evidence that it will work: on AV-33-22 the -20 C error has sd 0.075 around -0.264, so a leave-one-out constant offset gives MAE 0.061 ln (from 0.264); at -40 C 0.110 (from 0.637). The offset transfers across sources: AJFTD bias is -0.280 at -20 C and -0.664 at -40 C. Better form: Refutas at 40 C (unbiased, +0.014) plus a Walther slope, with the slope correction fitted on AV-33-22 and tested on AJFTD GCxGC (leave-source-out, reading declared first). Also diagnose the cause on Boehm pure compounds at -20/-40 C: if the pure fits are already low there, fix the provider instead.
3. **Re-run S4 heat on measured (D4809/D240) rows only, with within-method constants; stratify flash by method family** (half a day). Decides whether T3 has any positive result (B3).
4. **Pooled T2 trained model with leave-source-out validation** (2-3 days). Labelled fuels are far fewer than compositions: about 79 for flash (AJFTD 37, AV-33-22 38, FuelLib 4), 54-56 deduplicated groups in the S3 headline, fewer for density and 40 C viscosity. That supports a low-capacity residual (ridge on class and carbon-number moments, as the v4 structured ridge does) for flash, density and viscosity, not a deep model. Train on AV-33-22, test on AJFTD and the reverse. This is the "physics baseline + learned residual" the project set out to build and currently only exists on reconstructed CRC 647 compositions.
5. **Build a fresh holdout from the 75 composition-only AJFTD fuels and POSF 4658** (1-2 days, may need downloads). These are the only T2 fuels not yet read. Labels are POSF numbers (3895-4461 series and 10326-13718); link to NISTIR 6659 (121 density rows for 4658, on disk), AFRL/Edwards reports (download, user approval), CRC AV-27-18 rows already ingested. Declare readings before joining properties.
6. **Supervised reconstruction** (2 days). R flash is +11.6 C because the D86 route has no C5-C8 tail. With the 71 fuels that have both bulk descriptors and a measured composition (33 AJFTD + 38 AV-33-22), fit the map from (D86 points, aromatics, density) to the light-end fraction and n-paraffin share, leave-source-out, then apply to CRC 647 and re-score the 45 dev fuels. Also use AV-33-22 simdis (D2887, 21 points per fuel, unused) as a better carbon-number source than D86. This would upgrade T3 and CRC 647 evidence at once.
7. **Freezing: follow the R-route lead and calibrate SvL** (1 day). Unfitted R already beats the constant on AV-33-22 (3.5 vs 5.9 C). Test whether D86 T90/FBP plus n-paraffin share explains it; fit F1 with the complete Coutinho n-alkane table (F1x MAE 3.07) and add the 57 wax-appearance points and SLE2 n-alkane systems as T1 evidence for the dHfus/solid-solution term. Needs declared readings; AJFTD n is 13-15, so the claim will stay source-specific.
8. **Better representatives** (1-2 days). Add 1-methylnaphthalene / 2,6-dimethylnaphthalene for diaromatics, a perhydrophenanthrene or perhydrofluorene for tricycloparaffins, alkyl-tetralins and alkyl-decalins at C12-C14, and replace mesitylene as the single C9 aromatic. Requires verified Antoine/LFL/viscosity for each (partly in ThermoML; fusion data for dimethylnaphthalenes would need the Acree-Chickos download). Recovers the 5 FACE diesels now excluded and removes the tetralin/decalin carbon-number bias. One shared `map_bins()` in `src/fuelblend` at the same time.
9. **Jet-A/HEFA blends from FuelLib** (half a day plus confirming the measurement temperatures from the Vozka paper; library access or user confirmation). 18 blends of two fuels with measured GCxGC: a real-fuel blending-rule test, which is the paper's stated subject. First check overlap with the Vozka 2019 tables already in `extras/.../raw/vozka2019/` (not checked here).
10. **Adopt or reject the LFL overlay** (user decision). If adopted, regenerate flash margins once, together with the mesitylene fix, so flash numbers change a single time.
11. **AV-33-22 vapour pressure (266 rows)** as a direct test of the Raoult/Antoine half of the flash model on measured compositions (half a day). Separates Psat error from LFL error.
12. **NISTIR 6646 RP-1** link (composition covers 80.6 wt %, below the 0.95 gate): supplementary only. **Leenhouts published splits**: cheap comparability check for T1.

Dev/locked structure: the CRC 647 locked 12 are exhausted (five reads) and should be retired as "locked"; report them as a second dev fold. S4's locked set is one cell heavy; for any new T3 run use leave-source-group-out as the primary and keep the existing locked result as already-read. For T2 use leave-source-out (AV-33-22 / AJFTD / FuelLib) plus the fresh holdout in item 5.

## 4. Not verified

* Whether any README was written before its scores (except S4's hash guard).
* The meaning of AJFTD method codes 400C/D/E, 231A/B and 220x.
* The S4 correlation coefficients, the S6 transcription against the PDFs, the Boehm volume-to-mass conversion, and the S5/SLE2 fusion constants against their sources.
* Whether the pure-compound cold viscosities in the provider are measured or extrapolated per compound (I only confirmed that real fuels follow Walther and that the provider returns finite values).
* S4's treatment of the 64 censored freezing rows.
* No script or test was executed, so "383 passed" is taken from the handoff note.
