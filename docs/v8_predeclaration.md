# v8 predeclaration: ML-filled component constants and design certification

Date: 2026-11-23. Written before any v8 script or result exists.

## User decisions (recorded 2026-11-23, before any v8 result)

* G1 tolerances: flash 2 C, freezing 2 C, density 2 kg/m3, viscosity 0.02 ln.
* `splits/locked.json`: restore only this file to its committed version in Task 0 (BASE files untouched); locked readings (G4) run.
* Git: commit after each task; push to `origin feat/external-datasets-ingest` once, after Task 7.
* Execution: native, no subagents.

## Gates

| gate | question | pass rule (fixed) |
|---|---|---|
| G1 substitution | Does replacing a measured constant by its out-of-fold ML value keep blend accuracy? Per property: density (`rho15`), viscosity (Andrade from ML log viscosity at 293 K plus the v7 descriptor-ridge slope), flash (LFL from ML pure flash point via `Psat(Tf)=LFL*P`, only where measured Antoine exists), freezing (`Tm`, `dHfus`) | upper end of the paired group-bootstrap 90 % CI of dMAE (ML minus measured) at or below: density 2 kg/m3, viscosity 0.02 ln, flash 2 C, freezing 2 C |
| G2 real-fuel margin | Can a margin be earned for represented blendstocks on real GCxGC fuels (AJFTD + AV-33-22)? | all three: (a) the route beats the leave-one-group-out constant on the same fuels (paired group bootstrap, CI of the gain excludes 0); (b) margin finite; (c) the WS6 R-WS6-3 leave-source-out coverage reading holds in both directions. Coverage alone never passes |
| G3 no regression | With the v8 margin set off, are v5 outputs unchanged? | `reproduce_base.py` passes and re-evaluating the committed v5 design weights with `evaluate_recipe` gives identical constraint statuses and values |
| G4 locked | Do G1 passes hold on locked systems? | same rule as G1, through `src/fuelblend/v7/guard.py`. Viscosity and density lab locked systems were already read by v5 R-MLV-2 / R-MLD-2, so their G4 rows are labelled `already_read` |
| G5 design | How many of the 16 v5 rows (plus P4 x 2 specs if P4 forms) certify under v8? | descriptive count; no target. Each certified row lists which v8 margins it used |

## Reading text

* Pass: "G<k> passed on <set>: <statistic> <= <tolerance>." The constant kind may be used in later tasks.
* Fail: "G<k> failed on <set>: <statistic> > <tolerance>." The kind is not used; affected constraints stay `uncertified`. Failed gates are reported, not tuned.
* Real-fuel freezing and cold viscosity keep "status: no skill shown" unless the Task 3 margin passes G2.

## Unattended execution

A failed gate is a reading, not a stop. Stop only on a broken invariant (`reproduce_base.py` failing, a test that cannot pass without changing a gate, any write into `results/v5/BASE/`); then write `results/v8/STOPPED.md` and do not push.

## Ruling (Task 0)

`git restore splits/locked.json` was a no-op (file already matches the commit, LF). `test_committed_splits_match_rebuild` still fails because the `input_sha256` of the AJFTD long/wide inputs in `results/v5/S4/splits_ajftd_v5.json` differs from the working-tree files (line-ending or byte difference in the inputs). This is pre-existing, unrelated to v8 and not fixed here (no change to data or BASE). The Task 7 full-suite expectation is therefore 741 + new tests passed with that one known failure.
