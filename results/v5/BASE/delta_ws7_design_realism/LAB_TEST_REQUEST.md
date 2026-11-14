# Lab-test request (WS7 final check): NEEDS USER DECISION

Status: draft only. Nobody has been contacted. Whether to test, which lab, and whether these blends can be made at all are the user's decisions.
No design here is certified. Every number below is a model prediction (seed 42) from `results/v5/BASE/delta_ws7_design_realism/r_designs.csv`.

## Why test

No inverse design in v5 is certified. The two spec flags (aromatics vol % minimum, smoke point) are unverified, and the real-blendstock designs also
lack validated margins. One measured designed blend would test the model chain directly (composition -> properties), which no margin can do.

## Candidate 1 (recommended): pool C Jet A design (verified pure compounds only)

* Recipe (mass fraction): butylcyclohexane 0.5568, n-decane 0.1597, n-undecane 0.1357, 1,2,4-trimethylbenzene 0.0934, n-dodecane 0.0283,
  isooctane 0.0261 (trans-decalin 0). This is the 6dc66f8 / R-W7-5 lambda-0 recipe. The lambda-0.25 variant (BCH 0.579, TMB 0.093) is equally usable.
* Predicted: density 784.0 kg/m3, flash 40.9 C (D56 limit 38; margin 2.21 C set by one lab system), freezing -57.7 C, viscosity at -20 C
  2.90 cSt, NHOC 43.80 MJ/kg, aromatics 8.32 vol %, smoke 24.9 mm (needs the D1655 18 mm + naphthalenes alternative).
  D86 (this delta, unverified limits): T90-T10 is far below the D7566 40 C minimum (narrow-boiling recipe), so it is **not a D7566-conforming
  synthetic blend**. Test it to check the property models, not as a fuel candidate.
* Why first: every constant is verified and the counterfactual says only the spec flags block it. Pure compounds can be bought.
  The flash, freezing and NHOC slacks are thin (about 0.019 normalised), so a measurement tells us directly whether the margins hold.

## Candidate 2 (optional, harder to source): P2 SAF design for Jet A-1

* Recipe (mass fraction): Jet A base POSF-10325 0.490, ATJ-SPK (POSF11498 type) 0.403, HEFA-SPK (camelina type) 0.098,
  1,2,4-trimethylbenzene 0.010 (synthetic share 0.50, at the D7566 cap).
* Predicted: density 775.2 kg/m3 (edge), flash 50.0 C, freezing -46.96 C (**fails the -47 C limit by 0.04 C at the point estimate**),
  viscosity at -20 C 3.56 cSt, aromatics 7.99 vol % (**at the 8 vol % minimum, fails by 0.01**), smoke 25.0 mm, D86 T10 178, FBP 292,
  T50-T10 35.5, T90-T10 85.8 C.
* Caveat: it needs the same three real blendstocks (or close equivalents, with their own GCxGC). The prediction uses GCxGC representatives with
  non-verified constants. Use it only to test the blendstock route, and only if the user can obtain the streams.

## Tests to ask for (per candidate)

| property | method | why |
|---|---|---|
| flash point | ASTM D56 (Tag closed cup) | thinnest-margin safety property; over-prediction is the unsafe side |
| kinematic viscosity at -20 C | ASTM D445 | cold viscosity; under-prediction is the unsafe side |
| density at 15 C | ASTM D4052 | composition basis, spec 775-840 |
| freezing point | ASTM D5972 | margin measured on binaries only |
| smoke point | ASTM D1322 | smoke route via the real-fuel margin and the alternative limit |
| (optional) D86 distillation | ASTM D86 | tests the new D86 baseline and margins (this delta) |

## Questions for the user (needs user decision)

1. Can the university lab run D56, D445 at -20 C, D4052, D5972, D1322 (and D86)? What sample volume is needed?
2. Test Candidate 1 only, or both?
3. Who prepares the blend gravimetrically, and what purity is acceptable (record the purities of the bought compounds)?
