# v6 failures: literature-based causes, fixes, and ML options per property

Status: 2026-11-21. Companion to `docs/v6_summary.md` and `docs/literature_review_v6.md` (which covers v6 motivation; this file does not repeat it).
Sources: Consensus searches run 2026-11-21 (6 queries, about 110 records screened). Only papers returned by those searches are cited;
claims about our own data are hypotheses until tested, and each fix is a proposal with a declared gate, not a result. No repo number is
introduced here; all figures are from `docs/v6_summary.md`. Citations lack volume and page numbers; add before manuscript use.

## 1. Why v6 failed, by failure

### F1. Freezing point (Track A): dev gain, no holdout skill; UNIFAC-SvL biased +5.9 C
Observed: LOGO dev MAE 3.31 vs constant 6.64, holdout 9.30 vs 6.86; gain driven by AV-33-22 fuels (65 % of dev).

Literature reading:
- Hall et al. (2021, Energy & Fuels, M-QSPR) found that for freezing point and viscosity "the presence of fuels in the training data was
  found to be essential". Pure-compound or tile-constant models cannot carry them. We have only 64 fuel groups, so any freezing model is
  data-limited by design, not by method.
- Freezing point depends on the heaviest, highest-melting species and on solid-solution behaviour. Watanasiri et al. (2026, Fuel) improved
  the modified weighted average by weighting heavier compounds for freezing point and lighter for flash point. Boehm et al. (2022,
  Energy & Fuels) and Coburn et al. (2022) show mixture freeze points are strongly non-linear in the mole fraction of the highest-melting
  component and that a Gibbs-energy model is conservative only at some compositions. A mean-type model (our Tm average, bias +24 C) is
  structurally wrong; an SvL form with a wrong tail is wrong in sign (+5.9 C bias).
- Reported good results (Shi et al. 2017: MAE 0.82 C on 17 fuels; Shao et al. 2024: 104 jet fuels) use small, similar, mostly Jet A/RP-3
  fuel sets and random splits; they are not comparable to grouped holdouts across fuel families (Hall et al. and Dong et al. 2026 make
  the same family-dependence point). Our dev set mixed one dominant family (AV-33-22) with others, which is the classic setup where a
  random or LOGO score overstates generalisation.
- Compositional resolution matters. Mehta et al. (2023) convert GCxGC to UNIFAC group composition; Oh et al. (2022, AIAA) use an
  autoencoder novelty score that "preemptively" flags fuels with high prediction error.

Probable causes (to test): (a) too few independent fuels and family confounding; (b) tile constants average isomers whose Tm spans tens
of degrees (our own tile-fill test moved only 31 tiles and only fuel 64); (c) wrong functional form (mean vs. solid-liquid equilibrium).

Fixes, in order of cost:
1. Report freezing as interval-only: the live 13.84 C margin is already the honest output. Optional: a one-sided "over-prediction risk"
   flag for fuels whose heavy-n-alkane tail is outside the training range (novelty flag, Oh 2022 idea).
2. Replace mean-type baselines with a Hall-style M-QSPR on the tile-isomer probability space, trained jointly on pure compounds and the
   64 fuel groups (the method trains without mixing rules; its MC-dropout uncertainty gives the novelty signal). This is the one method in
   the literature built for exactly our data layout.
3. Use a heavy-tail feature set (n-alkane C14+ mass, naphthalene and dicycloparaffin content) with monotone constraints inside a Gaussian
   process, not a boosted residual; GPs gave the most stable results in small-data comparisons (section 4).
4. Evaluate by fuel family (leave-one-family-out), not only LOGO, and pre-declare which families count as dev.

### F2. Track C mixture ML: one of five heads passed dev and locked, then did not replicate
Observed: Leenhouts MON MolPool 5.28 vs descriptor GBM 4.97 on 51 fresh rows; flash oxygenate arm gate not established.

Literature reading:
- Mixture QSPR has a known validation problem: Yeghyan et al. (2026, Mol. Inform.) list "mixture interactions, data scarcity,
  representation difficulties, and validation issues" and say very few studies address real fuel mixtures. Positive reports on binary
  mixtures (Qu et al. 2026: R2 0.973 viscosity, 0.996 density; Song et al. 2026: 332 compositions from 33 binaries and 20 components)
  use random or composition-level splits where the same component pair appears in train and test. Our system-disjoint splits are
  stricter, so lower skill is expected, not a bug.
- Leenhouts et al. (2024, Fuel) report MolPool gains mainly when oxygenated species cause non-linear mixing; for near-ideal hydrocarbon
  mixtures the baseline is already hard to beat. Our hydrocarbon-only focus removes the regime where MolPool helps.
- Saldana et al. (2013) found that ML on pure properties combined with Le Chatelier beat direct mixture ML on the available databank, and
  Gaudin et al. (2015) found mixing rules plus pure-compound QSPR reached 2.9 to 4.4 C on binaries. That is the structure of our v5 route;
  v6 head tuning cannot beat it without more mixture data.
- Passing one of five heads without multiplicity control is expected by chance (our own summary says so).

Fixes:
1. Learn the non-ideality, not the property. Replace the residual head with a thermodynamically constrained excess-Gibbs-energy model
   (Section 3) so predictions stay physical on unseen systems and the learning target is smaller and smoother.
2. Pre-train the mixture head on large synthetic or quantum-derived mixture data (TeNNet-SAC pre-trains on one million synthetic
   points then fine-tunes; DeepThermoMix trains on COSMO-RS data). We have 23 to 241 measured rows per head; this is the missing scale.
3. Control multiplicity: one pre-registered primary head, Holm or Bonferroni across heads, minimum replication set sized in advance.

### F3. Track B: wins vs constant but absolute errors too large (flash MAE 19.1 C, density 21.9 kg/m3, Tm no better than v5)
- Constant is a weak bar (our summary). The right comparators in the literature are group contribution plus descriptor models. Ogbomo et
  al. (2026, Digital Discovery) found XGBoost on molecular descriptors consistently among the most accurate for density and viscosity of
  835 C10-C50 molecules, and that descriptor-space distance flags risky regions only at the population level. Li et al. (2021, Fuel)
  report R2 0.98 across 15 properties with a UOB functional-group QSPR and 1,797 pure compounds. A 1,208-compound D-MPNN with 5-fold
  scaffold CV should be compared with such a group-contribution GBM, not with a median.
- Schweidtmann et al. (2020) and Larsson et al. (2023) needed multi-task, transfer and ensembles to work at a few hundred labels; our 10
  targets have very uneven coverage (oxygenates octane only), so some tasks see few labels.
- Physics-structured outputs beat pure data-driven GNNs: ExPUFFIN (Rebello et al. 2026) lowers viscosity RMSE by 37 % against a
  data-driven baseline by predicting Andrade-type parameters from the graph. Our log viscosity at one T (0.225) discards the T dependence
  the literature says to model.

Fixes: add the GC-descriptor GBM as comparator; predict equation parameters (Antoine, Andrade, Rackett/COSTALD) per molecule instead of
point values; add a GP or ensemble-variance head and use it as an applicability filter before any downstream use.

### F4. Flash oxygenate arm: not established, Le Chatelier not evaluable
No Antoine/LFL coverage for oxygenates in the repo makes the physics comparator drop out. Paricaud et al. (2019, Fuel) predict
flash points of gasoline, diesel, biodiesel and jet fuels with Catoire-Naudet plus COSMO-SAC activity coefficients, which needs no
experimental Antoine or LFL per compound; Liaw et al. (2010) do the same with UNIFAC-Dortmund. Fix: build a COSMO-SAC or ML-activity
comparator for oxygenates (Section 3) so the gate compares like with like.

## 2. Cross-cutting causes

1. **Evaluation geometry.** Reports of high skill in the literature mostly use random or within-system splits. Our grouped designs are
   stricter and the small groups give wide CIs. This is a finding, not a defect; it should be stated in the manuscript.
2. **Data scale.** Real-fuel and mixture sets are 10^1 to 10^2 rows. Methods that work in the literature either use 10^5 rows
   (HANNA: 317,421 activity coefficients) or pre-train on synthetic/quantum data (TeNNet-SAC, DeepThermoMix).
3. **Residual-on-physics is not free.** In a small-data comparison (Grewe et al. 2026), residual learning on a compact physics baseline
   helped Gaussian processes but degraded tree-based models; model rankings from a 15-point holdout flipped under 10-fold CV. Our
   boosted residuals and small locked sets match both warnings.
4. **Shift between dev and lock.** Conformal intervals are valid only under exchangeability. Tibshirani et al. (2019) and group-weighted
   or KMM-based variants (Bhattacharyya 2026; Laghuvarapu 2026, preprint) restore coverage under covariate shift when overlap exists.
   Our freezing margin covered 0.67 of holdout fuels against 0.90 nominal, a shift symptom.
5. **Label quality.** Boehm et al. (2022) fit an error model for the Arrhenius viscosity rule from 675 measurements and report error
   growing with mixture complexity (13.2 % at -40 C for pure-component mixtures); our dedup flags help, but label noise sets a floor
   that more modelling cannot beat.

## 3. Ways to put ML into the pipeline, per property
Principle from the literature (Jirasek and Hasse 2023, Annu. Rev. Chem. Biomol. Eng.): ML should enter at the level of a thermodynamic
quantity or an equation parameter, so physics supplies structure and extrapolation, and ML supplies the molecular and interaction
terms. Pure residual-on-output is the weakest of the hybrid options.

| Property | Current route | ML insertion (literature) | What it replaces or adds | Gate to pre-declare |
|---|---|---|---|---|
| Viscosity | Refutas + LightGBM pair residual | ExPUFFIN-style Andrade parameters from graph for pure; mixture: MolPool or a pretrained 3D encoder (Qu 2026) only as a residual on the Arrhenius rule; use Boehm (2022) error model as the baseline interval | pure-component log-viscosity vs T | grouped by system, hydrocarbon non-inferiority, cold (-20 C) separately; keep "no skill shown" if no beat of constant |
| Flash point | Le Chatelier + Raoult + Antoine + LFL | (a) ML activity coefficients (HANNA, GE-GNN, SPT-NRTL) in place of UNIFAC gamma inside the same Le Chatelier or Liaw model; (b) ML for missing Antoine/LFL (Saldana 2013); (c) Bayesian NN from GCxGC for real fuels (Oh 2024, already in our review) | UNIFAC gamma, missing pure inputs | one-sided over-prediction error reported (Trap 6), not MAE alone |
| Density | Ideal volume + LightGBM residual | Excess volume from a GE-style thermodynamic head; pure density via COSTALD or Rackett with ML parameters; GP for real fuels with uncertainty (Hall 2021 GPR density, which also shows synthetic fuels in training are required for hold-out fuels) | excess volume term | beat ideal mixing on binaries from different families |
| Freezing point | Schroeder-van Laar + UNIFAC | ML activity coefficients in SvL (same HANNA family, solid-liquid data permitting); Hall M-QSPR over isomer-probability tiles; novelty flag from autoencoder or MC dropout | gamma in SvL; tile stand-ins | leave-one-family-out, coverage of margin 0.90 +/- tolerance |
| Octane, cetane, YSI (pure) | not a v5 route | Track B D-MPNN already beats constant (YSI holdout passes); add GC-GBM comparator; Schweidtmann 2020 transfer from larger sets | pure-component lookup | against GC-GBM, not constant |
| Octane (blends) | linear molar baseline | latent-space mixing (Kuzhagaliyeva 2022, Xin 2026) | only if oxygenate or polar blends are in scope | system-disjoint replication set fixed in advance |
| All | per-property margins | Conformal with covariate-shift weights (Tibshirani 2019; group-weighted) and family-conditional coverage (Gibbs 2025) | single global margin | coverage by fuel family, not only marginal |

Thermodynamically consistent activity-coefficient models are the common component across flash, freezing and density excess terms:
HANNA (Specht 2024, Chem. Sci.; extended to multicomponent, Hoffmann 2026, Nat. Commun.), GE-GNN (Rittig 2024), SPT-NRTL (Winter
2023), Gibbs-Duhem-informed NNs (Rittig 2023) and TeNNet-SAC (Yang 2025). All are trained on large Dortmund-type or synthetic data,
so they inherit coverage we cannot generate ourselves. Caveat: Zhang et al. (2026, PEGAT) note that unphysical behaviour can still occur
despite hard constraints, so each must be tested on hydrocarbon-hydrocarbon pairs relevant to jets and diesels before adoption.

## 4. Recommended order (one recommendation, not a menu)
1. **Test an ML gamma (HANNA) against UNIFAC on hydrocarbon-hydrocarbon systems first, inside the v5 flash route only** (the full-text read,
   section 4b, shows HANNA's advantage is demonstrated on random, mostly non-hydrocarbon systems, binary VLE and infinite-dilution gamma
   only, and not for solid-liquid equilibrium). Adopt it for flash only if it passes the v5 gates on the Leenhouts hydrocarbon binaries;
   do not use it in the freezing route until a solid-liquid-capable model is verified. Run pre-trained, no refit on our data.
2. Add the group-contribution descriptor GBM as a Track B comparator (cheap, closes the weak-constant critique).
3. Replace freezing residual with M-QSPR plus GP and a novelty flag; evaluate by family.
4. Pre-train a mixture head on COSMO-SAC or ML-generated mixture data before touching measured mixtures again.
5. Re-calibrate margins with family-weighted conformal; report coverage by family.

Items not recommended: further tuning of MolPool heads on the same 23 to 241 rows; any new head scored on a set already seen.

## 4b. Full-text verification (added 2026-11-21, after the abstract-only pass)
Three sources were read in full (downloaded PDF or EPUB, converted to text); everything else in this file stays abstract-level.

**HANNA, Specht et al. 2024 (arXiv 2407.18011 text and SI, the Chem. Sci. version).**
- Data: 317,421 ln-gamma points in 35,012 binary systems from the Dortmund Data Bank; VLE-derived gamma and gamma at infinite dilution only.
  No solid-liquid, liquid-liquid or excess-enthalpy data in this version (the paper lists them as future work), so it does not by itself
  give a freezing-point route. The Nat. Commun. 2026 extension claims LLE, excess enthalpy and multicomponent; I could not retrieve it.
- Split: random, system-wise 80/10/10 (test 3,502 systems). A test system can reuse components seen in training. UNIFAC could be
  evaluated on only 1,658 of the 3,502; on that shared set HANNA's mean system MAE (of ln gamma) is about one third of UNIFAC's and the
  median falls from 0.09 to 0.05. The authors state the comparison is biased in UNIFAC's favour (its training data are undisclosed
  and probably overlap the test set).
- Unknown components: only 56 test systems have one unseen component (40 on the UNIFAC horizon); both unseen, 1 system. Scores there are
  "significantly worse" than the main result; 70 % of those systems contain water and 80 % are infinite-dilution points, so the subset
  is not representative of hydrocarbon pairs.
- Not reported: any hydrocarbon-hydrocarbon breakdown, temperature-range limits, and calibrated uncertainty. Multicomponent mixtures are
  handled only by suggested binary-subsystem plus pair-interaction models.
- Consequence for the plan: the gain over UNIFAC is shown on random, mostly non-hydrocarbon systems. Whether it holds for the C8-C16
  alkane, cycloalkane and aromatic pairs in jet and diesel is untested, so step 1 of section 4 must start with that check and cannot be
  assumed to help. It also should not be used inside Schroeder-van Laar until a solid-liquid-capable model is verified.

**Hall dissertation (DLR/Stuttgart 2024), chapters 3 and 4, covering the M-QSPR and direct-correlation (DC) models.**
- Freezing-point data (Table 3): 402 pure compounds (1,053 points), 74 conventional fuels (80 points), 50 synthetic fuels (58 points).
  Fuel-level evaluation comes from 4-fold cross-validation on fuels, not family-grouped; three additional hold-out fuels (Jet A-1,
  FT-SPK, ATJ-SPK) were used for the screening test.
- Freezing MAE (Table 12): DC 1.51 C and M-QSPR 2.24 C on conventional fuels; 7.14 C and 7.3 C on synthetic fuels. Interval coverage
  (PICP, target 95 %): 89.7 % and 79 % conventional, 48.65 % and 50 % synthetic. Interval width is several times the D5972
  reproducibility (0.8 C). Large errors below -60 C, so the author limits the domain to above -60 C.
- Mechanism stated by the author: iso-alkane C10 freezing points span -110.15 to -12 C across isomers, so isomer treatment dominates
  freezing error; a sampled-isomer QSPR model failed on synthetic fuels, and the DC model needed similar fuels in training. M-QSPR
  depended least on synthetic training fuels and helped most where few fuel labels exist (cetane); only M-QSPR classified the ATJ-SPK
  hold-out correctly against the spec limits. The author calls freezing point the exception to otherwise accurate predictions.
- Consequence: our dev MAE (3.31 C, LOGO) and holdout MAE (9.30 C) match the pattern in Hall's results (about 1.5 to 2.2 C
  on conventional fuels in random CV, about 7 C on fuels unlike the training set). That supports the reading that our failure is a
  similarity-to-training problem, not an implementation bug, and that PICP below nominal is expected for dissimilar fuels. It is
  consistent with, but does not prove, the family-driven dev gain. Hall's fuel counts for freezing (about 124) are about twice our 64
  groups. Note: the Energy & Fuels paper is dated 2021 online and 2022 in print (lit review uses 2022).

**Boehm et al. 2022, Front. Energy Res. (EPUB).**
- 675 viscosity and density measurements of binary to hundreds-of-component hydrocarbon mixtures. Authors state error depends more on
  mixture complexity and temperature than on size or saturation differences, and log viscosity is linear in mole fraction within
  noise. For mixtures whose components are all identified pure molecules the Arrhenius-rule error is 13.2 % of nominal viscosity at -40 C
  and 9.6 % at -20 C, multiplied by the root mean square of component mole fractions (so small for many-component fuels). Our spec
  condition is -20 C, so 9.6 % is the relevant figure. This is an error model for the blending rule, not an ML model, and is the
  correct baseline interval for any viscosity ML claim.

**Not retrieved (still abstract-only):** Hoffmann 2026 (Nature redirect), Rittig 2024 (RSC 403), Leenhouts 2024, Yeghyan 2026,
Schweidtmann 2020 and 2023, Saldana 2013, Paricaud 2019, Gaudin 2015, Li 2021, Oh 2022, Ogbomo 2026, Rebello 2026, Qu 2026, and the
other Energy & Fuels and Fuel papers (paywalled or not found open). Three subagent fetches were also denied by the permission
classifier, so they were not run. Statements in sections 1 to 3 that rest on those abstracts, notably the "positive reports use random
splits" reading and the MolPool oxygenate-only gain, remain unverified inferences.

## 5. Limits of this review
- Only HANNA (Chem. Sci. text), the Hall dissertation and Boehm 2022 were read in full; see section 4b. Consensus returns abstracts for
  the rest, so the claim that positive reports use random splits is an inference from abstracts and standard practice, to be checked per paper.
- Preprints (KMM-CP, Laghuvarapu 2026; Lee 2026) are not peer reviewed.
- Several 2026 items (Dong, Watanasiri, Yeghyan, Ogbomo, Rebello, Qu, Song) have 0 to 5 citations; treat as emerging.

## 6. References (as returned by Consensus, 2026-11-21)
- Bhattacharyya et al. (2026). Group-weighted conformal prediction. Electron. J. Stat.
- Boehm et al. (2022). Blend prediction model for the freeze point of jet fuel range hydrocarbons. Energy Fuels.
- Boehm et al. (2022). Error quantification of the Arrhenius blending rule for viscosity of hydrocarbon mixtures. Front. Energy Res.
- Coburn et al. (2022). Determination of a freeze point blend prediction model. AIAA SciTech.
- Dobbelaere et al. (2022). ML for physicochemical property prediction of complex hydrocarbon mixtures. Ind. Eng. Chem. Res.
- Dong et al. (2026). SISSO-based interpretable ML for jet fuel properties from GCxGC. Adv. Theory Simul.
- Gaudin et al. (2015). Combining mixing rules with QSPR models for flash points of binary mixtures. Fire Saf. J.
- Gibbs et al. (2025). Conformal prediction with conditional guarantees. J. R. Stat. Soc. B.
- Grewe et al. (2026). Physics-informed ML under small-data constraints. LNCS.
- Hall et al. (2021). Probabilistic mean QSPR modeling of jet fuel properties. Energy Fuels.
- Hall et al. (2021). Probabilistic ML models for density of conventional and synthetic jet fuels (GPR). Energy Fuels.
- Hoffmann et al. (2026). Thermodynamically consistent ML model for excess Gibbs energy. Nat. Commun.
- Jirasek and Hasse (2023). Combining ML with physical knowledge in thermodynamic modeling of fluid mixtures. Annu. Rev. Chem. Biomol. Eng.
- Laghuvarapu et al. (2026). KMM-CP. arXiv (preprint).
- Larsson et al. (2023). ML for fuel property predictions: multi-task and transfer learning. SAE.
- Leenhouts et al. (2024). Property prediction of fuel mixtures using pooled GNNs. Fuel.
- Li et al. (2021). ML-QSPR for fuel physicochemical properties of multiple fuel types. Fuel.
- Liaw et al. (2010). Prediction of miscible mixtures flash point from UNIFAC. Fluid Phase Equilib.
- Mehta et al. (2023). Chemical functional group composition of jet fuels using GCxGC. Energy Fuels.
- Ogbomo et al. (2026). Benchmarking ML and MD for thermophysical properties of organic liquids. Digital Discovery.
- Oh et al. (2022). Deep neural networks for assessing sustainable jet fuels from GCxGC. AIAA SciTech.
- Paricaud et al. (2019). Flash points of multicomponent systems. Fuel.
- Qu et al. (2026). 3D molecular representation learning for organic mixtures. Chem. Eng. Sci.
- Rebello et al. (2026). ExPUFFIN. Chem. Eng. J. Adv.
- Rittig et al. (2023). Gibbs-Duhem-informed neural networks. Digital Discovery. Rittig et al. (2024). Thermodynamics-consistent GNNs. Chem. Sci.
- Saldana et al. (2013). Prediction of flash points for fuel mixtures using ML and a novel equation. Energy Fuels.
- Schweidtmann et al. (2020). GNNs for prediction of fuel ignition quality. Energy Fuels.
- Shao et al. (2024). Prediction of critical jet fuel properties with multiple ML models. Chem. Eng. Sci.
- Shi et al. (2017). Composition-property relationship of aviation fuel from GCxGC. Fuel.
- Song et al. (2026). QSPR flash point of binary miscible mixtures, DE-SVR. J. Mol. Graph. Model.
- Specht et al. (2024). HANNA. Chem. Sci.
- Tibshirani et al. (2019). Conformal prediction under covariate shift. arXiv.
- Watanasiri et al. (2026). Evaluation of methods for specification properties of aviation fuels. Fuel.
- Winter et al. (2023). SPT-NRTL. Fluid Phase Equilib.
- Yang et al. (2025). TeNNet-SAC. J. Chem. Inf. Model.
- Yeghyan et al. (2026). Mixture-QSPR approaches for jet fuel-component blends: a review. Mol. Inform.
- Zhang et al. (2026). Physics-informed GNN for multicomponent activity coefficients (PEGAT). Syst. Control Trans.
