# Literature synthesis: gaps and fixes for v6

Status: 2026-11-21. Scope: whether the literature shows anything v6 missed, and what to fix in `docs/v6_summary.md`. Sources: Consensus
searches run 2026-11-21 (freezing point, mixture ML, flash point, small-data pretraining, conformal prediction, melting point) plus the
verified list in `docs/literature_review_v6.md`. Only abstracts were read, not full texts, so every claim below is limited to what an
abstract states. No v6 number is changed. Numbers quoted from papers are as reported there and are not comparable with our grouped CV
unless stated. One figure (binomial tail, section 3.4) is derived here and marked as such.

## 1. What the literature says about each v6 outcome

| v6 outcome | What the literature adds | Reading |
|---|---|---|
| Freezing: dev gain did not replicate on holdout | Hall et al. report that for freezing point and viscosity, fuels must be in training, and the effect is strong ([M-QSPR](https://consensus.app/papers/details/w:ce2baf7e-35d8-4cc9-beb3-c728f38b6eab/?utm_source=claude_desktop)). Bell et al. get MAE 4.4 C from single-species control curves and note that models fitted on conventional fuels fail on dissimilar compositions, one by 27.9 C ([Bell 2025](https://consensus.app/papers/details/w:ea051352-6bd0-4588-b8bf-5ba5ae4b9720/?utm_source=claude_desktop)). | Our failure is the expected regime: small, family-clustered fuel sets. Not a bug. |
| Tm head loses to v5 fusion estimator (30.4 vs 28.4 K) | Reported organic melting-point MAE is about 21 to 27 K on 100k+ compound sets ([Huang 2025](https://consensus.app/papers/details/w:7fa8f0ff-1180-4112-8b20-6e4576770fcd/?utm_source=claude_desktop), [Dubonos 2026](https://consensus.app/papers/details/w:10c57f67-84c7-4749-a5a8-eab73f6e58a5/?utm_source=claude_desktop)); experimental noise near 32 C is argued to be the floor for drug-like sets ([Tetko 2016](https://consensus.app/papers/details/w:df91db18-607f-4e7c-af2f-a5cd6b98a41c/?utm_source=claude_desktop)). | 30.4 K on 1,208 compounds is at the published scale. Do not expect a GNN Tm to beat measured or thermodynamic Tm. |
| Track B flash 19.1 C | Graph models report MAE 17.8 K on 10,575 molecules ([Sun 2020](https://consensus.app/papers/details/w:ee0de57d-a6d7-4fe5-bf5e-06e17eb3543e/?utm_source=claude_desktop)) and 14.42 C on 14,696 ([Lee 2026](https://consensus.app/papers/details/w:29eecbdf-8c75-4c55-bcaa-a2040078480a/?utm_source=claude_desktop)). | Same order. Our CV is scaffold-grouped and smaller, so no claim of being better or worse. |
| Mixture heads fail replication | Boehm et al. conclude that blending-rule error is small relative to error from incomplete composition and component property data ([Boehm 2024](https://consensus.app/papers/details/w:0ef9d397-c5fd-4775-8f22-09b3ec607c4c/?utm_source=claude_desktop)). Dobbelaere et al. find mixture ML helps where mixing is non-linear ([Dobbelaere 2022](https://consensus.app/papers/details/w:e798b7e9-ddd3-41e2-a1a8-68f811f0ac7e/?utm_source=claude_desktop)). A graph-model OOD benchmark finds large generalisation gaps versus random splits ([Omee 2024](https://consensus.app/papers/details/w:77f59d8f-ef04-4e74-aba7-6528789dad2d/?utm_source=claude_desktop)). | Consistent with our negative results. The ceiling on mixtures is data and composition quality, not architecture. |
| Flash oxygenate arm: points favour MolPool, CI crosses 0 | Leenhouts et al. report MolPool gains over blending laws for DCN, flash point and viscosity, largest when oxygenates are present ([Leenhouts, Fuel](https://consensus.app/papers/details/w:f0fc82d4-e0ae-48d5-9f78-06c9021b29b0/?utm_source=claude_desktop)). | Our dev and locked point estimates agree in direction with the original claim. The gate was not met, so the correct wording is "direction consistent, not established". |

## 2. Gaps: things the literature suggests v6 did not do

Ordered by expected value per effort. Each is a hypothesis to gate under the existing protocol, not a result.

1. **No classical baseline on Track B beyond a constant (and ridge for YSI).** The summary itself calls the constant "a weak bar".
   Foundation-model work states that D-MPNNs often fail to beat classical methods on small real benchmarks, and that Random Forest on
   descriptors is a strong comparator ([CheMeleon, Burns 2026](https://consensus.app/papers/details/w:162c27f6-8d23-435f-b916-c890e816ca03/?utm_source=claude_desktop)).
   Fix: add RF or GBM on RDKit descriptors, same scaffold folds, all 10 targets. Cheap and decisive for whether the GNN earns its place.
2. **Hand-written D-MPNN not cross-checked against Chemprop.** Chemprop v2 is the reference implementation
   ([Heid 2023](https://consensus.app/papers/details/w:6e0f07bb-4dc5-4694-993c-f17498a6cbc1/?utm_source=claude_desktop),
   [Graff 2025](https://consensus.app/papers/details/w:6e7f5bf0-a650-4d81-b1ac-f07ce48d4d3c/?utm_source=claude_desktop)).
   Fix: one parity run of `src/fuelblend/ml/dmpnn.py` against Chemprop on the same folds for two targets. If they agree, the encoder
   is trustworthy; if not, Track C inherits the defect.
3. **Encoder pretrained only on 1,208 labelled compounds.** Descriptor-based pretraining
   ([CheMeleon](https://consensus.app/papers/details/w:162c27f6-8d23-435f-b916-c890e816ca03/?utm_source=claude_desktop)) and
   QM-to-experiment transfer ([Vermeire 2021](https://consensus.app/papers/details/w:0cc3a6b2-eeb7-4927-a165-0074318905d6/?utm_source=claude_desktop))
   report gains in small-data and out-of-sample settings. QM descriptors help D-MPNN mostly on small sets, and only when they correlate
   with the target ([Li 2024](https://consensus.app/papers/details/w:26f36bfb-1095-47c1-9ce3-1989e677dfa3/?utm_source=claude_desktop)).
   Fix: test a descriptor-pretrained encoder as one extra arm. Check licences before using any released checkpoint.
4. **Freezing: restrict to crucial components.** Bell et al. and Boehm et al. both say freezing is set by the first species to freeze,
   so a full-composition average is the wrong structure ([Bell 2025](https://consensus.app/papers/details/w:ea051352-6bd0-4588-b8bf-5ba5ae4b9720/?utm_source=claude_desktop),
   [Boehm 2022](https://consensus.app/papers/details/w:38b20121-671f-468f-a6be-e74e86d0eefc/?utm_source=claude_desktop)).
   This matches our bias of +24 C for the mass-weighted Tm average. Check that baseline (d) was run in its limiting-species form and
   report it beside the others; the summary lists it but quotes no number for it.
5. **Freezing: joint pure plus fuel training (M-QSPR style).** Hall et al. train one model on pure compounds and fuels together using
   a probability-weighted mean structure, with no mixing rule. Track A used a LightGBM residual on coarse composition, which is a
   different design. Fix: MolPool over the GC x GC species list (Track B encoder) trained on fuel freezing labels with Tm as auxiliary.
   Risk: 64 groups is the same order as their 82 plus 50 fuels, so expect the same instability.
6. **Flash oxygenate arm: restore the physics comparator.** Le Chatelier dropped out because Antoine and LFL coverage was zero on
   oxygenate rows. GNNs now predict vapour pressure and flash point jointly ([Lee 2026](https://consensus.app/papers/details/w:29eecbdf-8c75-4c55-bcaa-a2040078480a/?utm_source=claude_desktop)),
   and thermodynamically consistent activity-coefficient GNNs exist ([Rittig 2024, GE-GNN](https://consensus.app/papers/details/w:9394a059-8c70-479a-be59-5fd3259234b6/?utm_source=claude_desktop)).
   Fix: predicted Antoine and LFL inputs plus GE-GNN gamma in a Liaw-type flash model as a physics comparator. Otherwise "MolPool
   beats descriptor GBM" is compared only against a statistical model.
7. **Stronger descriptor comparator for Track C.** Tabular foundation models with CheMeleon or RDKit descriptors are reported to
   beat fine-tuned GNNs at low data ([Ben Hicham 2026](https://consensus.app/papers/details/w:a3d9da01-6765-4694-bc39-3da6d41aeea3/?utm_source=claude_desktop);
   Zenodo-hosted, treat as non-peer-reviewed until checked). If the descriptor GBM was already hard to beat (4.97 vs 5.28), this is the
   next comparator to expect from a reviewer.
8. **Conformal margins.** See section 3.4.

## 3. Fixes for `docs/v6_summary.md`

3.1 **Soften the Track C verdict (section 5).** "The earlier pass is best read as chance" is stronger than the data. The replication
CI is [-1.24, 0.58], which includes both a real gain and a loss. Suggested wording: "the replication cannot distinguish the earlier
pass from chance; no skill shown over the strongest comparator."

3.2 **Add literature context to section 4.** The sentence "Every 90 % bootstrap gain interval excludes zero" sits beside a weak
constant. Add: published graph-model flash MAE is 14 to 18 K on 10k+ molecules ([Sun 2020](https://consensus.app/papers/details/w:ee0de57d-a6d7-4fe5-bf5e-06e17eb3543e/?utm_source=claude_desktop),
[Lee 2026](https://consensus.app/papers/details/w:29eecbdf-8c75-4c55-bcaa-a2040078480a/?utm_source=claude_desktop)); our 19.1 C is
the same order but on different splits, so not comparable.

3.3 **Reconcile citation metadata in `docs/literature_review_v6.md`.** Consensus lists Leenhouts et al. (Fuel) as 2024 and Hall et al.
(M-QSPR) as 2021; the review says 2025 vol 381 and 2022 vol 36. Online versus print year probably explains both. Check the DOIs
before any manuscript use. Zhang et al. (MolSets) is 2024 in both lists.

3.4 **Conformal margins (sections 3 and 4).** Two points.
- The freezing 90 % margin covered 8 of 12 holdout fuels. *Derived here:* if true coverage were 0.90 and fuels exchangeable, the chance of
  8 or fewer of 12 is about 0.026. Fuels are not exchangeable (AV-33-22 dominates dev), which is the point: the guarantee is marginal
  and the shift between dev and holdout violates exchangeability. Say so in the summary.
- The YSI margin came from CV folds that had seen the holdout compounds (already stated). Literature fixes to try, on dev only:
  group-conditional coverage ([Gibbs 2025](https://consensus.app/papers/details/w:ac8b186d-1a7a-48e8-b4ff-426919b3f8cb/?utm_source=claude_desktop)),
  covariate-shift weighting ([Tibshirani 2019](https://consensus.app/papers/details/w:9eada008-6d6d-444e-99a2-2db044b6f8c1/?utm_source=claude_desktop)),
  and the finite-group case ([Bhattacharyya 2026](https://consensus.app/papers/details/w:a932f0dd-8e2f-41f4-ba26-a7a6a8c3c004/?utm_source=claude_desktop)).
  With 12 holdout groups, none of these can be validated tightly. Say that.

3.5 **Multiplicity policy (section 7).** The summary notes that 1 of 5 heads passed without correction. Add a prospective rule for
v7: declare the family of heads and the correction (for example Holm) in the lock commit. This is methodology, not from the abstracts.

3.6 **Link the new document** from the header line, next to `docs/literature_review_v6.md`.

## 4. Limits of this review

- Abstracts only; no full texts. Methods, splits and leakage controls in the cited papers were not audited, so published MAEs are
  context, not benchmarks.
- Consensus returns relevance-ranked results, not an exhaustive search. No PRISMA-style screening was done. Several 2026 entries have
  zero citations and are preprints or very recent; treat them as leads.
- The assistant did not run WebSearch cross-checks beyond the Consensus index, so metadata mismatches in 3.3 are unresolved.
- Rate limits stopped three parallel searches; they were rerun one at a time.
