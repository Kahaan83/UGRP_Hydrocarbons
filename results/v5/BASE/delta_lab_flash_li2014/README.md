# delta_lab_flash_li2014: more lab flash systems (n-alkane mixtures)

Why: BASE Task 10 found that the lab one-sided flash margin of 3.70 C is the group score of one system (mesitylene) out of 31, and that it falls to 0.31 C without it. More measured hydrocarbon mixture systems are the honest way to settle Q2, instead of choosing between 3.70 and 0.31 C after seeing them.

## Sources (numbers transcribed by hand; no PDF, no quoted text in the repo)

| Source | doi | Tables used | Page(s) | Licence (first page) |
|---|---|---|---|---|
| Li, Duan, Sun, Yan, Adv. Mater. Sci. Eng. 2014, 239023 | 10.1155/2014/239023 | 5, 6, 7 (binary), 8, 9, 10 (ternary), column "Exp." | 3-4 | Creative Commons Attribution (open access) |
| Ha, Lee, J. Korean Inst. Gas 20(6), 31-36 (2016) | 10.7842/kigas.2016.20.6.31 | 2 (n-nonane + n-decane + n-dodecane), column "Exp." | 34 | not an open licence: copyright held by the Korean Institute of Gas; only the measured numbers are transcribed, with the citation |

Transcribed: only measured flash points of mixtures, mole-fraction basis (printed as "Mole fraction"). Not transcribed: calculated columns (model, Raoult, empirical), pure-component rows (they are not mixtures), Table 3 (binary mixtures at a 1:1 volume ratio, no mole fractions printed). Li et al. use a closed cup tester (WFY-111C); Ha and Lee use a Seta flash closed cup tester (ASTM D3278). Mole fractions as printed sum to 0.99-1.01 in a few rows; the raw CSV keeps the printed values and the loader renormalises. The raw file is `data/raw/flash_lab_li2014/li2014_flash_mixtures.csv`; the loader is `src/fuelblend/data/loaders_li2014.py`.

## De-duplication against `results/v4/F/flash_lab_loso.csv`

Rule: same `system_id` (md5[:12] of the sorted InChIKey set), component-aligned x within 0.005 and flash within 0.5 C = duplicate, counted and dropped.

* Transcribed mixture points: 33 (24 Li et al. + 9 Ha and Lee), in 6 n-alkane systems.
* Duplicates dropped: 0 (checked before any score was computed).
* New points kept: 33. New systems: 5; one system (octane + dodecane, `314cc4c34772`) already exists in the lab file with a single point (x_dodecane 0.927), so 4 new points join an existing system.
* The ternary nonane + decane + dodecane has 13 new points (4 Li et al. Table 10, 9 Ha and Lee Table 2), all with distinct compositions from each other and from the lab file; they share one system_id and count as one group.
* Lab file before: 126 points, 31 systems. After: 159 points, 36 systems.

## Declared readings (committed before any score on the new points is computed)

New points are labelled `dev` (first read), never `held_out`. Predictions use the same per-row route as `delta_flash_inputs.lab_predictions` (ideal Le Chatelier through `flash_point_unifac_detail`, `use_unifac=False`), for both the designated variant and the LFL overlay.

* **R-LI-1.** On the new points, the overlay flash has |bias| <= 2 C and an over-prediction rate (error > 5 C) <= 10 %.
* **R-LI-2.** The lab one-sided margin, `one_sided_margin(direction="over", alpha=0.10)`, over the old plus new systems, computed with and without the mesitylene system, differs by less than 0.5 C.
* If R-LI-2 holds, that margin is proposed for `configs/uncertainty_v5.yaml`. Writing that file stays the user's decision (Q2); this task does not write it.
* If no new point survives de-duplication, both readings are recorded as "not evaluable: no new points" (not the case here: 33 points survive).

## Result (script `scripts/v5/base/delta_lab_flash_li2014.py`; files `delta.csv`, `lab_margin.csv`, `points.csv`, `readings.json`)

Note on commit order: the loader and its fixture test were staged by mistake in the same commit as this README's declared readings; they contain no scores and no result existed at that point.

**R-LI-1: fails.** On the 33 new points the overlay flash has bias -3.54 C (|bias| > 2 C), MAE 4.22 C. The over-prediction rate (> 5 C) is 0 % (passes its half), because the model under-predicts: the five binary systems sit at -4.7 to -6.0 C mean error (safe side), the ternary nonane + decane + dodecane (13 points) at -0.6 C mean, with a maximum over-prediction of +2.2 C. The overlay and the designated variant give the same values for these n-alkane points (LFLs of the n-alkanes were already measured/identical), so the overlay changes nothing here.

**R-LI-2: fails.** Lab margin (overlay, one-sided over, alpha 0.10): old 31 systems 3.70 C; old without mesitylene 0.31 C; old plus new 36 systems 3.70 C; old plus new without mesitylene 2.21 C. The difference with versus without mesitylene is 1.49 C (> 0.5 C), so no margin is proposed for `configs/uncertainty_v5.yaml`. The margin with mesitylene removed rises from 0.31 to 2.21 C because the new ternary system (max over-prediction +2.2 C) becomes the 90 % quantile group; the margin is still determined by a single system, now two different ones depending on the scope. The designated variant has identical margins.

Honest reading: the extra systems do not settle Q2. The lab one-sided margin remains sensitive to which single system sits at the quantile; the new n-alkane data say the route is not unsafe on binary n-alkane mixtures (error on the conservative side), and they add one ternary system with small positive error. The new points are `dev` (first read), not `held_out`. Q2 stays the user's decision.
