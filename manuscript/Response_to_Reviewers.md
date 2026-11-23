# Point-by-Point Response to Reviewers

**Manuscript Title**: *Predicting Hydrocarbon Fuel Blend Properties from Chemical Composition: A Physics-First Residual Learning Framework with Calibrated Conformal Safety Margins*  
**Authors**: Kahaan Shah, Namyaa Parmar, Arijit Ganguli, Jayendra Bhalodiya, A. S. K. Sinha  
**Affiliation**: School of Engineering and Applied Science, Ahmedabad University  
**Corresponding Author**: `kahaan.s2@ahduni.edu.in`

---

We thank the Editor and the five Reviewers for their thoughtful, rigorous, and constructive feedback. Below, we provide our detailed point-by-point response outlining the revisions made to the manuscript. All updates have been incorporated into [`main.tex`](./main.tex) and recompiled in the final PDF [`main.pdf`](./main.pdf).

---

### Response to the Editor-in-Chief & Journal-Fit Reviewer

**Comment 1 (SAF Transition Pathways Context)**:  
*“The Introduction mentions SAF (HEFA, ATJ) in §1, but §5 could briefly outline how the 3 certified blends compare with approved ASTM D7566 annexes.”*

**Response**:  
We appreciate this suggestion. In the revised manuscript (§5.3), we have expanded the discussion on aviation fuel blending streams. Specifically, we discuss how the P0y and P0y-b candidate blends—which utilize synthetic paraffinic kerosene (SPK) backbones augmented by decalins and aromatics—relate to synthetic kerosene blending components (such as HEFA-SPK and ATJ-SPK) cleared under ASTM D7566.

---

### Response to Reviewer 1 (Methodology & Statistical Machine Learning)

**Comment 1 (Binary-to-Multicomponent Conformal Transfer)**:  
*“Conformal prediction guarantees $1-\alpha$ coverage only when calibration and test points are exchangeable. Calibrating margins on binary lab systems and applying them to 6-component blends is an extrapolation transfer. Explicitly defend why molecular interaction in non-polar liquid hydrocarbons is dominated by pairwise potentials, justifying why binary-calibrated margins provide conservative bounds for multicomponent fuels.”*

**Response**:  
We have added a rigorous physical and thermodynamic defense in §5.1 and §5.4. In liquid-state perturbation theory and statistical mechanics, non-polar and weakly polar hydrocarbon mixtures (normal/iso-alkanes, cycloalkanes, and alkylbenzenes) interact predominantly through London dispersion forces. Two-body pairwise dispersion potentials ($r^{-6}$) account for over 95\% of the cohesive energy. Three-body and higher-order multi-body dispersion energies (such as Axilrod--Teller triple-dipole interactions) decay asymptotically as $r^{-9}$ or faster, contributing under $2\text{--}5\%$ of total excess mixing properties. Therefore, binary pairwise interactions capture the vast majority of non-ideal behavior, validating the transfer of binary-calibrated conformal bounds to 6-component blends.

**Comment 2 (Operational Behavior of Infinite Margins)**:  
*“The manuscript notes that cold viscosity has a margin of $\infty$, but does not state how the simplex optimizer handles an active limit when $\text{Margin} = \infty$.”*

**Response**:  
We have clarified this operational protocol in §3.5. When an optimization constraint involves a property whose conformal margin is infinite (such as real-fuel viscosity below $263\text{ K}$, where models exhibit no proven skill over a constant baseline), the inverse design engine strictly refuses to relax or drop the constraint. Instead, it immediately halts execution and raises a `ModelUnavailable` exception, ensuring that uncertifiable fuel formulations can never be falsely emitted.

---

### Response to Reviewer 2 (Chemical Thermodynamics & Petroleum Blending)

**Comment 1 (NIST Antoine Out-of-Range Extrapolation & Re-anchoring)**:  
*“Explain that standard WebBook Antoine constants for $\text{C}_{10+}$ alkanes were fitted to high-temperature boiling data ($95\text{--}190^\circ\text{C}$) and produce substantial errors when extrapolated down to $0\text{--}60^\circ\text{C}$ flash point temperatures, necessitating re-fitting to measured vapor pressures.”*

**Response**:  
We have incorporated this exact thermodynamic explanation in §2.2.3. Historically, standard Antoine parameters in the NIST WebBook for $\text{C}_{10}\text{--}\text{C}_{16}$ hydrocarbons were derived from reduced-pressure or atmospheric boiling measurements in the $95\text{--}190^\circ\text{C}$ range. Extrapolating these parameters downward into ambient and sub-ambient flash regimes ($0\text{--}60^\circ\text{C}$) leads to severe vapor pressure underestimation, which previously caused flash point under-prediction biases as severe as $-21.7^\circ\text{C}$. In this work, Antoine constants were systematically re-anchored and refitted against experimental vapor pressure measurements spanning the $0\text{--}100^\circ\text{C}$ window, reducing real-fuel flash bias to $+0.059^\circ\text{C}$.

---

### Response to Reviewer 3 (Aviation Fuels & Industrial Standards Specialist)

**Comment 1 (Dicycloparaffin Concentration in Jet A-1 Design)**:  
*“The certified P0y-b Jet A-1 candidate blend contains $37.8\text{ wt}\%$ dicycloparaffins (Table 4). In conventional aviation kerosene, dicycloparaffins rarely exceed $15\text{--}20\%$. In §5.3, discuss whether high dicycloparaffin fractions could influence seal compatibility (e.g., nitrile O-ring swelling) or smoke point, and recommend compositional penalty functions for future iterations.”*

**Response**:  
We have added an in-depth discussion of this chemical reality in §5.3. While bicyclohexyl and decalin isomers provide exceptional volumetric energy density and depress freezing point down to $-61.0^\circ\text{C}$, conventional commercial kerosenes typically contain $\le 20\text{ wt}\%$ dicycloparaffins. Abnormally elevated dicycloparaffin concentrations may alter elastomeric seal swelling (e.g., nitrile O-rings in legacy aircraft) and slightly reduce smoke point. We explicitly recommend that future iterations of the inverse design optimizer implement compositional penalty constraints (e.g., bounding dicycloparaffins $\le 30\text{ wt}\%$ in accordance with AJFTD distributions).

**Comment 2 (Distillation Offset Model Specification)**:  
*“The paper references distillation limits in ASTM D1655 but does not specify which distillation model was used in the optimization.”*

**Response**:  
We have updated §3.5 to state explicitly that distillation boiling percentiles ($T_{10}, T_{50}, T_{90}$) are evaluated using the empirical Engler D86 model, and note that the alternative Rayleigh distillation model was evaluated during development but rejected because its $+11.6^\circ\text{C}$ bias on $T_{10}$ failed pre-declared gating criteria.

---

### Response to the Devil's Advocate Reviewer

**Comment 1 (The Razor-Thin Slack Vulnerability)**:  
*“The three certified blends have minimum slacks ranging from $0.003$ to $0.013$ (Table 4, Figure 4). The optimizer has converged exactly onto the conformal margin boundary. The authors must explicitly state in §5.3 that the certified blends represent the theoretical mathematical frontier of compliance, and that practical industrial implementation requires introducing a deliberate operational buffer ($\Delta_{\text{safety}} \ge 1.0^\circ\text{C}$) beyond the statistical margin.”*

**Response**:  
We fully agree with this crucial engineering insight. In §5.3, we have explicitly noted that the tight slacks ($0.003$ to $0.013$) demonstrate convergence onto the active mathematical constraint boundary. We have added a clear operational recommendation: in physical refinery and blending operations, an engineering safety buffer (e.g., $\Delta_{\text{oper}} \ge 1.0^\circ\text{C}$ on flash point or $0.5\text{ cSt}$ on viscosity) must be superimposed upon the conformal margin to account for physical metering tolerances, tank stratification, and ambient temperature swings.

**Comment 2 (Reconstructed CRC Compositions)**:  
*“Compositions for the 57 CRC 647 jet fuels were reconstructed from bulk hydrocarbon classes (Table G.1), not direct molecular $\text{GC}\times\text{GC}$ speciation. Acknowledge in §5.4 that the density and flash errors observed on the CRC 647 cohort represent an upper bound on error.”*

**Response**:  
We have updated §5.4 (Limitation 4) to explicitly state that the 57 CRC 647 jet fuel compositions were reconstructed from bulk hydrocarbon classes rather than full molecular $\text{GC}\times\text{GC}$ speciation. Consequently, reported errors on the CRC cohort represent a compound upper bound that combines surrogate reconstruction uncertainty with intrinsic thermodynamic model error.
