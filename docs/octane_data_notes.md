# Octane Data Quality and Inference Limits (`docs/octane_data_notes.md`)

Per Task T2.19:
> "Document data quality limits, since many octane values are inferred rather than measured."

## 1. Experimental Reality of CFR Engine Testing
Research Octane Number (RON, ASTM D2699) and Motor Octane Number (MON, ASTM D2700) are evaluated using a standardized variable-compression single-cylinder Cooperative Fuel Research (CFR) engine:
- **Reproducibility Limit:** ASTM standard test reproducibility for CFR engines is typically **0.6 to 1.0 octane units**.
- **Scale Extrapolation Artifacts:** The classical octane scale is defined by primary reference fuels: $n$-heptane ($	ext{ON} = 0$) and isooctane ($	ext{ON} = 100$). For pure hydrocarbons with octane $> 100$ (e.g., toluene $	ext{RON} pprox 120$, benzene $	ext{RON} pprox 115$), values are determined by extrapolation or by adding tetraethyl lead (TEL) equivalents, introducing systematic scaling deviations between literature sources.

## 2. Inferred vs Measured Values in the Dataset
In `data/processed/measurements_octane.csv`, every single observation is explicitly labeled in `octane_data_type`:
1. **`pure` (361 compounds, Tier A):** Directly measured experimental octane numbers from vom Lehn et al. (2021) and CombustDB (Kuzhagaliyeva et al. 2022).
2. **`blending_value` (85 compounds, Tier C):** Inferred / neural-network predicted blending octane numbers (ANN predictions from vom Lehn 2021 for complex branched alkanes). Marked with `quality_flag = 'suspect'`.
3. **`measured_mixture` (292 fuels, Tier A):** Laboratory CFR engine measurements of primary reference fuel blends (PRF), toluene reference fuels (TRF), and Fuels for Advanced Combustion Engines (FACE gasolines) from the Leenhouts/Published MON repository.

## 3. Practical Guidance for Machine Learning Models
- Models training on octane numbers should use `quality_flag == 'ok'` to restrict training to true CFR engine experimental measurements.
- Octane sensitivity ($	ext{RON} - 	ext{MON}$) reflects fuel aromatic and olefin content; paraffin-rich fuels exhibit near-zero sensitivity ($	ext{RON} pprox 	ext{MON}$), while aromatic fuels exhibit high sensitivity ($\Delta 	ext{ON} = 8\dots 15$).
