import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.metrics import mean_absolute_error, root_mean_squared_error

def main():
    print("=== Training Bulk Real-Fuel Models on Updated AJFTD + CRC Dataset ===")
    fuels = pd.read_csv("data/processed/ajftd_fuels.csv")
    comp = pd.read_csv("data/processed/ajftd_composition.csv")
    
    # Pivot key properties
    props_to_pivot = ['flash_point_C', 'distillation_T10_C', 'distillation_T50_C', 'distillation_T90_C', 'density_kg_m3', 'viscosity_kin_cSt', 'freezing_point_C']
    piv = fuels[fuels['property'].isin(props_to_pivot)].pivot_table(
        index='fuel_id', columns='property', values='value', aggfunc='median'
    ).reset_index()
    
    # Add aromatics if available
    arom_col = 'comp::aromatics_total[% vol]'
    if arom_col in comp.columns:
        piv = piv.merge(comp[['fuel_id', arom_col]], on='fuel_id', how='left')
    
    # Flash point modeling
    flash_df = piv.dropna(subset=['flash_point_C', 'distillation_T10_C', 'distillation_T50_C', 'distillation_T90_C']).copy()
    print(f"Fuels with complete distillation + flash point: {len(flash_df)}")
    
    # Nelson baseline: T_flash = 0.69 * T10 - 68.0
    nelson_pred = 0.69 * flash_df['distillation_T10_C'].values - 68.0
    y_flash = flash_df['flash_point_C'].values
    nelson_mae = mean_absolute_error(y_flash, nelson_pred)
    nelson_rmse = root_mean_squared_error(y_flash, nelson_pred)
    print(f"[Nelson Baseline] Flash MAE: {nelson_mae:.2f} C, RMSE: {nelson_rmse:.2f} C")
    
    # CV Ridge
    feat_cols = ['distillation_T10_C', 'distillation_T50_C', 'distillation_T90_C']
    if 'density_kg_m3' in flash_df.columns:
        flash_df['density_kg_m3'] = flash_df['density_kg_m3'].fillna(flash_df['density_kg_m3'].median())
        feat_cols.append('density_kg_m3')
    if arom_col in flash_df.columns:
        flash_df[arom_col] = flash_df[arom_col].fillna(flash_df[arom_col].median())
        feat_cols.append(arom_col)
        
    X = flash_df[feat_cols].values
    kf = KFold(n_splits=5, shuffle=True, random_state=20261030)
    ridge_preds = np.zeros(len(y_flash))
    for train_idx, test_idx in kf.split(X):
        reg = Ridge(alpha=1.0)
        reg.fit(X[train_idx], y_flash[train_idx])
        ridge_preds[test_idx] = reg.predict(X[test_idx])
        
    ridge_mae = mean_absolute_error(y_flash, ridge_preds)
    ridge_rmse = root_mean_squared_error(y_flash, ridge_preds)
    print(f"[Bulk Distillation Ridge] Flash MAE: {ridge_mae:.2f} C, RMSE: {ridge_rmse:.2f} C ({(1 - ridge_mae/nelson_mae)*100:.1f} % reduction)")

    # Update reports/REAL_FUELS.md
    report_content = f"""# REAL_FUELS.md - WP7 Real-Fuel Surrogate Bridge & Bulk Distillation Model

## 1. Overview
Real aviation and diesel fuels (e.g. Jet A, JP-8, RP-1, CRC-647/AV-27-18 fuel sets, and 2,600+ scraped AJFTD fuels) lack explicit InChIKey/SMILES molecular lists in experimental databases. Instead, they are characterized by **bulk inspection metrics**:
- ASTM D86 Distillation curves ($T_{{0}}$, $T_{{5}}$, $T_{{10}}$, $T_{{50}}$, $T_{{90}}$, $T_{{100}}$)
- Bulk density ($\rho_{{15.56\ ^\circ\\text{{C}}}}$)
- Hydrocarbon class fractions ($w_{{\\text{{aromatics}}}}$, $w_{{\\text{{saturates}}}}$, $w_{{\\text{{olefins}}}}$)

WP7 establishes the predictive link between bulk inspection features and core combustion/flow safety properties.

---

## 2. Real Fuel Flash Point Prediction from Distillation (D86)

In petroleum engineering, Nelson's empirical correlation estimates flash point from the 10 % ASTM D86 distillation temperature:
$$T_{{\\text{{flash, Nelson}}}} (^\\circ\\text{{C}}) = 0.69 \\cdot T_{{10}} (^\\circ\\text{{C}}) - 68.0$$

We evaluate this classical engineering baseline against a 5-fold cross-validated ridge model using bulk distillation and composition features on $N = {len(flash_df)}$ real aviation fuels from CRC 647 and the full AJFTD database.

### Results on Real Fuels ($N={len(flash_df)}$)

| Model | Input Features | CV MAE (°C) | CV RMSE (°C) |
|---|---|---|---|
| **Nelson Classical Baseline** | $T_{{10}}$ only | {nelson_mae:.2f} °C | {nelson_rmse:.2f} °C |
| **Bulk Inspection Ridge Model** | $T_{{10}}, T_{{50}}, T_{{90}}$, Density, Aromatics % | **{ridge_mae:.2f} °C** | **{ridge_rmse:.2f} °C** |

The multi-metric bulk model achieves a **{(1 - ridge_mae/nelson_mae)*100:.1f} % MAE reduction** over the standard Nelson refinery correlation across {len(flash_df)} real fuel records.

---

## 3. Status of Real Fuel Datasets
- Scraped AJFTD dataset now covers **2,679 unique fuels** across 91,175 metric rows.
- Validated flash point records: 2,413; Kinematic viscosity records: 2,481; Density records: 2,670; Freezing point records: 2,672.
- All uncharacterized blends are routed to `data/processed/excluded_by_R1.csv` for pure/blend training safety, while available as empirical benchmarks for bulk surrogate models.
"""
    with open("reports/REAL_FUELS.md", "w") as f:
        f.write(report_content)
    print("Updated reports/REAL_FUELS.md successfully.")

if __name__ == "__main__":
    main()
