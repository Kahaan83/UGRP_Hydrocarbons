"""Extract critical properties and derive acentric factor and COSTALD V* for hydrocarbons.
Tasks T39 and T40 from DATA_ACQUISITION_TASKS.md.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd
import chemicals
from chemicals import critical, acentric, identifiers

TIER_A_METHODS = {'WEBBOOK', 'HEOS', 'IUPAC', 'MATTHEWS', 'CRC', 'PD', 'PSRK', 'PINAMARTINES'}
PRIORITY_CRIT = [
    'WEBBOOK', 'HEOS', 'IUPAC', 'MATTHEWS', 'CRC', 'PD', 'PSRK', 'PINAMARTINES',
    'YAWS', 'WILSON_JASPERSON', 'JOBACK', 'FEDORS'
]
PRIORITY_OMEGA = ['HEOS', 'PSRK', 'PD', 'YAWS', 'ACENTRIC_DEFINITION']


def get_cas_resolver():
    """Build fast index for resolving InChI / InChIKey / stereoisomers."""
    db = identifiers.pubchem_db
    db.autoload_main_db()
    prefix_index = {}
    for k, meta in db.InChI_key_index.items():
        prefix = k.split('-')[0]
        if prefix not in prefix_index:
            prefix_index[prefix] = []
        prefix_index[prefix].append(meta)
    return db, prefix_index


def resolve_cas(row: pd.Series, db, prefix_index) -> str | None:
    """Resolve CAS via inchi, smiles, inchikey or skeleton stereoisomer match."""
    ik = row['inchikey']
    smi = row['smiles']
    cas = None
    for val in [ik, row.get('inchi'), smi]:
        if pd.isna(val):
            continue
        try:
            cas = identifiers.CAS_from_any(str(val))
            if cas:
                return cas
        except Exception:
            pass

    # Try InChIKey skeleton (first 14-char block)
    if pd.notna(ik):
        prefix = str(ik).split('-')[0]
        if prefix in prefix_index:
            meta = prefix_index[prefix][0]
            return identifiers.int_to_CAS(meta.CAS)
    return None


def build_critical_constants():
    compounds_path = "data/unified/compounds.csv"
    refit_path = "data/processed/antoine_refit.csv"
    out_path = "data/processed/critical_constants_hc.csv"

    if not os.path.exists(compounds_path):
        raise FileNotFoundError(f"Missing {compounds_path}")

    df_comp = pd.read_csv(compounds_path)

    # Load refit Antoine parameters if present
    refit_map = {}
    if os.path.exists(refit_path):
        df_refit = pd.read_csv(refit_path)
        for _, r in df_refit.iterrows():
            if r['status'] == 'ok' and pd.notna(r['A']):
                refit_map[r['inchikey']] = (r['A'], r['B'], r['C'], r['T_min'], r['T_max'])

    db, prefix_index = get_cas_resolver()

    rows = []
    for _, r in df_comp.iterrows():
        ik = r['inchikey']
        smi = r['smiles']
        cas = resolve_cas(r, db, prefix_index)

        tc, pc, vc, om = np.nan, np.nan, np.nan, np.nan
        stc, spc, svc, som = '', '', '', ''

        if cas:
            tc_m = critical.Tc_methods(cas)
            for m in PRIORITY_CRIT:
                if m in tc_m:
                    try:
                        val = critical.Tc(cas, method=m)
                        # Sanity check: hydrocarbon critical temperature must be in [150, 1500] K
                        if 150.0 <= val <= 1500.0:
                            tc = val
                            stc = m
                            break
                    except Exception:
                        pass

            pc_m = critical.Pc_methods(cas)
            for m in PRIORITY_CRIT:
                if m in pc_m:
                    try:
                        pc = critical.Pc(cas, method=m)
                        spc = m
                        break
                    except Exception:
                        pass

            vc_m = critical.Vc_methods(cas)
            for m in PRIORITY_CRIT:
                if m in vc_m:
                    try:
                        vc = critical.Vc(cas, method=m)
                        svc = m
                        break
                    except Exception:
                        pass

            om_m = acentric.omega_methods(cas)
            for m in PRIORITY_OMEGA:
                if m in om_m:
                    try:
                        om = acentric.omega(cas, method=m)
                        som = m
                        break
                    except Exception:
                        pass

        # T40: derived omega from Antoine at Tr = 0.7 if available
        derived = False
        ant_params = None
        ant_src = None
        if ik in refit_map:
            ant_params = refit_map[ik]
            ant_src = 'antoine_refit'
        elif pd.notna(r.get('antoine_A')) and pd.notna(r.get('antoine_B')) and pd.notna(r.get('antoine_C')):
            ant_params = (r['antoine_A'], r['antoine_B'], r['antoine_C'], r.get('antoine_Tlo'), r.get('antoine_Thi'))
            ant_src = 'compounds_csv'

        if ant_params and pd.notna(tc) and pd.notna(pc) and pc > 0:
            A, B, C, T_lo, T_hi = ant_params
            T_07 = 0.7 * tc
            in_range = True
            if pd.notna(T_lo) and pd.notna(T_hi):
                if T_07 < T_lo - 30.0 or T_07 > T_hi + 30.0:
                    in_range = False
            if in_range:
                try:
                    log10_P_kPa = A - B / (T_07 + C)
                    Psat_Pa = (10.0 ** log10_P_kPa) * 1000.0
                    if Psat_Pa > 0:
                        om_calc = -np.log10(Psat_Pa / pc) - 1.0
                        if 0.0 < om_calc < 2.0:
                            if pd.isna(om) or som in ['YAWS', 'ACENTRIC_DEFINITION']:
                                om = om_calc
                                som = f'Psat_T07_{ant_src}'
                                derived = True
                except Exception:
                    pass

        # T40: COSTALD characteristic volume V* = Vc * (0.29056 - 0.08775 * omega)
        v_star = np.nan
        if pd.notna(vc) and pd.notna(om) and vc > 0:
            v_star = vc * (0.29056 - 0.08775 * om)
            derived = True

        # Tier assignment
        if stc in TIER_A_METHODS and spc in TIER_A_METHODS:
            tier = 'A'
        elif (stc in TIER_A_METHODS or spc in TIER_A_METHODS) or (stc and spc):
            tier = 'B'
        elif stc or spc or svc:
            tier = 'C'
        else:
            tier = 'none'

        pc_kPa = pc / 1000.0 if pd.notna(pc) else np.nan

        rows.append({
            'inchikey': ik,
            'smiles': smi,
            'formula': r['formula'],
            'mw': r['mw'],
            'tc_K': round(tc, 2) if pd.notna(tc) else np.nan,
            'pc_kPa': round(pc_kPa, 2) if pd.notna(pc_kPa) else np.nan,
            'vc_m3_mol': vc if pd.notna(vc) else np.nan,
            'omega': round(om, 4) if pd.notna(om) else np.nan,
            'v_star_m3_mol': v_star if pd.notna(v_star) else np.nan,
            'source_tc': stc,
            'source_pc': spc,
            'source_omega': som,
            'tier': tier,
            'derived': derived
        })

    out_df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    out_df.to_csv(out_path, index=False)
    print(f"Saved {len(out_df)} rows to {out_path}")
    print(f"Tc non-null: {out_df['tc_K'].notna().sum()}, Pc non-null: {out_df['pc_kPa'].notna().sum()}")
    print(f"Vc non-null: {out_df['vc_m3_mol'].notna().sum()}, Omega non-null: {out_df['omega'].notna().sum()}")
    print(f"V* non-null: {out_df['v_star_m3_mol'].notna().sum()}, Derived True: {out_df['derived'].sum()}")
    print("Tier counts:\n", out_df['tier'].value_counts().to_dict())
    return out_df


if __name__ == '__main__':
    build_critical_constants()
