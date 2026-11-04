"""T24 loader: GCV pure-compound from chemicals (tier A stoichiometric) + Dulong fallback (tier C).
Lazy: uses compounds.csv CAS lookup; writes data/processed/gcv_MJ_kg.csv.
No modeling, no build_unified.py change."""
import csv, re
from chemicals import search_chemical, Hfl, Hfg, HHV_stoichiometry, combustion_stoichiometry, HHV_modified_Dulong

def parse_formula(s):
    elems = re.findall(r'([A-Z][a-z]?)(\d*)', s)
    return {el: int(num) if num else 1 for el, num in elems}

in_path = 'fuelblend_data/data/unified/compounds.csv'
out_path = 'fuelblend_data/data/processed/gcv_MJ_kg.csv'

rows = []
with open(in_path, newline='') as f:
    for r in csv.DictReader(f):
        inchikey = r['inchikey']
        smiles = r['smiles']
        cas = r.get('cas', '').strip() or ''
        # Try CAS first, then SMILES
        meta = None
        if cas:
            try: meta = search_chemical(cas)
            except Exception: pass
        if meta is None and smiles:
            try: meta = search_chemical(smiles)
            except Exception: pass
        if meta is None:
            # try common_name / inchikey not supported; skip to Dulong with manual formula from SMILES impossible
            # For minimal load, use a hardcoded formula lookup for a few known compounds? No — skip.
            # We'll attempt formula from SMILES using rdkit only if available; skip for now.
            continue
        # Use liquid formation enthalpy Hfl for 25 C liquid fuel basis
        try:
            H = Hfl(meta.CAS)
        except Exception:
            try:
                H = Hfg(meta.CAS)
            except Exception:
                H = None
        if H is not None:
            stoich = combustion_stoichiometry(parse_formula(meta.formula))
            hhv = HHV_stoichiometry(stoich, H)
            value = abs(hhv) / meta.MW / 1000.0  # J/mol / (g/mol) /1000 = MJ/kg
            tier = 'A'
            method = 'stoichiometric'
            derived = False
        else:
            # fallback: Dulong from formula mass fractions
            f = parse_formula(meta.formula)
            mw = meta.MW
            # simplified mass fractions of C,H,O,S; ignore others for Dulong
            # compute directly using atomic weights from simple table
            aw = {'C':12.011,'H':1.008,'O':15.999,'S':32.06,'N':14.007}
            mass_frac = {el: f[el]*aw.get(el,0)/mw for el in f}
            # Dulong only needs C,H,O,S; others ignored in formula but fractions must sum ~1
            val = HHV_modified_Dulong({'C':mass_frac.get('C',0), 'H':mass_frac.get('H',0), 'O':mass_frac.get('O',0), 'S':mass_frac.get('S',0)})
            value = abs(val) / 1000.0
            tier = 'C'
            method = 'modified_Dulong'
            derived = True
        rows.append({
            'inchikey': inchikey,
            'smiles': smiles,
            'property': 'gcv_MJ_kg',
            'value': f"{value:.4f}",
            'unit': 'MJ/kg',
            'method': method,
            'tier': tier,
            'derived': str(derived),
            'T_ref_K': '298.15',
            'ref': f"chemicals_1.5.2_{method}; meta.CAS={meta.CAS}",
        })

with open(out_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['inchikey','smiles','property','value','unit','method','tier','derived','T_ref_K','ref'])
    writer.writeheader()
    writer.writerows(rows)

print('Wrote', len(rows), 'rows to', out_path)
# sanity
for r in rows:
    try:
        v = float(r['value'])
    except: continue
    # check well-known
    if r['inchikey'].startswith('DRDJHWIZZEZHHV') or 'octane' in r['smiles'] or r['smiles']=='CCCCCCCC':
        print('SANITY octane-like:', r['smiles'], 'GCV=', r['value'], 'tier=', r['tier'])
    if r['smiles']=='c1ccccc1' or r['inchikey']=='UHOVQNZJYSORNB-UHFFFAOYSA-N':
        print('SANITY benzene:', r['smiles'], 'GCV=', r['value'], 'tier=', r['tier'])
