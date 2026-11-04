"""Offline tests (no network). Run:  python tests/test_offline.py"""
import io, os, sys, tarfile, tempfile, subprocess
from pathlib import Path
HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import numpy as np
import xml.etree.ElementTree as ET

from common import mole_to_mass, mass_to_mole, mass_to_volume, mol_from, is_hydrocarbon, hc_class
import webbook_antoine as wb
import thermoml_extract as tx
import fetch

# ---- composition conversions
M = [86.18, 92.14]
w = mole_to_mass([0.5, 0.5], M)
assert abs(w.sum() - 1) < 1e-12 and abs(w[0] - 86.18 / (86.18 + 92.14)) < 1e-12
assert np.allclose(mass_to_mole(w, M), [0.5, 0.5])
assert abs(mass_to_volume([0.5, 0.5], [660, 867]).sum() - 1) < 1e-12

# ---- hydrocarbon filter / classes
cases = {"CCCCCC": "n_paraffin", "CC(C)CC(C)(C)C": "iso_paraffin", "C1CCCCC1": "cycloalkane",
         "c1ccccc1": "aromatic", "C=CC": "alkene_alkyne"}
for s, c in cases.items():
    m = mol_from(smiles=s)
    assert is_hydrocarbon(m) and hc_class(m) == c, (s, hc_class(m))
assert not is_hydrocarbon(mol_from(smiles="CCO")) and not is_hydrocarbon(mol_from(smiles="[H][H]"))

# ---- Antoine: benzene (A=4.01814 B=1203.835 C=-53.226, 287.7-354.07 K) -> ~24.3 kPa at 40 C
p = wb.psat_kpa(4.01814, 1203.835, -53.226, 313.15)
assert 24.0 < p < 24.7, p
html = """<html><h1 id="Top">Benzene</h1><table>
<tr><th>Temperature (K)</th><th>A</th><th>B</th><th>C</th><th>Reference</th><th>Comment</th></tr>
<tr><td>287.7 to 354.07</td><td>4.01814</td><td>1203.835</td><td>-53.226</td><td>Willingham, 1945</td><td></td></tr>
<tr><td>422.11 to 535.9</td><td>4.19827</td><td>1094.254</td><td>-111.603</td><td>Ambrose, 1963</td><td></td></tr>
</table></html>"""
rows = wb.parse_antoine(html)
assert len(rows) == 2 and rows[0]["T_lo"] == 287.7 and rows[1]["C"] == -111.603
assert wb.page_name(html) == "Benzene"
b = max(rows, key=wb.coverage)
assert b["T_hi"] == 354.07 and wb.coverage(b) < 0.999   # does not reach 373 K -> flagged

# ---- WebBook page name must be the compound, not the site banner
banner = '<html><title>Hexane</title><body><h1>NIST Chemistry WebBook, SRD 69</h1><h1 id="Top">Hexane</h1></body></html>'
assert wb.page_name(banner) == "Hexane"
assert wb.page_name('<html><title>Decane</title><h1>NIST Chemistry WebBook, SRD 69</h1></html>') == "Decane"

# (ThermoML parsing is covered by tests/test_thermoml.py, which uses the schema-correct structure and a real record)

# ---- Zenodo file-list parsing
new = {"id": 1, "files": [{"key": "a.csv", "size": 5, "links": {"self": "https://zenodo.org/api/records/1/files/a.csv/content"}}]}
old = {"id": 2, "files": [{"filename": "b.csv", "size": 7, "links": {"download": "https://zenodo.org/record/2/files/b.csv"}}]}
assert fetch.zenodo_files(new)[0][:2] == ("a.csv", "https://zenodo.org/api/records/1/files/a.csv/content")
assert fetch.zenodo_files(old)[0][:2] == ("b.csv", "https://zenodo.org/record/2/files/b.csv")
print("ALL OFFLINE TESTS PASSED")
