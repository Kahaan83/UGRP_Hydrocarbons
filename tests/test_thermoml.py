"""ThermoML extractor regression tests (real record + synthetic cases). Run: python tests/test_thermoml.py"""
import collections, io, os, subprocess, sys, tarfile, tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import thermoml_extract as tx

def parse(root, max_p=1000.0, wanted=None):
    st = collections.Counter()
    comps, rows = tx.parse_file(root, "t.xml", wanted or tx.WANTED, max_p, st)
    return comps, rows, st

# ---- 1) REAL record (je8006138.xml, TRC): every point must have T, P, liquid phase and a complete composition
comps, rows, st = parse(ET.parse(HERE / "tests/data/je8006138.xml").getroot())
df = pd.DataFrame(rows)
assert len(df) == 150 and df.T_K.notna().all() and df.P_kPa.notna().all()
assert set(df.prop) == {"Mass density, kg/m3", "Viscosity, Pa*s"} and set(df.prop_phase) == {"Liquid"}
assert all(abs(sum(r["comp"].values()) - 1) < 1e-6 for r in rows) and set(df.comp_type) == {"pure", "mole fraction"}

# ---- synthetic builders (structure copied from the schema: Constraint > ConstraintID > ConstraintType ...)
NS = 'xmlns="http://www.iupac.org/namespaces/ThermoML"'
HEX = "InChI=1S/C6H14/c1-3-5-6-4-2/h3-6H2,1-2H3"; TOL = "InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3"
HEP = "InChI=1S/C7H16/c1-3-5-7-6-4-2/h3-7H2,1-2H3"; ETH = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
def compound(i, inchi): return f"<Compound><RegNum><nOrgNum>{i}</nOrgNum></RegNum><sStandardInChI>{inchi}</sStandardInChI><sCommonName>c{i}</sCommonName></Compound>"
def prop(n, name, phase="Liquid"):
    return (f"<Property><nPropNumber>{n}</nPropNumber><Property-MethodID><PropertyGroup><VolumetricProp><ePropName>{name}</ePropName>"
            f"</VolumetricProp></PropertyGroup></Property-MethodID><PropPhaseID><ePropPhase>{phase}</ePropPhase></PropPhaseID></Property>")
def var(n, tag, text, org=None, phase=None):
    reg = f"<RegNum><nOrgNum>{org}</nOrgNum></RegNum>" if org else ""
    ph = f"<VarPhaseID><ePhase>{phase}</ePhase></VarPhaseID>" if phase else ""
    return f"<Variable><nVarNumber>{n}</nVarNumber><VariableID><VariableType><{tag}>{text}</{tag}></VariableType>{reg}</VariableID>{ph}</Variable>"
def const(tag, text, value, org=None):
    reg = f"<RegNum><nOrgNum>{org}</nOrgNum></RegNum>" if org else ""
    return (f"<Constraint><nConstraintNumber>1</nConstraintNumber><ConstraintID><ConstraintType><{tag}>{text}</{tag}></ConstraintType>{reg}"
            f"</ConstraintID><nConstraintValue>{value}</nConstraintValue></Constraint>")
def nv(vars_, pvals):
    v = "".join(f"<VariableValue><nVarNumber>{n}</nVarNumber><nVarValue>{x}</nVarValue></VariableValue>" for n, x in vars_)
    p = "".join(f"<PropertyValue><nPropNumber>{n}</nPropNumber><nPropValue>{x}</nPropValue></PropertyValue>" for n, x in pvals)
    return f"<NumValues>{v}{p}</NumValues>"
def block(ids, props, body):
    cs = "".join(f"<Component><RegNum><nOrgNum>{i}</nOrgNum></RegNum></Component>" for i in ids)
    return f"<PureOrMixtureData>{cs}{props}{body}</PureOrMixtureData>"
def doc(*parts): return ET.parse(io.StringIO(f'<DataReport {NS}>' + "".join(parts) + "</DataReport>")).getroot()
CMP = compound(1, HEX) + compound(2, TOL) + compound(3, HEP) + compound(4, ETH)
DENS = "Mass density, kg/m3"

# ---- 2) FIXED temperature/pressure in Constraint blocks (the bug that dropped 61 % of temperatures)
r = doc(CMP, block([1], prop(1, DENS), const("eTemperature", "Temperature, K", 298.15) + const("ePressure", "Pressure, kPa", 101.3)
                   + var(1, "eMiscellaneous", "x") + nv([(1, 1)], [(1, 654.8)])))
_, rows, _ = parse(r); assert len(rows) == 1 and rows[0]["T_K"] == 298.15 and rows[0]["P_kPa"] == 101.3

# ---- 3) excess viscosity must be ignored (it produced NEGATIVE 'viscosities'); critical density ignored
r = doc(CMP, block([1], prop(1, "Excess viscosity, Pa*s") + prop(2, "Critical density, kg/m3") + prop(3, "Viscosity, Pa*s"),
                   const("eTemperature", "Temperature, K", 298.15) + nv([], [(1, -0.0001), (2, 232.0), (3, 0.0003)])))
_, rows, st = parse(r); assert [x["prop"] for x in rows] == ["Viscosity, Pa*s"] and st["other property (ignored)"] == 2

# ---- 4) VLE: vapour-phase composition must NOT overwrite the liquid composition
r = doc(CMP, block([1, 2], prop(1, "Vapor or sublimation pressure, kPa"),
                   var(1, "eTemperature", "Temperature, K") + var(2, "eComponentComposition", "Mole fraction", 1, "Liquid")
                   + var(3, "eComponentComposition", "Mole fraction", 1, "Gas") + nv([(1, 330.0), (2, 0.30), (3, 0.55)], [(1, 55.0)])))
_, rows, _ = parse(r); assert len(rows) == 1 and abs(rows[0]["comp"]["1"] - 0.30) < 1e-12 and abs(rows[0]["comp"]["2"] - 0.70) < 1e-12

# ---- 5) ternary with the last component implied; incomplete / molality compositions are dropped, not guessed
r = doc(CMP, block([1, 2, 3], prop(1, DENS), const("eTemperature", "Temperature, K", 298.15)
                   + var(1, "eComponentComposition", "Mole fraction", 1) + var(2, "eComponentComposition", "Mole fraction", 2)
                   + nv([(1, 0.2), (2, 0.3)], [(1, 700.0)]) + nv([(1, 0.2), (2, 0.9)], [(1, 700.0)])))
_, rows, st = parse(r); assert len(rows) == 1 and abs(rows[0]["comp"]["3"] - 0.5) < 1e-12 and st["mixture: inconsistent composition"] == 1
r = doc(CMP, block([1, 2], prop(1, DENS), const("eTemperature", "Temperature, K", 298.15)
                   + var(1, "eComponentComposition", "Molality, mol/kg", 1) + nv([(1, 0.5)], [(1, 700.0)])))
_, rows, st = parse(r); assert rows == [] and st["mixture: unsupported/ambiguous composition scale"] == 1

# ---- 6) high-pressure liquid, gas-phase density and solid-phase vapour pressure are all dropped
r = doc(CMP,
        block([1], prop(1, DENS), const("eTemperature", "Temperature, K", 298.15) + var(1, "ePressure", "Pressure, kPa") + nv([(1, 100000)], [(1, 700.0)]),),
        block([1], prop(1, DENS, "Gas"), const("eTemperature", "Temperature, K", 400.0) + nv([], [(1, 2.5)])),
        block([1], prop(1, "Vapor or sublimation pressure, kPa", "Crystal"), const("eTemperature", "Temperature, K", 150.0) + nv([], [(1, 0.01)])))
_, rows, st = parse(r); assert rows == []
assert st["density: pressure > 1000 kPa (dropped)"] == 1 and st["density: non-liquid phase (dropped)"] == 1 and st["vapour pressure of a solid (dropped)"] == 1

# ---- 7) end to end through a tarball: ethanol mixture dropped by the hydrocarbon filter, rows land in the CSV with the new columns
xml = (f'<DataReport {NS}>' + CMP + block([1, 2], prop(1, DENS), const("eTemperature", "Temperature, K", 298.15) + const("ePressure", "Pressure, kPa", 101.0)
       + var(1, "eComponentComposition", "Mass fraction", 1) + nv([(1, 0.5)], [(1, 760.0)]))
       + block([1, 4], prop(1, DENS), const("eTemperature", "Temperature, K", 298.15) + var(1, "eComponentComposition", "Mole fraction", 1)
               + nv([(1, 0.5)], [(1, 740.0)])) + "</DataReport>")
with tempfile.TemporaryDirectory() as d:
    d = Path(d); (d / "a.xml").write_text(xml)
    with tarfile.open(d / "t.tgz", "w:gz") as t: t.add(d / "a.xml", arcname="a.xml")
    out = subprocess.run([sys.executable, str(HERE / "thermoml_extract.py"), str(d / "t.tgz")], capture_output=True, text=True, cwd=d,
                         env={**os.environ, "PYTHONPATH": str(HERE)})
    assert out.returncode == 0 and "hydrocarbon rows kept: 1" in out.stdout and "not all-hydrocarbon" in out.stdout, out.stdout + out.stderr
    df = pd.read_csv(d / "data/processed/thermoml_hc.csv")
    assert list(df.columns[:3]) == ["file", "prop", "prop_phase"] and df.T_K.iloc[0] == 298.15 and df.comp_type.iloc[0] == "mass fraction"
print("ALL THERMOML TESTS PASSED")
