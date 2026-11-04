"""Offline AJFTD HTML parser + ingest, using pages built like the real ones your scraper parses. Run: python tests/test_ajftd_html.py"""
import shutil, subprocess, sys, tempfile
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent.parent
W = Path(tempfile.mkdtemp()); [shutil.copy(f, W / f.name) for f in HERE.glob("*.py")]
(W / "debug").mkdir()
HDR = ["Value", "Unit", "temperature_value", "temperature_unit", "specification_limits", "info", "Test Date", "Test Method"]

def table(tid, rows):
    h = "".join(f"<th>{c}</th>" for c in HDR)
    b = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<table id="{tid}"><tr>{h}</tr>{b}</table>'
def btn(tid, label): return f"<button onclick=\"myFunction('{tid}')\">{label}</button>"
def page(label, props):
    comp = btn("tbl_1", "aromatics_total") + table("tbl_1", [["18.2", "% vol", "", "", "25 Max", "", "7/3/2016 00:00", "D1319"]])
    body = "".join(btn(t, m) + table(t, r) for t, m, r in props)
    return ("<html><head><title>View Fuel | AJF:TD | U of I</title></head><body><h1>Alternative Jet Fuels</h1><h2>View Fuel</h2>"
            "<h3>Select Metrics to View</h3><h3>Fuel Details for " + label + "</h3><table><tr><th>Producer</th><th>Feedstock</th></tr>"
            "<tr><td>ACME</td><td>Camelina</td></tr></table><h3>Composition Metrics:</h3>" + comp + "<h3>Property Metrics:</h3>" + body + "</body></html>")
(W / "debug/detail_365.html").write_text(page("(02239-070416)T108", [
    ("tbl_15", "density", [["No Test", "kg/m3", "No Test", "C", "775-840", "", "7/3/2016 00:00", "D4052"]]),
    ("tbl_16", "density_api_by_meter", [["41.8", "API", "", "", "37-51", "", "7/3/2016 00:00", "D1298"]]),
    ("tbl_40", "freezing_point", [["-52.5", "C", "", "", "-40 Max", "", "7/3/2016 00:00", "D5972"]])]), encoding="utf-8")
(W / "debug/detail_429.html").write_text(page("(JP-5 batch 7)", [
    ("tbl_24", "flash_point", [["43.5", "C", "", "", "38 Min", "", "4/17/2015 00:00", "D56"]]),
    ("tbl_33", "viscosity_kinematic", [["4.885", "mm2/s", "-20", "C", "8 Max", "", "4/17/2015 00:00", "D445"]])]), encoding="utf-8")
(W / "debug/detail_999.html").write_text("<html><body><h1>Alternative Jet Fuels</h1><p>no tables</p></body></html>", encoding="utf-8")

# a scraper CSV for the cross-check (what ajftd_scraper.py would have written for fuel 429)
pd.DataFrame([dict(fuel_id="429", metric="flash_point", measurement_index=1, Value="43.5"),
              dict(fuel_id="429", metric="viscosity_kinematic", measurement_index=1, Value="4.885")]).to_csv(W / "ajftd_data.csv", index=False)

r = subprocess.run([sys.executable, "ajftd_parse_html.py", "--compare", "ajftd_data.csv"], cwd=W, capture_output=True, text=True)
assert r.returncode == 0, r.stdout + r.stderr
out = r.stdout
assert "pages parsed: 3" in out and "pages with no metric tables: 1" in out and "0 differ" in out, out
a = pd.read_csv(W / "data/raw/ajftd/ajftd_all_metrics.csv", dtype=str)
assert set(a.fuel_label) == {"(02239-070416)T108", "(JP-5 batch 7)"}                         # real names, not the site banner
assert "freezing_point" in set(a.metric)                                                    # a metric the old scraper skipped
inv = pd.read_csv(W / "data/raw/ajftd/ajftd_metrics_inventory.csv")
assert inv.set_index("metric").loc["density", "fuels_with_number"] == 0 and inv.set_index("metric").loc["flash_point", "fuels_with_number"] == 1
det = pd.read_csv(W / "data/raw/ajftd/ajftd_fuel_details.csv"); assert {"ACME", "Camelina"} <= set(det.value)   # description table captured

# ...and it flows through the ingester
(W / "data/raw/ajftd/ajftd_data.csv").unlink(missing_ok=True)
r = subprocess.run([sys.executable, "ingest_ajftd.py"], cwd=W, capture_output=True, text=True); assert r.returncode == 0, r.stdout + r.stderr
f = pd.read_csv(W / "data/processed/ajftd_fuels.csv"); g = lambda p: f[f.property == p]
assert abs(g("freezing_point_C").value.iloc[0] + 52.5) < 1e-9 and abs(g("flash_point_C").value.iloc[0] - 43.5) < 1e-9
v = g("viscosity_kin_cSt").iloc[0]; assert abs(v.value - 4.885) < 1e-9 and abs(v.T_K - 253.15) < 1e-6          # -20 C parsed
assert abs(g("density_kg_m3").value.iloc[0] - 141.5 / (41.8 + 131.5) * 999.016) < 1e-6                       # API gravity -> kg/m3
assert set(f.fuel_label) == {"(02239-070416)T108", "(JP-5 batch 7)"}
print("ALL AJFTD HTML TESTS PASSED")
