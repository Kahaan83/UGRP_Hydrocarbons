"""Patched ajftd_scraper.py: name fix + resume helper + unchanged extraction (no browser needed). Run: python tests/test_scraper_patch.py"""
import os, sys, tempfile
from pathlib import Path
HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
_cwd0 = os.getcwd(); os.chdir(tempfile.mkdtemp())   # save_html writes into CWD; keep it out of the repo
import ajftd_scraper as sc

try:
    HDR = ["Value", "Unit", "temperature_value", "temperature_unit", "specification_limits", "info", "Test Date", "Test Method"]
    def table(tid, rows):
        return f'<table id="{tid}"><tr>' + "".join(f"<th>{c}</th>" for c in HDR) + "</tr>" + "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows) + "</table>"
    html = ("<html><head><title>View Fuel | AJF:TD | U of I</title></head><body><h1>Alternative Jet Fuels</h1><h2>View Fuel</h2>"
            "<h3>Select Metrics to View</h3><h3>Fuel Details for (02239-070416)T108</h3><h3>Property Metrics:</h3>"
            "<button onclick=\"myFunction('tbl_24')\">flash_point</button>" + table("tbl_24", [["43.5", "C", "", "", "38 Min", "", "4/17/2015 00:00", "D56"]]) + "</body></html>")
    rec = sc.extract_fuel_record(html, "https://altjetfuels.illinois.edu/fuels/fuel-viewcompare.asp?id=429&action=view_fuel")
    assert rec["name"] == "(02239-070416)T108" and rec["fuel_label"] == "(02239-070416)T108" and rec["fuel_id"] == "429"
    fp = rec["properties"]["flash_point"][0]["measurements"][0]
    assert fp["Value"] == "43.5" and fp["Test Method"] == "D56"                      # extraction itself unchanged
    # a page with no 'Fuel Details for' heading falls back to the old behaviour instead of crashing
    rec2 = sc.extract_fuel_record(html.replace("Fuel Details for (02239-070416)T108", "Something else"), "https://x/?id=1")
    assert rec2["name"] == "Alternative Jet Fuels" and rec2["fuel_label"] is None

    # resume helper: only complete, previously saved pages are reused
    assert sc.cached_detail_html("429") is None
    sc.save_html(html + " " * 3000, "detail_429.html"); assert sc.cached_detail_html("429") is not None
    sc.save_html("<html>Login required</html>" + " " * 3000, "detail_430.html"); assert sc.cached_detail_html("430") is None      # not a fuel page
    sc.save_html("tiny", "detail_431.html"); assert sc.cached_detail_html("431") is None
    sc.RESUME = False; assert sc.cached_detail_html("429") is None                                                              # --fresh
finally:
    os.chdir(_cwd0)
print("ALL SCRAPER PATCH TESTS PASSED")
