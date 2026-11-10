"""The Task 17 verifier must stay independent of the code it audits."""
import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN = ("fuelblend",)


def test_verifier_imports_no_fuelblend_code():
    tree = ast.parse((ROOT / "scripts/v4/verify_designs_v4.py").read_text(encoding="utf-8"))
    mods = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            mods += [a.name for a in n.names]
        elif isinstance(n, ast.ImportFrom):
            mods.append(n.module or "")
    bad = [m for m in mods if m.split(".")[0] in FORBIDDEN]
    assert not bad, f"verifier imports audited code: {bad}"
    assert any(m.split(".")[0] in ("thermo", "chemicals") for m in mods)
    src = (ROOT / "scripts/v4/verify_designs_v4.py").read_text(encoding="utf-8")
    assert "sys.path.insert" not in src            # no path hack that would make fuelblend importable
