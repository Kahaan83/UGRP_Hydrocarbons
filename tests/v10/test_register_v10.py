import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "v10"))
from register_pdfs import register  # noqa: E402


def test_register_rejects_html(tmp_path):
    (tmp_path / "a.pdf").write_bytes(b"<!DOCTYPE html><html></html>")
    df = register(tmp_path, expected=[])
    r = df.iloc[0]
    assert not r.is_pdf and r.status == "not_pdf"


def test_register_hash(tmp_path):
    b = b"<!DOCTYPE x>"
    (tmp_path / "a.pdf").write_bytes(b)
    df = register(tmp_path, expected=[])
    assert df.iloc[0].sha256 == hashlib.sha256(b).hexdigest()


def test_register_not_provided(tmp_path):
    df = register(tmp_path, expected=["AIAA 2017-0146"])
    assert df.iloc[0].status == "not provided"
