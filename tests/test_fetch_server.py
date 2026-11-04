"""Network-behaviour test for fetch.py using a local mock Zenodo. Run: python tests/test_fetch_server.py"""
import json, os, shutil, signal, subprocess, sys, tempfile, threading, time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
HERE = Path(__file__).resolve().parent.parent
PAYLOAD = os.urandom(3_000_000)
PDF = b"%PDF-1.4\n" + os.urandom(100_000)
hits = {"b": 0, "c": 0}

class H(BaseHTTPRequestHandler):
    def log_message(self, *a): pass
    def _file(self, data):
        rng = self.headers.get("Range")
        if rng:
            start = int(rng.split("=")[1].split("-")[0]); body = data[start:]
            self.send_response(206); self.send_header("Content-Length", str(len(body)))
        else:
            body = data; self.send_response(200); self.send_header("Content-Length", str(len(body)))
        self.end_headers(); self.wfile.write(body)
    def do_GET(self):
        p = self.path
        if p.startswith("/api/records/1"):
            base = f"http://127.0.0.1:{self.server.server_port}/f"
            files = [{"key": k, "size": len(PAYLOAD), "links": {"self": f"{base}/{k}"}} for k in ("a.csv", "b.csv")]
            files.append({"key": "c_stall.csv", "size": len(PAYLOAD), "links": {"self": f"{base}/c_stall.csv"}})
            body = json.dumps({"id": 1, "files": files}).encode()
            self.send_response(200); self.send_header("Content-Length", str(len(body))); self.end_headers(); self.wfile.write(body)
        elif p == "/f/a.csv": self._file(PAYLOAD)
        elif p == "/f/b.csv":
            hits["b"] += 1
            if hits["b"] == 1:
                self.send_response(429); self.send_header("Retry-After", "1"); self.send_header("Content-Length", "0"); self.end_headers()
            else: self._file(PAYLOAD)
        elif p == "/f/consent.pdf":
            body = b"<html>Please accept the terms of use</html>" * 30
            self.send_response(200); self.send_header("Content-Type", "text/html"); self.send_header("Content-Length", str(len(body))); self.end_headers(); self.wfile.write(body)
        elif p == "/f/dtic.pdf":
            if "Mozilla" not in self.headers.get("User-Agent", ""):
                self.send_response(403); self.send_header("Content-Length", "0"); self.end_headers()
            else: self._file(PDF)
        elif p == "/f/c_stall.csv":
            hits["c"] += 1
            self.send_response(200); self.send_header("Content-Length", str(len(PAYLOAD))); self.end_headers()
            self.wfile.write(PAYLOAD[:1000]); self.wfile.flush(); time.sleep(120)   # stalls mid-body
        else:
            self.send_response(404); self.send_header("Content-Length", "0"); self.end_headers()

srv = ThreadingHTTPServer(("127.0.0.1", 0), H); threading.Thread(target=srv.serve_forever, daemon=True).start()
api = f"http://127.0.0.1:{srv.server_port}/api/records"
work = Path(tempfile.mkdtemp())
env = {**os.environ, "ZENODO_API": api, "PYTHONPATH": str(HERE)}
code = ("import fetch; fetch.ZENODO={'t':'1'}; fetch.READ_TO=2; fetch.CONNECT_TO=2; "
        "fetch.fetch_zenodo(**__import__('json').loads(__import__('sys').argv[1])); fetch.summary()")
def run(kw): return subprocess.run([sys.executable, "-c", code, json.dumps(kw)], cwd=work, env=env,
                                   capture_output=True, text=True, timeout=120)
# make fetch importable from work dir
shutil.copy(HERE / "fetch.py", work / "fetch.py")

# 1) --list downloads nothing
r = run({"list_only": True}); assert "3 files" in r.stdout and not (work / "data/raw/zenodo").exists(), r.stdout + r.stderr
# 2) include filter + 429 handling (b.csv gets 429 once, retries after Retry-After)
r = run({"include": "^(a|b)\\.csv$"}); out = r.stdout
assert (work / "data/raw/zenodo/t/a.csv").read_bytes() == PAYLOAD and (work / "data/raw/zenodo/t/b.csv").read_bytes() == PAYLOAD, out
assert "HTTP 429, waiting 1s" in out, out
# 3) re-run skips existing files
r = run({"include": "^(a|b)\\.csv$"}); assert r.stdout.count("[skip]") == 2, r.stdout
# 4) stalled file times out instead of hanging, and is reported in the summary (bounded: 3 tries x ~2s + backoff)
t0 = time.time(); r = run({"include": "c_stall"}); dt = time.time() - t0
assert dt < 60 and "c_stall.csv" in r.stdout.split("FAILED")[1], (dt, r.stdout)
# 5) resume: pre-seed a partial file, server honours Range -> final file identical
part = work / "data/raw/zenodo/t/a2.csv.part"; part.write_bytes(PAYLOAD[:1_000_000])
sys.path.insert(0, str(work)); import importlib, fetch; fetch.READ_TO = 5; fetch.RAW = work / "data/raw"
fetch.MANIFEST = fetch.RAW / "MANIFEST.json"
ok = fetch.download(f"http://127.0.0.1:{srv.server_port}/f/a.csv", work / "data/raw/zenodo/t/a2.csv")
assert ok and (work / "data/raw/zenodo/t/a2.csv").read_bytes() == PAYLOAD
# 6) max-gb guard
r = run({"max_gb": 0.0001}); assert "exceeds --max-gb" in r.stdout, r.stdout
# 7) Ctrl-C mid-download: group is skipped, summary still printed, process exits cleanly
code2 = ("import fetch,_thread,threading;fetch.ZENODO={'t':'1'};fetch.READ_TO=3;fetch.CONNECT_TO=5;"
         "threading.Timer(1,_thread.interrupt_main).start();"  # os.kill(pid,SIGINT)=TerminateProcess on Windows;
         "fetch.main(['zenodo','--include','c_stall'])")       # and a blocked recv only checks signals when it returns
r = subprocess.run([sys.executable, "-c", code2], cwd=work, env=env, capture_output=True, text=True, timeout=60)
assert r.returncode == 0 and "interrupted by Ctrl-C" in r.stdout and "ALWAYS MANUAL" in r.stdout, r.stdout + r.stderr
# 8) 403 for scripts -> succeeds with browser-like headers; sha256 verification works both ways
import hashlib; good = hashlib.sha256(PAYLOAD).hexdigest()
u = f"http://127.0.0.1:{srv.server_port}/f/dtic.pdf"
assert fetch.download(u, work / "data/raw/pdfs/dtic.pdf")
assert fetch.download(f"http://127.0.0.1:{srv.server_port}/f/a.csv", work / "data/raw/x/ok.csv", expect_sha256=good)
before = len(fetch.failed)
assert not fetch.download(f"http://127.0.0.1:{srv.server_port}/f/a.csv", work / "data/raw/x/bad.csv", retries=1, expect_sha256="0" * 64)
assert not (work / "data/raw/x/bad.csv").exists() and len(fetch.failed) == before + 1
# 9) HTML page served as .pdf must be rejected (not saved, not recorded) and reported with a manual hint
before = len(fetch.failed)
assert not fetch.download(f"http://127.0.0.1:{srv.server_port}/f/consent.pdf", work / "data/raw/pdfs/consent.pdf")
assert not (work / "data/raw/pdfs/consent.pdf").exists() and len(fetch.failed) == before + 1 and "normal browser" in fetch.failed[-1]
# 10) an earlier-saved bad PDF is removed and refetched (valid content this time)
bad = work / "data/raw/pdfs/dtic.pdf"; bad.write_bytes(b"<html>x</html>")
assert fetch.download(u, bad) and fetch.valid_pdf(bad) and bad.stat().st_size == len(PDF)
print("ALL FETCH NETWORK TESTS PASSED")
