"""Verify current working-tree data files against archive/baseline_manifest.json checksums."""
import hashlib
import json
import os

m = json.load(open("archive/baseline_manifest.json"))
ok = bad = miss = 0
for f in m["files"]:
    p = f["path"]
    if not os.path.exists(p):
        print("MISSING", p)
        miss += 1
        continue
    h = hashlib.sha256(open(p, "rb").read()).hexdigest()
    if h == f["sha256"]:
        ok += 1
    else:
        print("CHANGED", p, os.path.getsize(p), "vs manifest", f["bytes"])
        bad += 1
print(f"matched={ok} changed={bad} missing={miss} of {len(m['files'])} manifest files")
print("manifest paths:")
for f in m["files"]:
    print("  ", f["path"])
