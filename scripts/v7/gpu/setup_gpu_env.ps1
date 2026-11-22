# Creates C:\fuelblend-gpu\.venv (Python 3.12) with CUDA 12.8 torch. Downloads about 3 GB: run only after approval.
$ErrorActionPreference = "Stop"
$root = "C:\fuelblend-gpu"
New-Item -ItemType Directory -Force "$root\runs" | Out-Null
py -3.12 -m venv "$root\.venv"
$py = "$root\.venv\Scripts\python.exe"
$repo = Split-Path -Parent (Split-Path -Parent (Split-Path -Parent $PSScriptRoot))
(Get-Content "$repo\requirements.txt") | Where-Object { $_ -notmatch '^torch' } | Set-Content "$root\req_no_torch.txt"
& $py -m pip install --upgrade pip
& $py -m pip install -r "$root\req_no_torch.txt"
& $py -m pip install torch --index-url https://download.pytorch.org/whl/cu128
& $py -m pip install optuna tensorboard
[Environment]::SetEnvironmentVariable("FUELBLEND_RUNS", "$root\runs", "User")
& $py "$repo\scripts\v7\gpu\check_gpu.py"
