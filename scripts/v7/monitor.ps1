# Live curves: opens TensorBoard on the run directory (http://localhost:6006).
$runs = if ($env:FUELBLEND_RUNS) { $env:FUELBLEND_RUNS } else { "C:\fuelblend-gpu\runs" }
New-Item -ItemType Directory -Force $runs | Out-Null
$tb = if (Test-Path "C:\fuelblend-gpu\.venv\Scripts\tensorboard.exe") { "C:\fuelblend-gpu\.venv\Scripts\tensorboard.exe" } else { "tensorboard" }
& $tb --logdir $runs --port 6006
