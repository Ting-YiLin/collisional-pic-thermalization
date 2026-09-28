"""Compile the recovered C++ source and execute a deliberately tiny end-to-end run."""
from pathlib import Path
import shutil, subprocess, tempfile, sys, math

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "pic_2010_recovered.cpp"
CFG = ROOT / "config" / "smoke_input.txt"

compiler = shutil.which("g++")
if compiler is None:
    raise SystemExit("SMOKE_SKIP: g++ not found")

with tempfile.TemporaryDirectory(prefix="pic_smoke_") as td:
    td = Path(td)
    exe = td / ("pic_thermalization.exe" if sys.platform.startswith("win") else "pic_thermalization")
    subprocess.run([compiler, "-std=c++17", "-O2", str(SRC), "-o", str(exe)], check=True)
    shutil.copy2(CFG, td / "input.txt")
    cp = subprocess.run([str(exe)], cwd=td, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=60)
    if cp.returncode != 0:
        print(cp.stdout)
        print(cp.stderr, file=sys.stderr)
        raise SystemExit(f"SMOKE_FAIL: exit={cp.returncode}")
    required = ["space.txt", "time.txt", "Hermite_polynomials.txt", "output.txt", "fitting_number(omega_p*tau_R).txt"]
    missing = [name for name in required if not (td/name).exists() or (td/name).stat().st_size == 0]
    if missing:
        raise SystemExit("SMOKE_FAIL: missing outputs: " + ", ".join(missing))
    print("SMOKE_PASS: compile + reduced end-to-end run + required outputs")
