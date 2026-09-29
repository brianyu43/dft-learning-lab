"""Read-only installation inventory. ORCA may legitimately be pending."""
from datetime import datetime, timezone
import importlib
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
packages = {
    "ase": "ase", "rdkit": "rdkit.Chem", "pymatgen": "pymatgen.core",
    "spglib": "spglib", "pyscf": "pyscf", "geometric": "geometric",
    "numpy": "numpy", "scipy": "scipy", "pandas": "pandas",
    "matplotlib": "matplotlib", "jupyterlab": "jupyterlab",
    "ipykernel": "ipykernel", "cclib": "cclib", "py3dmol": "py3Dmol",
}
result = {"checked_at": datetime.now(timezone.utc).isoformat(),
          "python": sys.version, "architecture": platform.machine(),
          "macOS": platform.mac_ver()[0], "packages": {}, "engines": {}}
for package, module in packages.items():
    try:
        importlib.import_module(module)
        result["packages"][package] = {"version": importlib.metadata.version(package), "import": "ok"}
    except Exception as exc:
        result["packages"][package] = {"import": "failed", "error": str(exc)}
for name in ["pw.x", "ph.x", "bands.x", "dos.x", "projwfc.x", "pp.x", "q2r.x", "matdyn.x", "dynmat.x", "xtb", "orca"]:
    path = shutil.which(name)
    result["engines"][name] = {"path": path, "status": "available" if path else "pending"}
    if path:
        result["engines"][name]["architecture"] = subprocess.check_output(["file", "-b", str(Path(path).resolve())], text=True).strip()
for name in ["pyscf", "qe", "xtb", "orca"]:
    report = ROOT / f"reports/smoke-{name}.json"
    if report.exists():
        result.setdefault("latest_smoke_tests", {})[name] = json.loads(report.read_text())["status"]
(ROOT / "reports/doctor.json").write_text(json.dumps(result, indent=2) + "\n")
for name, value in result["packages"].items():
    print(f"{name:14} {value.get('version', '-'):<14} {value['import']}")
for name, value in result["engines"].items():
    print(f"{name:14} {value['status']}")
print("Saved reports/doctor.json")
if any(x["import"] != "ok" for x in result["packages"].values()):
    raise SystemExit(1)
if any(result["engines"][x]["status"] != "available" for x in ["pw.x", "xtb"]):
    raise SystemExit(1)
