"""Run bounded real engine checks and keep raw inputs/outputs plus provenance."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import resource
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("engine", choices=["pyscf", "xtb", "qe", "orca"])
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("timeout must be positive")
    name = args.engine
    executable = {"pyscf": sys.executable, "xtb": str(ROOT / ".tools/xtb/bin/xtb"),
                  "qe": str(ROOT / ".tools/bin/pw.x"),
                  "orca": shutil.which("orca")}.get(name)
    if not executable or not Path(executable).is_file():
        raise SystemExit(f"{name} executable unavailable; see docs/ORCA_SETUP.ko.md for ORCA.")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    run = ROOT / "runs" / f"{stamp}-{name}"
    run.mkdir(parents=True)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=env.get("DFT_THREADS", "4"), OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    meta = {"purpose": "installation smoke test, not a scientific benchmark",
            "engine": name, "charge": 0, "spin_polarized": False,
            "geometry_optimization": "not_run", "minimum_verification": "not_run",
            "thread_environment": {k: env[k] for k in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]}}
    sources = []
    if name == "pyscf":
        sources = [ROOT / "scripts/pyscf_water.py"]
        cmd = [executable, "pyscf_water.py"]
        meta.update(version=importlib.metadata.version("pyscf"), method="PBE/def2-SVP", system="H2O", solvent="vacuum", multiplicity=1)
    elif name == "xtb":
        sources = [ROOT / "structures/molecules/water.xyz"]
        cmd = [executable, "water.xyz", "--gfn", "2", "--chrg", "0", "--uhf", "0"]
        meta.update(version="6.7.1", method="GFN2-xTB (semiempirical, not DFT)", system="H2O", multiplicity=1, solvent="vacuum")
    elif name == "orca":
        sources = [ROOT / "inputs/orca/water_sp.inp"]
        cmd = [str(Path(executable).resolve()), "water_sp.inp"]
        meta.update(method="PBE/def2-SVP", system="H2O", multiplicity=1, solvent="vacuum", execution="serial")
    else:
        sources = [ROOT / "inputs/qe/silicon_scf.in"]
        cmd = [executable, "-in", "silicon_scf.in"]
        pseudo = ROOT / "inputs/qe/pseudo/Si_ONCV_PBE-1.1.upf"
        (run / "pseudo").mkdir()
        shutil.copy2(pseudo, run / "pseudo" / pseudo.name)
        (run / "scratch").mkdir()
        meta.update(version="7.6", method="PBE/ONCV", system="Si2 diamond", ecutwfc_Ry=30,
                    ecutrho_Ry=120, k_grid=[2, 2, 2], k_offset=[0, 0, 0],
                    lattice_parameter_angstrom=5.43, SOC=False,
                    pseudo_sha256=sha(pseudo), accuracy="unconverged demonstration settings")
    for src in sources:
        shutil.copy2(src, run / src.name)
    meta.update(input_sha256={p.name: sha(p) for p in sources},
                executable=str(Path(executable).resolve()), executable_sha256=sha(Path(executable)),
                command=cmd, cwd=str(run), timeout_seconds=args.timeout,
                started_at=datetime.now(timezone.utc).isoformat())
    (run / "manifest.json").write_text(json.dumps(meta, indent=2) + "\n")
    start = time.monotonic()
    timed_out = False
    with (run / "output.log").open("w") as out, (run / "stderr.log").open("w") as err:
        proc = subprocess.Popen(cmd, cwd=run, env=env,
                                stdout=out, stderr=err, start_new_session=True)
        try:
            code = proc.wait(timeout=args.timeout)
        except (subprocess.TimeoutExpired, KeyboardInterrupt):
            timed_out = True
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            code = proc.returncode
    text = (run / "output.log").read_text(errors="replace")
    errtext = (run / "stderr.log").read_text(errors="replace")
    if name == "pyscf":
        result = json.loads((run / "result.json").read_text()) if (run / "result.json").exists() else {}
        converged = result.get("scf_converged", False)
        energy, unit = result.get("energy_hartree"), "Hartree"
    else:
        if name == "qe":
            converged = "convergence has been achieved" in text and "JOB DONE." in text
            pattern, unit = r"!\s+total energy\s+=\s+([-\d.]+)\s+Ry", "Ry"
            version = re.search(r"Program PWSCF v\.([^\s]+)", text)
            meta["version"] = version.group(1) if version else None
        elif name == "xtb":
            converged = "convergence criteria satisfied" in text and "normal termination of xtb" in text + errtext
            pattern, unit = r"TOTAL ENERGY\s+([-\d.]+)\s+Eh", "Hartree"
            version = re.search(r"xtb version\s+([^\s]+)", text)
            meta["version"] = version.group(1) if version else None
        else:
            converged = "ORCA TERMINATED NORMALLY" in text and "SCF CONVERGED" in text
            pattern, unit = r"FINAL SINGLE POINT ENERGY\s+([-\d.]+)", "Hartree"
            version = re.search(r"Program Version\s+(\S+)", text)
            meta["version"] = version.group(1) if version else None
        matches = re.findall(pattern, text)
        energy = float(matches[-1]) if matches else None
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    rss_bytes = int(rss if sys.platform == "darwin" else rss * 1024)
    status = "passed" if code == 0 and converged and energy is not None and not timed_out else "failed"
    meta.update(status=status, ended_at=datetime.now(timezone.utc).isoformat(),
                elapsed_seconds=round(time.monotonic() - start, 3), exit_code=code,
                timed_out=timed_out, scf_converged=bool(converged), energy=energy, energy_unit=unit,
                max_rss_bytes=rss_bytes,
                memory_measurement="getrusage(RUSAGE_CHILDREN).ru_maxrss; bytes on macOS; maximum child RSS, not summed process-tree or system unified-memory peak",
                output_sha256=sha(run / "output.log"))
    (run / "manifest.json").write_text(json.dumps(meta, indent=2) + "\n")
    report = ROOT / "reports" / f"smoke-{name}.json"
    report.write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps({k: meta[k] for k in ["engine", "status", "elapsed_seconds", "max_rss_bytes", "energy", "energy_unit"]}, indent=2))
    print(f"Raw run: {run}")
    raise SystemExit(0 if status == "passed" else 1)


if __name__ == "__main__":
    main()
