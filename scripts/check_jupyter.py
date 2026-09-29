"""Start a temporary localhost Jupyter server, verify it, then stop it."""
from pathlib import Path
import json
import os
import secrets
import socket
import subprocess
import sys
import time
from urllib.request import Request, urlopen

root = Path(__file__).resolve().parents[1]
token = secrets.token_urlsafe(32)
env = dict(os.environ, JUPYTER_TOKEN=token)
with socket.socket() as s:
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
log = root / "reports/jupyter-startup.log"
result = {"status": "failed", "host": "127.0.0.1", "persistent_server": False}
with log.open("w") as output:
    proc = subprocess.Popen([sys.executable, "-m", "jupyterlab", "--no-browser",
                             "--ServerApp.ip=127.0.0.1", f"--ServerApp.port={port}",
                             "--ServerApp.port_retries=0", "--ServerApp.log_level=ERROR"],
                            cwd=root, env=env, stdout=output, stderr=subprocess.STDOUT)
    try:
        deadline = time.monotonic() + 45
        while time.monotonic() < deadline:
            if proc.poll() is not None:
                raise RuntimeError("Jupyter exited before responding")
            try:
                req = Request(f"http://127.0.0.1:{port}/api/status", headers={"Authorization": "token " + token})
                with urlopen(req, timeout=1) as response:
                    status = json.loads(response.read())
                req = Request(f"http://127.0.0.1:{port}/lab", headers={"Authorization": "token " + token})
                with urlopen(req, timeout=3) as response:
                    html = response.read().decode()
                    assert response.status == 200 and "jupyter" in html.lower()
                result.update(status="passed", api_status="HTTP 200", lab_page="HTTP 200", kernels=status["kernels"])
                break
            except (OSError, ValueError):
                time.sleep(0.5)
        else:
            raise RuntimeError("Jupyter did not become ready within 45 seconds")
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
        result["server_stopped"] = proc.poll() is not None
log.write_text(log.read_text().replace(token, "[REDACTED]"))
(root / "reports/jupyter-check.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
