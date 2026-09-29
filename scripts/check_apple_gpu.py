"""A bounded MPS installation check; no chemistry, weights or timing benchmark."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    target = args.output.resolve()
    if not target.is_relative_to((ROOT / 'reports/completion').resolve()) or target.exists():
        raise SystemExit('Use a new report under reports/completion; previous attempts are preserved')
    started = time.monotonic()
    result = dict(started_at=datetime.now(timezone.utc).isoformat(), status='checking',
        python=sys.version, executable=sys.executable,
        packages={p: importlib.metadata.version(p) for p in ('torch', 'mace-torch', 'numpy', 'ase')},
        fallback_environment=os.environ.get('PYTORCH_ENABLE_MPS_FALLBACK'),
        actual_model_inference_executed=False, dft_acceleration_verified=False,
        scope='Dimensionless tensor and gradient installation check only; no molecular model or speedup claim')
    try:
        import torch
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
        result.update(mps_built=torch.backends.mps.is_built(), mps_available=torch.backends.mps.is_available())
        if not result['mps_available']:
            raise RuntimeError('MPS is unavailable in this process environment')
        torch.manual_seed(17)
        x = torch.randn(128, 128, dtype=torch.float32)
        w = torch.randn(128, 128, dtype=torch.float32)
        cpu = x.double().requires_grad_()
        cpu_loss = (cpu @ w.double()).square().mean()
        cpu_loss.backward()
        gpu = x.to('mps').requires_grad_()
        gpu_loss = (gpu @ w.to('mps')).square().mean()
        gpu_loss.backward()
        torch.mps.synchronize()
        value = gpu_loss.detach().cpu().double()
        gradient = gpu.grad.detach().cpu().double()
        torch.testing.assert_close(value, cpu_loss.detach(), atol=2e-5, rtol=2e-5)
        torch.testing.assert_close(gradient, cpu.grad, atol=2e-5, rtol=2e-5)
        result.update(actual_result_device=str(gpu_loss.device),
            loss_absolute_error=float(abs(value - cpu_loss.detach())),
            maximum_gradient_absolute_error=float((gradient - cpu.grad).abs().max()),
            installation_check_atol=2e-5, installation_check_rtol=2e-5,
            allocated_gpu_bytes=torch.mps.current_allocated_memory())
        try:
            torch.ones(1, dtype=torch.float64, device='mps')
            result['float64_allocation_supported'] = True
        except (TypeError, RuntimeError) as exc:
            result.update(float64_allocation_supported=False, float64_allocation_error=str(exc))
        result['status'] = 'mps_tensor_and_gradient_check_passed'
    except Exception as exc:
        result.update(status='unavailable_or_check_failed', error=f'{type(exc).__name__}: {exc}')
    result['elapsed_seconds_including_import'] = time.monotonic() - started
    result['source_sha256'] = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
        for p in ('scripts/check_apple_gpu.py', 'environment/mps-requirements.in', 'environment/mps-requirements.lock')}
    with target.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))
    return 0 if result['status'] == 'mps_tensor_and_gradient_check_passed' else 2


if __name__ == '__main__':
    raise SystemExit(main())
