"""Bounded closed-loop Boltz load test with per-request and GPU telemetry."""
import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import statistics
import subprocess
import threading
import time
import urllib.error
import urllib.request


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs', type=Path, default=Path('inputs/load_test_20260914_100'))
    p.add_argument('--concurrency', type=int, required=True)
    p.add_argument('--count', type=int, default=100, help='Total requests; cycles through the manifest if larger than the sample')
    p.add_argument('--timeout', type=float, default=180, help='Socket timeout in seconds; not a server-side cancellation deadline')
    p.add_argument('--endpoint', action='append', help='Repeat for independent workers; concurrency slots are assigned round-robin')
    p.add_argument('--label', default='sweep')
    a = p.parse_args()
    a.endpoint = a.endpoint or ['http://127.0.0.1:8001/biology/mit/boltz2/predict']
    if a.concurrency < 1 or a.count < 1 or not math.isfinite(a.timeout) or a.timeout <= 0:
        p.error('count, concurrency, and timeout must be positive and finite')
    if not a.label or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in a.label):
        p.error('label must contain only letters, numbers, underscores, and hyphens')
    manifest = json.loads((a.inputs / 'manifest.json').read_text())
    if not manifest.get('jobs'):
        p.error('Manifest must contain jobs')
    # Load before timing so missing files fail before any requests are submitted.
    bodies = [(a.inputs / job['request']).read_bytes() for job in manifest['jobs']]
    for body in bodies:
        json.loads(body)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out = Path('results/load_test') / f'{stamp}_{a.label}_c{a.concurrency}'
    out.mkdir(parents=True)
    config = vars(a).copy()
    config['inputs'] = str(a.inputs)
    (out / 'config.json').write_text(json.dumps(config, indent=2))
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    (out / 'benchmark_boltz.py').write_bytes(Path(__file__).read_bytes())
    stop = threading.Event()
    samples = []
    telemetry_errors = []
    started = time.monotonic()

    def monitor():
        with (out / 'gpu.jsonl').open('w') as f:
            while not stop.is_set():
                command = ['nvidia-smi', '--id=1', '--query-gpu=utilization.gpu,utilization.memory,memory.used,power.draw,temperature.gpu', '--format=csv,noheader,nounits']
                try:
                    result = subprocess.run(command, capture_output=True, text=True, timeout=5, check=True)
                    values = [float(x.strip()) for x in result.stdout.strip().split(',')]
                    if len(values) != 5:
                        raise ValueError('Unexpected nvidia-smi output')
                    sample = dict(zip(['gpu_util_pct', 'memory_util_pct', 'vram_mib', 'power_w', 'temperature_c'], values))
                    sample['elapsed_s'] = time.monotonic() - started
                    samples.append(sample)
                    f.write(json.dumps(sample) + '\n')
                    f.flush()
                except (ValueError, OSError, subprocess.SubprocessError) as error:
                    telemetry_errors.append(str(error))
                stop.wait(1)

    def request(i, endpoint):
        job = manifest['jobs'][i % len(manifest['jobs'])]
        body = bodies[i % len(bodies)]
        begin = time.monotonic()
        record = {'index': i, **job, 'endpoint': endpoint, 'start_s': begin - started}
        req = urllib.request.Request(endpoint, data=body, headers={'Content-Type': 'application/json'})
        try:
            with urllib.request.urlopen(req, timeout=a.timeout) as response:
                raw = response.read()
                record['status'] = response.status
            payload = json.loads(raw)
            if not payload.get('structures') or not payload.get('affinities'):
                raise ValueError('Response missing structures or affinities')
            record['success'] = True
            record['response_bytes'] = len(raw)
        except Exception as error:
            record['success'] = False
            record['error'] = str(error)
            record['transport_error'] = isinstance(error, OSError) and not isinstance(error, urllib.error.HTTPError)
            if isinstance(error, urllib.error.HTTPError):
                record['status'] = error.code
                (out / f'error_{i:04d}.txt').write_bytes(error.read())
        record['end_s'] = time.monotonic() - started
        record['latency_s'] = time.monotonic() - begin
        return record

    thread = threading.Thread(target=monitor, daemon=True)
    thread.start()
    records = []
    aborted = False
    print(f'Output: {out}', flush=True)
    try:
        with (out / 'requests.jsonl').open('w') as f, ThreadPoolExecutor(max_workers=a.concurrency) as pool:
            next_index = min(a.concurrency, a.count)
            pending = {pool.submit(request, i, a.endpoint[i % len(a.endpoint)]): a.endpoint[i % len(a.endpoint)]
                       for i in range(next_index)}
            while pending:
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                available_slots = []
                for future in done:
                    available_slots.append(pending.pop(future))
                    record = future.result()
                    records.append(record)
                    f.write(json.dumps(record) + '\n')
                    f.flush()
                    print(f"[{len(records)}/{a.count}] {record['request']} success={record['success']} latency={record['latency_s']:.2f}s", flush=True)
                    # A timed-out client does not cancel server compute. Stop feeding
                    # new jobs and drain existing clients before a separate recovery.
                    aborted = aborted or record.get('transport_error', False)
                for endpoint in available_slots:
                    if aborted or next_index >= a.count:
                        break
                    pending[pool.submit(request, next_index, endpoint)] = endpoint
                    next_index += 1
    finally:
        elapsed = time.monotonic() - started
        stop.set()
        thread.join()
    good = [r for r in records if r['success']]
    latencies = sorted(r['latency_s'] for r in good)
    summary = {'output': str(out), 'concurrency': a.concurrency, 'count': len(records),
               'requested_count': a.count, 'aborted_after_transport_error': aborted,
               'telemetry_samples': len(samples), 'telemetry_errors': telemetry_errors,
               'successes': len(good), 'failures': len(records) - len(good),
               'elapsed_s': elapsed, 'successful_requests_per_s': len(good) / elapsed,
               'latency_mean_s': statistics.mean(latencies) if latencies else None,
               'latency_p95_s': latencies[math.ceil(.95*len(latencies))-1] if latencies else None}
    if samples:
        summary.update(gpu_util_mean_pct=statistics.mean(s['gpu_util_pct'] for s in samples),
                       gpu_util_p95_pct=sorted(s['gpu_util_pct'] for s in samples)[int(.95*len(samples))],
                       vram_peak_mib=max(s['vram_mib'] for s in samples),
                       power_mean_w=statistics.mean(s['power_w'] for s in samples))
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)
    if aborted or len(good) != a.count:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
