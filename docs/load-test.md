# Boltz GPU 1 load test

## Multiple inference workers on GPU 1

The installed NIM creates one inference worker per visible GPU. Increasing HTTP
concurrency against one container queues work behind that worker. To test actual
overlapping inference, start an additional container using the same GPU:

```bash
docker compose -f compose.yaml -f compose.load-test.yaml up -d boltz-worker-2
curl --fail http://127.0.0.1:8002/v1/health/ready
python3 scripts/benchmark_boltz.py --concurrency 2 --count 100 --label two_workers \
  --endpoint http://127.0.0.1:8001/biology/mit/boltz2/predict \
  --endpoint http://127.0.0.1:8002/biology/mit/boltz2/predict
```

Wait for readiness before benchmarking. Each concurrency slot remains attached
to its endpoint; with two endpoints and concurrency 2, each worker has at most
one active client request. Whichever finishes first receives the next job.
Both services inherit GPU 1 and the same model cache from `boltz2`. Port 8002
is bound only to localhost. The override requires Compose support for `!override`.
Stop the additional worker with
`docker compose -f compose.yaml -f compose.load-test.yaml stop boltz-worker-2`.

The overlay also defines workers 3 and 4 on ports 8003 and 8004, both on GPU 1.
Start them by explicitly naming those services, then add their URLs with
`--endpoint` and use `--concurrency 4`. Stop them after the comparison to free
their model memory. Always use both Compose files when operating these services.

## Client and measurements

`scripts/benchmark_boltz.py` sends the prepared random BindingDB requests to
the existing Boltz API on port 8001. Compose assigns that server to GPU 1.
The script does not start containers, change GPU assignments, or tune inference
parameters. GPU telemetry is explicitly collected from GPU 1. If you override
the endpoint, you must verify that server's GPU assignment separately.

Run one stage at a time from the repository root:

```bash
python3 scripts/benchmark_boltz.py --concurrency 1 --count 100 --label baseline
python3 scripts/benchmark_boltz.py --concurrency 2 --count 100 --label sweep
```

Inspect each stage before increasing concurrency to 4, 8, 16, and beyond.
The client maintains at most `--concurrency` active HTTP requests and immediately
replaces each completed request until `--count` have been submitted. It creates
no additional inference workers. This measures the existing server's ability
to handle concurrent requests; queued HTTP requests may not compute concurrently.
The same manifest order is used at every level. Counts above 100 cycle through
the sample, allowing a longer run at a promising concurrency level:

```bash
python3 scripts/benchmark_boltz.py --concurrency 2 --count 500 --label sustained
```

Each timestamped directory under `results/load_test/` contains:

- `config.json` and `manifest.json`: run configuration and workload provenance.
- `requests.jsonl`: individual start/end times, latency, status, success, and input identity.
- `error_*.txt`: HTTP error response bodies. Failures are not retried.
- `gpu.jsonl`: approximately one-second GPU utilization, memory activity, VRAM,
  power, and temperature samples for GPU 1.
- `summary.json`: successful requests per second over the entire run, success
  and failure counts, successful-request mean/p95 latency, GPU utilization,
  peak VRAM, and mean power. Telemetry failures are explicitly recorded.

Successful responses must include nonempty structures and affinities. Full
prediction responses are discarded to avoid including output-storage costs.
All request latencies, including failures, remain available in `requests.jsonl`.
Summary latency percentiles cover successful requests only. End-to-end run
throughput includes initial fill, final drain, and any first-request warmup.
GPU utilization is the `nvidia-smi` activity metric, not SM occupancy or achieved
FLOPS. Memory utilization is memory activity, not a percentage of VRAM capacity.

`--timeout` is a socket timeout, not a total prediction deadline or server-side
cancellation. On a transport error the script stops submitting new work and
waits for existing client requests to finish. HTTP errors are recorded and the
stage continues. Any failed request produces exit status 1. After a timeout,
inspect the server before starting another stage: outstanding inference may
still be running even after the client has exited.

This is a closed-loop concurrency test, not a fixed-arrival-rate test. A rising
latency with flat throughput indicates queueing. A short successful stage does
not establish sustainability: repeat the best level for a longer run, inspect
throughput over time and error bodies, and distinguish invalid input failures
from overload. To establish a sustainable offered requests/second rate, follow
up with a paced arrival test below measured capacity and check for queue growth.
