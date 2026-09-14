"""Summarize completed Boltz load tests and plot sustained GPU 1 telemetry."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/load_test'


def main():
    runs = []
    for path in sorted(OUT.glob('*/summary.json')):
        summary = json.loads(path.read_text())
        config = json.loads((path.parent / 'config.json').read_text())
        endpoints = config.get('endpoint', [])
        workers = len(endpoints) if isinstance(endpoints, list) else 1
        runs.append((path.parent, summary, config, workers))
    lines = ['# GPU 1 Boltz load test', '',
             'NVIDIA H100 NVL; Boltz-2 NIM 1.9.0. All inference workers use GPU 1. '
             'GPU 0 was unused during these tests.', '',
             '| Run | Workers | Concurrent clients | Success / attempted | Successful req/s | Predictions/min | Mean latency (s) | Mean GPU activity | Peak VRAM (GiB) |',
             '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for path, s, config, workers in runs:
        lines.append(f"| [{config['label']}]({path.name}/summary.json) | {workers} | {s['concurrency']} | "
                     f"{s['successes']}/{s['count']} | {s['successful_requests_per_s']:.4f} | "
                     f"{60*s['successful_requests_per_s']:.2f} | {s['latency_mean_s'] or 0:.2f} | "
                     f"{s.get('gpu_util_mean_pct', 0):.1f}% | {s.get('vram_peak_mib', 0)/1024:.2f} |")
    lines += ['', '## Interpretation', '',
              'The one-worker concurrency sweep and two/four-worker comparisons use the same first '
              '12 randomly sampled rows. These short runs establish a preliminary comparison; '
              'a 1–2% difference is not convincing evidence of a throughput improvement.', '',
              'The overload burst uses the first 32 rows. The sustained test uses all 100 rows. '
              'Different protein/ligand sizes change inference cost, so those run-level rates are '
              'not controlled worker-count comparisons.', '',
              'The installed server creates one worker per visible GPU per container. Additional '
              'HTTP clients queue for that worker. Its GPU-acquisition timeout is 120 seconds. '
              'The overload error bodies distinguish that queue timeout from GPU out-of-memory.', '',
              'Throughput includes startup fill and final drain. Latency summaries cover successful '
              'requests only; individual failed request timings are in requests.jsonl. GPU activity '
              'comes from nvidia-smi and does not measure SM occupancy or achieved FLOPS.', '']
    sustained = [run for run in runs if run[2]['label'].startswith('sustained')]
    if sustained:
        path, s, _, _ = sustained[-1]
        records = [json.loads(line) for line in (path / 'requests.jsonl').read_text().splitlines()]
        samples = [json.loads(line) for line in (path / 'gpu.jsonl').read_text().splitlines()]
        lines += ['## Sustained run', '',
                  f"Completed {s['count']} requests in {s['elapsed_s']/60:.2f} minutes: "
                  f"{s['successes']} successes, {s['failures']} failures. "
                  'This measures this finite workload, not a long-term availability guarantee.', '']
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
        times = [sample['elapsed_s']/60 for sample in samples]
        axes[0].plot(times, [sample['gpu_util_pct'] for sample in samples], linewidth=.8)
        axes[0].set(ylabel='GPU activity (%)', ylim=(0, 105), title='Two Boltz workers on GPU 1 — 100 BindingDB rows')
        axes[1].plot(times, [sample['vram_mib']/1024 for sample in samples], color='tab:orange')
        axes[1].set(ylabel='VRAM (GiB)')
        completions = sorted(record['end_s'] for record in records if record['success'])
        axes[2].step([0]+[t/60 for t in completions], [0]+list(range(1,len(completions)+1)), where='post', color='tab:green')
        axes[2].set(ylabel='Successful predictions', xlabel='Elapsed minutes')
        for axis in axes:
            axis.grid(alpha=.25)
        fig.tight_layout()
        fig.savefig(OUT / 'sustained.png', dpi=160)
        plt.close(fig)
        lines += ['![Sustained workload telemetry](sustained.png)', '']
    (OUT / 'report.md').write_text('\n'.join(lines))
    print(OUT / 'report.md')


if __name__ == '__main__':
    main()
