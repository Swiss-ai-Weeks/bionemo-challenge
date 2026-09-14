"""Download example data, prepare inputs, wait for NIMs, and run predictions."""
import argparse
import concurrent.futures
from pathlib import Path
import subprocess
import sys
import time
import urllib.request
root = Path(__file__).resolve().parents[1]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--with-openfold', action='store_true', help='Also run the optional OpenFold service')
args = parser.parse_args()
models = [('boltz2', 8001)]
if args.with_openfold:
    models.append(('openfold3', 8000))

subprocess.run([sys.executable, str(root / 'scripts/download_bindingdb.py')], check=True)
subprocess.run([sys.executable, str(root / 'scripts/prepare_example.py')], check=True)

def run(model, port):
    deadline = time.monotonic() + 3600
    while True:
        try:
            with urllib.request.urlopen(f'http://127.0.0.1:{port}/v1/health/ready', timeout=5) as r:
                if r.status == 200:
                    break
        except OSError:
            pass
        if time.monotonic() >= deadline:
            raise TimeoutError(f'{model} was not ready within one hour')
        time.sleep(15)
    print(f'{model} ready', flush=True)
    subprocess.run([sys.executable, str(root / 'scripts/predict.py'), model], check=True)

with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    futures = [pool.submit(run, name, port) for name, port in models]
    errors = []
    for f in concurrent.futures.as_completed(futures):
        try:
            f.result()
        except Exception as e:
            errors.append(str(e))
            print(e, flush=True)
subprocess.run([sys.executable, str(root / 'scripts/make_viewer.py')], check=True)
if errors:
    raise SystemExit('; '.join(errors))
