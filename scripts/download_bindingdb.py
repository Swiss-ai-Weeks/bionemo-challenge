"""Cache the pinned BindingDB Articles release and extract the example sample."""
import csv
import io
import itertools
from pathlib import Path
import shutil
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = 'BindingDB_BindingDB_Articles_202609_tsv.zip'
URL = f'https://www.bindingdb.org/rwd/bind/downloads/{ARCHIVE}'
MEMBER = 'BindingDB_BindingDB_Articles.tsv'


def prepare_sample(data_dir):
    archive = data_dir / 'raw' / ARCHIVE
    archive.parent.mkdir(parents=True, exist_ok=True)
    if not archive.exists():
        temporary = archive.with_suffix('.zip.part')
        print(f'Downloading {URL}', flush=True)
        try:
            with urllib.request.urlopen(URL, timeout=120) as response, temporary.open('wb') as output:
                shutil.copyfileobj(response, output)
            with zipfile.ZipFile(temporary) as source:
                source.getinfo(MEMBER)
                if source.testzip() is not None:
                    raise ValueError('BindingDB archive failed its CRC check')
            temporary.replace(archive)
        finally:
            temporary.unlink(missing_ok=True)
    else:
        print(f'Using cached {archive}', flush=True)

    # Read just the sample; avoid unpacking the entire dataset to disk.
    with zipfile.ZipFile(archive) as source, source.open(MEMBER) as member:
        reader = csv.reader(io.TextIOWrapper(member, encoding='utf-8'), delimiter='\t')
        rows = list(itertools.islice(reader, 6))
    if len(rows) != 6:
        raise ValueError('Expected a header and five BindingDB records')
    sample = data_dir / 'BindingDB_top5.tsv'
    temporary = sample.with_suffix('.tsv.part')
    try:
        with temporary.open('w', encoding='utf-8', newline='') as output:
            csv.writer(output, delimiter='\t', lineterminator='\n').writerows(rows)
        temporary.replace(sample)
    finally:
        temporary.unlink(missing_ok=True)
    print(f'Prepared {sample}', flush=True)


if __name__ == '__main__':
    prepare_sample(ROOT / 'data')
