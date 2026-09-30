"""A local archive of full-run folders the repository does not hold, with its checksum, checked member by
member against the files it was made from (EVALUATION_GUIDE.md, section 8).

    python -m full_run_28092026.backup_archive scores                          # full_run_scores_<date>.zip
    python -m full_run_28092026.backup_archive traces/repeat1 traces/repeat2 traces/repeat3 --name full_run_traces_repeats

The archive goes to ~/EngTrace_private_backup/, with a `.sha256` in sha256sum's format like the earlier
backups; members are named from full_run_28092026/. Python's zipfile reads a file another process holds
open, which PowerShell's Compress-Archive refuses (it failed so on the repeat traces, D-151). A file that
changes while it is archived, such as a reply store under a running stage, is reported as differing:
archive a stage's store after the stage ends.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEST = Path.home() / 'EngTrace_private_backup'


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('folders', nargs='+', help='folders under full_run_28092026/')
    ap.add_argument('--name', help='the archive name before the date (default: full_run_<last part of the first folder>)')
    a = ap.parse_args()
    files = sorted(f for d in a.folders for f in (HERE / d).rglob('*') if f.is_file())
    if not files:
        raise SystemExit(f'nothing to archive under {a.folders}')
    name = f"{a.name or 'full_run_' + Path(a.folders[0]).name}_{datetime.date.today():%Y-%m-%d}.zip"
    DEST.mkdir(exist_ok=True)
    dest = DEST / name
    with zipfile.ZipFile(dest, 'w', zipfile.ZIP_DEFLATED) as z:
        for f in files:
            z.write(f, f.relative_to(HERE).as_posix())
    digest = hashlib.sha256(dest.read_bytes()).hexdigest()
    (DEST / f'{name}.sha256').write_text(f'{digest}  {name}\n', encoding='ascii', newline='\n')
    sha = lambda b: hashlib.sha256(b).hexdigest()   # noqa: E731
    with zipfile.ZipFile(dest) as z:
        differ = [f for f in files if sha(z.read(f.relative_to(HERE).as_posix())) != sha(f.read_bytes())]
    for f in differ:
        print(f'  DIFFERS: {f.relative_to(HERE).as_posix()}')
    print(f'{len(files)} files, {len(differ)} differ from their source; {dest} {dest.stat().st_size:,} bytes; '
          f'sha256 {digest}')
    return 1 if differ else 0


if __name__ == '__main__':
    raise SystemExit(main())
