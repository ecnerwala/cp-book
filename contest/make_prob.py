#!/usr/bin/env python3
"""Create problem directories from layered .template directories.

Searches upward from the current directory for `.template` directories,
continuing past any that contain a `__PARENT__` marker file. Templates are
applied top-down, so deeper `.template` directories override shallower ones
file-by-file.

`__PROBLEM_NAME__` is substituted in file names and file contents; nothing
else is touched (in particular `$` has no special meaning), so Makefiles and
shell scripts template safely. If the resulting directory contains an
executable `setup` file, it is run from inside the new directory.

Usage:
  make_prob.py <name>...
"""

import os
import re
import shutil
import stat
import subprocess
import sys
from pathlib import Path

TEMPLATE_DIR = '.template'
PARENT_FILE = '__PARENT__'


def find_template_dirs(start=None):
    """Return template dirs from shallowest to deepest."""
    dirs = []
    cur = (Path.cwd() if start is None else Path(start)).resolve()
    while True:
        template = cur / TEMPLATE_DIR
        if template.is_dir():
            dirs.append(template)
            if not (template / PARENT_FILE).exists():
                break
        if cur.parent == cur:
            break
        cur = cur.parent
    return list(reversed(dirs))


def substitute(text, mapping):
    pattern = r'__(%s)__' % '|'.join(map(re.escape, mapping))
    return re.sub(pattern, lambda m: mapping[m.group(1)], text)


def instantiate(template_dirs, prob_dir, mapping):
    prob_dir.mkdir(parents=True)
    for template in template_dirs:
        for src in sorted(template.rglob('*')):
            rel = src.relative_to(template)
            if rel == Path(PARENT_FILE):
                continue
            dst = prob_dir / Path(*(substitute(part, mapping) for part in rel.parts))
            if src.is_dir():
                dst.mkdir(exist_ok=True)
            else:
                data = src.read_bytes()
                try:
                    data = substitute(data.decode(), mapping).encode()
                except UnicodeDecodeError:
                    pass  # binary file; copy as-is
                dst.write_bytes(data)
                shutil.copystat(src, dst)


def make_prob(name):
    prob_dir = Path(name)
    prob_name = prob_dir.name

    if prob_dir.exists():
        print(f"{name} already exists. Remove it and retry.")
        return False

    template_dirs = find_template_dirs(prob_dir.parent if prob_dir.parent != Path('') else None)
    if not template_dirs:
        print(f"No {TEMPLATE_DIR} directory found above {Path.cwd()}")
        return False

    instantiate(template_dirs, prob_dir, {'PROBLEM_NAME': prob_name})

    setup = prob_dir / 'setup'
    if setup.exists() and os.stat(setup).st_mode & stat.S_IXUSR:
        print("Running setup")
        try:
            subprocess.check_call(['./setup'], cwd=prob_dir)
        except subprocess.CalledProcessError as e:
            print(f"setup failed: {e}", file=sys.stderr)

    return True


def main():
    args = sys.argv[1:]
    if not args or args[0] in ('-h', '--help'):
        print(__doc__.strip(), file=sys.stderr)
        sys.exit(0 if args else 1)
    ok = True
    for name in args:
        ok &= make_prob(name)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
