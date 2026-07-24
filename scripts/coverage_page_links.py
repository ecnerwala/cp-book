#!/usr/bin/env python3
"""Map source files to their gcovr details pages for the docs site.

Writes <jekyll-dir>/_data/coverage_pages.yml mapping each covered source path
to its gcovr details page href (gcovr names the pages
``index.<basename>.<md5 of path>.html``), so templates can link to them.

Usage: coverage_page_links.py --coverage-dir coverage --jekyll-dir _jekyll
"""

import argparse
import hashlib
import json
import pathlib
import sys


def details_page(coverage_dir: pathlib.Path, filename: str) -> str | None:
    digest = hashlib.md5(filename.encode()).hexdigest()
    name = f"index.{pathlib.PurePosixPath(filename).name}.{digest}.html"
    return name if (coverage_dir / name).exists() else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverage-dir", type=pathlib.Path, required=True)
    parser.add_argument("--jekyll-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()

    summary = json.loads((args.coverage_dir / "summary.json").read_text())

    lines = []
    for f in summary["files"]:
        filename = f["filename"]
        name = details_page(args.coverage_dir, filename)
        if name:
            lines.append(f'"{filename}": "/coverage/{name}"\n')

    data_dir = args.jekyll_dir / "_data"
    data_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "coverage_pages.yml").write_text("".join(lines))
    print(f"coverage_page_links: mapped {len(lines)} files", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
