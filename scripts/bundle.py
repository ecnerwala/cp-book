#!/usr/bin/env -S uv run --script
# /// script
# dependencies = ["competitive-verifier"]
# ///
"""Inline cp-book headers to produce a single submittable file.

Usage:
    scripts/bundle.py path/to/solution.cpp > submission.cpp

Any `#include "foo.hpp"` resolved from src/ (or relative to the including
file) is expanded in place, like `oj-bundle -I src`.

With --site verify_files.json, instead rewrites .competitive-verifier/bundled/
for every file listed there and registers the results as "bundled" sources,
for use when building the documentation site.

Runs via `uv run` (or plain python3 with competitive-verifier installed).
"""

import argparse
import json
import pathlib
import sys

from competitive_verifier.oj.languages.cplusplus_bundle import Bundler

ROOT = pathlib.Path(__file__).resolve().parent.parent


def bundle(path: pathlib.Path) -> bytes:
    bundler = Bundler(iquotes=[ROOT / "src"])
    bundler.update(path)
    return bundler.get()


def update_site(verify_files_json: pathlib.Path) -> None:
    data = json.loads(verify_files_json.read_text())
    bundled_root = ROOT / ".competitive-verifier" / "bundled"
    for path_str, file in data["files"].items():
        path = pathlib.Path(path_str)
        if path.suffix not in (".cpp", ".hpp"):
            continue
        try:
            bundled = bundle(path)
        except Exception as e:  # e.g. #pragma GCC target confuses -fpreprocessed -dD
            print(f"bundle failed for {path}: {e!r}", file=sys.stderr)
            continue
        dest = bundled_root / path
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(bundled)
        sources = [
            s for s in file.get("additonal_sources", [])
            if s["name"] not in ("bundled", "bundle error")
        ]
        sources.append({"name": "bundled", "path": str(dest.relative_to(ROOT))})
        file["additonal_sources"] = sources
    verify_files_json.write_text(json.dumps(data))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=pathlib.Path)
    parser.add_argument("--site", action="store_true",
                        help="treat PATH as verify_files.json and rebundle for the docs site")
    args = parser.parse_args()
    if args.site:
        update_site(args.path)
    else:
        sys.stdout.buffer.write(bundle(args.path))


if __name__ == "__main__":
    main()
