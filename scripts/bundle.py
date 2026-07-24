#!/usr/bin/env -S uv run --script
# /// script
# dependencies = ["competitive-verifier"]
# ///
"""Inline cp-book headers to produce a single submittable file.

Usage:
    scripts/bundle.py path/to/solution.cpp > submission.cpp

Any `#include "foo.hpp"` resolved from src/ (or relative to the including
file) is expanded in place, like `oj-bundle -I src`.

Runs via `uv run` (or plain python3 with competitive-verifier installed).
"""

import argparse
import pathlib
import sys

from competitive_verifier.oj.languages.cplusplus_bundle import Bundler

ROOT = pathlib.Path(__file__).resolve().parent.parent


def bundle(path: pathlib.Path) -> bytes:
    bundler = Bundler(iquotes=[ROOT / "src"])
    bundler.update(path)
    return bundler.get()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=pathlib.Path)
    args = parser.parse_args()
    sys.stdout.buffer.write(bundle(args.path))


if __name__ == "__main__":
    main()
