#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["competitive-verifier @ git+https://github.com/ecnerwala/competitive-verifier.git@cp-book-integration"]
# ///
"""Inline cp-book headers to produce a single submittable file.

Usage:
    scripts/bundle.py path/to/solution.cpp > submission.cpp
    scripts/bundle.py fft/series.hpp ds/seg_tree.hpp | xclip -selection clipboard
    scripts/bundle.py --minify fft/series.hpp > fft_series.min.cpp
    scripts/bundle.py --all -o dist/       # pregenerate all headers

Any `#include "foo.hpp"` resolved from src/ (or relative to the including
file) is expanded in place, like `oj-bundle -I src`. Header names relative
to src/ (e.g. `fft/series.hpp`) are looked up in src/. Multiple inputs are bundled into
one output with shared includes deduplicated.

--minify additionally strips comments (compiler-directed, via
`g++ -fpreprocessed -dD -E`), collapses whitespace, and packs lines,
keeping `#line` markers at file boundaries.

--all writes bundled (and minified) copies of every src/ header to
`<outdir>/bundled/` and `<outdir>/minified/`.

Runs via `uv run` (or plain python3 with the competitive-verifier fork
installed).
"""

import argparse
import pathlib
import shlex
import sys
from typing import Literal

from competitive_verifier.oj.languages.cplusplus_bundle import Bundler
from competitive_verifier.oj.languages.cplusplus_minify import (
    minify,
    raw_token_stream,
)

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
REPO_URL = "https://github.com/ecnerwala/cp-book"


def header_comment(args: list[str] | None = None) -> bytes:
    if args is None:
        args = sys.argv[1:]
    cmd = shlex.join(["scripts/bundle.py", *args])
    return f"// {REPO_URL} (`{cmd}`)\n".encode()


def resolve_input(path: pathlib.Path) -> pathlib.Path:
    if path.exists():
        return path
    if not path.is_absolute() and (SRC / path).exists():
        return SRC / path
    raise SystemExit(f"error: no such file: {path}")


def bundle(
    paths: list[pathlib.Path],
    *,
    level: Literal["light", "medium", "full"] | None,
    line_markers: bool = False,
) -> bytes:
    bundler = Bundler(iquotes=[SRC])
    for path in paths:
        bundler.update(resolve_input(path))
    code = bundler.get()
    if level is None:
        return code
    return minify(code, level=level, line_markers=line_markers)


def bundle_all(outdir: pathlib.Path, *, check: bool) -> None:
    headers = sorted(
        p for p in SRC.rglob("*.hpp") if not p.name.endswith(".test.hpp")
    )
    failures = []
    for header in headers:
        rel = header.relative_to(SRC)
        outputs = {}
        for name, level in (("bundled", None), ("minified", "medium")):
            dest = outdir / name / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            outputs[name] = bundle([header], level=level)
            args = (["-m"] if level else []) + [str(rel)]
            dest.write_bytes(header_comment(args) + outputs[name])
        if check and raw_token_stream(outputs["bundled"]) != raw_token_stream(
            outputs["minified"]
        ):
            failures.append(rel)
        print(rel, file=sys.stderr)
    if failures:
        raise SystemExit(
            "error: minified token stream differs for: "
            + " ".join(map(str, failures))
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "paths",
        type=pathlib.Path,
        nargs="*",
        help="files to bundle together (bare header names resolve from src/)",
    )
    parser.add_argument(
        "-m",
        "--minify",
        action="store_true",
        help="minify the bundled output (at --minify-level)",
    )
    parser.add_argument(
        "--minify-level",
        choices=["light", "medium", "full"],
        help="light: strip comments and blank lines only; medium (default): "
        "also compress whitespace, one statement per line; full: also pack "
        "statements onto shared lines up to 120 columns. Implies --minify",
    )
    parser.add_argument(
        "-o", "--output", type=pathlib.Path, help="output file (--all: output dir)"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="pregenerate bundled+minified copies of every src/ header",
    )
    parser.add_argument(
        "--line-markers",
        action="store_true",
        help="when minifying: keep exact #line directives (for in-repo "
        "compiles) "
        "instead of the default // file comments (safe to copy-paste)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="with --all: verify (via clang's raw lexer) that minification "
        "preserves every token",
    )
    args = parser.parse_args()

    if args.all:
        if args.paths:
            parser.error("--all takes no positional paths")
        bundle_all(args.output or ROOT / "dist", check=args.check)
        return
    if not args.paths:
        parser.error("no input files")
    level = args.minify_level or ("medium" if args.minify else None)
    code = header_comment() + bundle(
        args.paths, level=level, line_markers=args.line_markers
    )
    if args.output:
        args.output.write_bytes(code)
    else:
        sys.stdout.buffer.write(code)


if __name__ == "__main__":
    main()
