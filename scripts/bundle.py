#!/usr/bin/env -S uv run
"""Inline cp-book headers to produce a single submittable file.

Usage:
    scripts/bundle.py path/to/solution.cpp > submission.cpp
    scripts/bundle.py fft/series.hpp ds/seg_tree.hpp | xclip -selection clipboard  # or wl-copy
    scripts/bundle.py --minify fft/series.hpp > fft_series.min.cpp
    scripts/bundle.py --all -o dist/       # pregenerate all headers

A thin wrapper over cpp-bundle and cpp-minify
(https://github.com/ecnerwala/cpp-bundle, installed into the uv environment
by pyproject.toml). Every `#include "foo.hpp"` (resolved relative to the
including file, then src/) is expanded in place; `#include <...>` lines stay,
deduplicated. Header names relative to src/ (e.g. `fft/series.hpp`) are
looked up in src/. Multiple inputs are bundled into one output.

--minify additionally puts `#include <bits/stdc++.h>` and `#include <cassert>`
(not part of `<bits/stdc++.h>` in recent g++) first, dropping the standard
includes they cover, strips comments and collapses whitespace. The minified
token stream is checked against the input.

--all writes bundled (and minified) copies of every src/ header to
`<outdir>/bundled/` and `<outdir>/minified/`.

The output is wrapped in a single fold so the pasted block can be
collapsed in an editor: an `#if 1` / `#endif` pair (treesitter and other
syntax-aware folding) carrying `// region ...` / `// endregion` comments
(IntelliJ region folding).

Runs via `uv run` (or plain python3 with cpp-bundle / cpp-minify on PATH).
"""

import argparse
import pathlib
import shlex
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
REPO_URL = "https://github.com/ecnerwala/cp-book"

CLANG_ARGS = ["-std=c++23", "-I", str(SRC)]
MINIFY_PRELUDE = ["bits/stdc++.h", "cassert"]


def wrap_fold(code: bytes, args: list[str] | None = None) -> bytes:
    if args is None:
        args = sys.argv[1:]
    cmd = shlex.join(["scripts/bundle.py", *args])
    head = f"#if 1 // region {REPO_URL} (`{cmd}`)\n".encode()
    tail = b"#endif // endregion\n"
    return head + code + tail


def resolve_input(path: pathlib.Path) -> pathlib.Path:
    if path.exists():
        return path
    if not path.is_absolute() and (SRC / path).exists():
        return SRC / path
    raise SystemExit(f"error: no such file: {path}")


def tool(name: str) -> str:
    venv_bin = pathlib.Path(sys.executable).parent / name
    found = str(venv_bin) if venv_bin.exists() else shutil.which(name)
    if found is None:
        raise SystemExit(f"error: {name} not found; run `uv sync` (see pyproject.toml)")
    return found


def bundle(paths: list[pathlib.Path], *, minify: bool) -> bytes:
    cmd = [tool("cpp-bundle"), *CLANG_ARGS]
    if minify:
        for header in MINIFY_PRELUDE:
            cmd += ["-include", header]
    cmd += [str(resolve_input(path)) for path in paths]
    code = subprocess.run(cmd, check=True, stdout=subprocess.PIPE).stdout
    if minify:
        code = subprocess.run(
            [tool("cpp-minify"), "--check"], input=code, check=True, stdout=subprocess.PIPE
        ).stdout
    return code


def bundle_all(outdir: pathlib.Path) -> None:
    headers = sorted(
        p for p in SRC.rglob("*.hpp") if not p.name.endswith(".test.hpp")
    )
    failures = []
    for header in headers:
        rel = header.relative_to(SRC)
        for name, minify in (("bundled", False), ("minified", True)):
            dest = outdir / name / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            try:
                code = bundle([header], minify=minify)
            except subprocess.CalledProcessError:
                failures.append(rel)
                continue
            args = (["-m"] if minify else []) + [str(rel)]
            dest.write_bytes(wrap_fold(code, args))
        print(rel, file=sys.stderr)
    if failures:
        raise SystemExit("error: bundling failed for: " + " ".join(map(str, failures)))


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
        "-m", "--minify", action="store_true", help="minify the bundled output"
    )
    parser.add_argument(
        "-o", "--output", type=pathlib.Path, help="output file (--all: output dir)"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="pregenerate bundled+minified copies of every src/ header",
    )
    args = parser.parse_args()

    if args.all:
        if args.paths:
            parser.error("--all takes no positional paths")
        bundle_all(args.output or ROOT / "dist")
        return
    if not args.paths:
        parser.error("no input files")
    try:
        code = wrap_fold(bundle(args.paths, minify=args.minify))
    except subprocess.CalledProcessError as e:
        raise SystemExit(e.returncode)
    if args.output:
        args.output.write_bytes(code)
    else:
        sys.stdout.buffer.write(code)


if __name__ == "__main__":
    main()
