#!/usr/bin/env -S uv run
"""Inline cp-book headers to produce a single submittable file.

Usage:
    scripts/bundle.py path/to/solution.cpp > submission.cpp
    scripts/bundle.py fft/series.hpp ds/seg_tree.hpp | xclip -selection clipboard  # or wl-copy
    scripts/bundle.py --minify fft/series.hpp > fft_series.min.cpp
    scripts/bundle.py --all -o dist/       # pregenerate all headers
    scripts/bundle.py --at 12 ds/seg_tree.hpp  # to paste at line 12
    :'<,'>!scripts/bundle.py --at 12 -     # in vim, on lines 12 onwards

Any `#include "foo.hpp"` resolved from src/ (or relative to the including
file) is expanded in place, like `oj-bundle -I src`. Header names relative
to src/ (e.g. `fft/series.hpp`) are looked up in src/, and inputs under src/
are recorded by that name. Multiple inputs are bundled into
one output with shared includes deduplicated and system includes hoisted
into one block at the top.

--minify additionally collapses the standard includes into
`#include <bits/stdc++.h>` plus `#include <cassert>` (not part of
`<bits/stdc++.h>` in recent g++), strips comments (compiler-directed, via
`g++ -fpreprocessed -dD -E`), collapses whitespace, and packs lines,
keeping `#line` markers at file boundaries.

--all writes bundled (and minified) copies of every src/ header to
`<outdir>/bundled/` and `<outdir>/minified/`.

The output is wrapped in a single fold so the pasted block can be
collapsed in an editor: an `#if 1` / `#endif` pair (treesitter and other
syntax-aware folding) carrying `// region ...` / `// endregion` comments
(IntelliJ region folding).

--at LINE says the output will be pasted with its first line at line LINE of
the main source file; the fold then ends with a `#line ... __BASE_FILE__`, so
code after the pasted block keeps its own file and line numbers.

A `-` input reads a selection from stdin, for use as an editor filter. The
selection holds bundle regions written by this script (their recorded
arguments are reused) and `#include "foo.hpp"` lines (added as inputs), so
adding a header to an existing bundle is: write its `#include` next to the
region, then filter both through `scripts/bundle.py -`.

Runs via `uv run` (or plain python3 with the competitive-verifier fork
installed).
"""

import argparse
import pathlib
import re
import shlex
import sys
from typing import Literal

from competitive_verifier.oj.languages.cplusplus_bundle import Bundler
from competitive_verifier.oj.languages.cplusplus_minify import (
    MinifyCheckError,
    minify,
)

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
REPO_URL = "https://github.com/ecnerwala/cp-book"


def wrap_fold(code: bytes, args: list[str], at: int | None = None) -> bytes:
    cmd = shlex.join(["scripts/bundle.py", *args])
    head = f"#if 1 // region {REPO_URL} (`{cmd}`)\n".encode()
    tail = b"#endif // endregion\n"
    if at is not None:
        # The #line directive itself is line `at + 1 + code lines`.
        line = at + 2 + code.count(b"\n")
        code += f"#line {line} __BASE_FILE__\n".encode()
    return head + code + tail


REGION_HEAD = re.compile(r"#if 1 // region \S+ \(`(.*)`\)")
REGION_TAIL = "#endif // endregion"
QUOTED_INCLUDE = re.compile(r'\s*#\s*include\s*"([^"]*)"\s*')


def read_selection(text: str) -> list[str]:
    args = []
    in_region = False
    for lineno, line in enumerate(text.splitlines(), 1):
        if in_region:
            in_region = line != REGION_TAIL
        elif m := REGION_HEAD.fullmatch(line):
            args += shlex.split(m.group(1))[1:]
            in_region = True
        elif m := QUOTED_INCLUDE.fullmatch(line):
            args.append(m.group(1))
        elif line.strip():
            raise SystemExit(
                f"error: stdin line {lineno}: expected a bundle region or "
                f'#include "...", got: {line}'
            )
    if in_region:
        raise SystemExit(f"error: stdin: unterminated region (no `{REGION_TAIL}`)")
    return args


def resolve_input(path: pathlib.Path) -> pathlib.Path:
    if path.exists():
        return path
    if not path.is_absolute() and (SRC / path).exists():
        return SRC / path
    raise SystemExit(f"error: no such file: {path}")


def input_name(path: pathlib.Path) -> pathlib.Path:
    """The name recorded in the region header: src/ headers by bare name."""
    resolved = resolve_input(path).resolve()
    if resolved.is_relative_to(SRC):
        return resolved.relative_to(SRC)
    return path


MINIFY_PRELUDE = ["bits/stdc++.h", "cassert"]


def bundle(
    paths: list[pathlib.Path],
    *,
    level: Literal["light", "medium", "full"] | None,
    line_markers: bool = False,
    check: bool = False,
) -> bytes:
    bundler = Bundler(
        iquotes=[SRC],
        prelude_includes=MINIFY_PRELUDE if level else [],
        hoist_system_includes=True,
    )
    for path in paths:
        bundler.update(resolve_input(path))
    code = bundler.get()
    if level is None:
        return code
    return minify(code, level=level, line_markers=line_markers, check=check)


def bundle_all(outdir: pathlib.Path, *, check: bool) -> None:
    headers = sorted(
        p for p in SRC.rglob("*.hpp") if not p.name.endswith(".test.hpp")
    )
    failures = []
    for header in headers:
        rel = header.relative_to(SRC)
        for name, level in (("bundled", None), ("minified", "medium")):
            dest = outdir / name / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            try:
                code = bundle([header], level=level, check=check)
            except MinifyCheckError:
                failures.append(rel)
                continue
            args = (["-m"] if level else []) + [str(rel)]
            dest.write_bytes(wrap_fold(code, args))
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
        "--at",
        type=int,
        metavar="LINE",
        help="the output will be pasted with its first line at LINE of the "
        "main source file: end with a #line back to it",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="with --minify: verify (via clang's raw lexer) that minification "
        "preserves every token",
    )
    argv = sys.argv[1:]
    if "-" in argv:
        argv = read_selection(sys.stdin.read()) + [a for a in argv if a != "-"]
    args = parser.parse_intermixed_args(argv)

    if args.all:
        if args.paths:
            parser.error("--all takes no positional paths")
        if args.at:
            parser.error("--all does not take --at")
        bundle_all(args.output or ROOT / "dist", check=args.check)
        return
    if not args.paths:
        parser.error("no input files")
    args.paths = list(dict.fromkeys(map(input_name, args.paths)))
    level = args.minify_level or ("medium" if args.minify else None)
    if args.minify_level:
        cmd = ["--minify-level", args.minify_level]
    else:
        cmd = ["-m"] if args.minify else []
    if args.line_markers:
        cmd.append("--line-markers")
    cmd += map(str, args.paths)
    code = wrap_fold(
        bundle(
            args.paths,
            level=level,
            line_markers=args.line_markers,
            check=args.check,
        ),
        cmd,
        at=args.at,
    )
    if args.output:
        args.output.write_bytes(code)
    else:
        sys.stdout.buffer.write(code)


if __name__ == "__main__":
    main()
