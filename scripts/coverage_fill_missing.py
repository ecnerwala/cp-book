#!/usr/bin/env python3
"""Add zero-coverage entries to a gcovr JSON report for code gcov never saw.

Uninstantiated templates emit no code even when compiled into an instrumented
stub TU, so gcov reports nothing for them: whole template-only headers are
absent from the report, and uninstantiated template functions in otherwise
covered headers are silently missing from the denominator. This fills both in
with count-0 line entries so the docs render them as uncovered.

Function extents come from universal-ctags; a function counts as
uninstantiated when the report has no line entries inside it. Line selection
within a function approximates gcov: blank lines, comments, preprocessor
directives, and punctuation-only lines (``}``, ``{``, ``else``, labels...)
are not counted, and ``assert(false)`` lines are excluded to match the
report's --exclude-lines-by-pattern.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

NONCODE_RE = re.compile(
    r"^(?:[{}();,]*|else|do|try|public:|private:|protected:|default:|case[^:]*:)$"
)
EXCLUDE_RE = re.compile(r".*assert\(false\).*")


def executable_lines(path: Path) -> list[int]:
    lines: list[int] = []
    in_block_comment = False
    for i, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.strip()
        if in_block_comment:
            end = line.find("*/")
            if end < 0:
                continue
            in_block_comment = False
            line = line[end + 2 :].strip()
        # strip block comments fully contained in the line
        line = re.sub(r"/\*.*?\*/", "", line).strip()
        if line.endswith("\\"):
            line = line[:-1].strip()
        start = line.find("/*")
        if start >= 0:
            in_block_comment = True
            line = line[:start].strip()
        line = line.split("//", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        if NONCODE_RE.fullmatch(line) or EXCLUDE_RE.fullmatch(line):
            continue
        lines.append(i)
    return lines


def function_ranges(path: Path) -> list[tuple[int, int]]:
    """Line ranges of every function definition, per universal-ctags."""
    tags = subprocess.run(
        [
            "ctags",
            "--output-format=json",
            "--fields=+ne-t",
            "--kinds-c++=f",
            "--language-force=c++",
            "-o",
            "-",
            str(path),
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    ranges = []
    for tag_line in tags.splitlines():
        tag = json.loads(tag_line)
        if tag.get("kind") == "function" and "end" in tag:
            ranges.append((tag["line"], tag["end"]))
    return ranges


def main() -> None:
    report_path = Path(sys.argv[1])
    report = json.loads(report_path.read_text())
    entries = {f["file"]: f for f in report["files"]}

    headers = subprocess.run(
        ["git", "ls-files", "src/**/*.hpp", "src/*.hpp"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split()

    files_added = lines_added = 0
    for header in headers:
        path = Path(header)
        entry = entries.get(header)
        if entry is None:
            # Header absent entirely: fill every executable line.
            report["files"].append(
                {
                    "file": header,
                    "lines": [
                        {"line_number": n, "count": 0, "branches": []}
                        for n in executable_lines(path)
                    ],
                    "functions": [],
                }
            )
            files_added += 1
            print(f"added zero-coverage entry: {header}")
            continue
        # Header present: fill executable lines of functions gcov reported
        # nothing for (uninstantiated templates).
        seen = {line["line_number"] for line in entry["lines"]}
        fill: set[int] = set()
        candidates = executable_lines(path)
        for start, end in function_ranges(path):
            if any(start <= n <= end for n in seen):
                continue
            fill.update(n for n in candidates if start <= n <= end)
        if fill:
            entry["lines"].extend(
                {"line_number": n, "count": 0, "branches": []} for n in sorted(fill)
            )
            entry["lines"].sort(key=lambda line: line["line_number"])
            lines_added += len(fill)
            print(f"added {len(fill)} zero-coverage lines: {header}")

    report_path.write_text(json.dumps(report))
    print(f"{files_added} headers and {lines_added} function lines added")


if __name__ == "__main__":
    main()
