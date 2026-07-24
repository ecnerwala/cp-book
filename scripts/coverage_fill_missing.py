#!/usr/bin/env python3
"""Add zero-coverage entries to a gcovr JSON report for headers gcov never saw.

Template-only headers emit no code even when compiled into an instrumented
stub TU, so they never appear in gcov's output. This fills them in with
count-0 line entries so the docs render them as 0% covered. Line selection
approximates gcov: blank lines, comments, preprocessor directives, and
punctuation-only lines (``}``, ``{``, ``else``, labels...) are not counted.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

NONCODE_RE = re.compile(
    r"^(?:[{}();,]*|else|do|try|public:|private:|protected:|default:|case[^:]*:)$"
)


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
        if NONCODE_RE.fullmatch(line):
            continue
        lines.append(i)
    return lines


def main() -> None:
    report_path = Path(sys.argv[1])
    report = json.loads(report_path.read_text())
    present = {f["file"] for f in report["files"]}

    headers = subprocess.run(
        ["git", "ls-files", "src/**/*.hpp", "src/*.hpp"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split()

    added = 0
    for header in headers:
        if header in present:
            continue
        report["files"].append(
            {
                "file": header,
                "lines": [
                    {"line_number": n, "count": 0, "branches": []}
                    for n in executable_lines(Path(header))
                ],
                "functions": [],
            }
        )
        added += 1
        print(f"added zero-coverage entry: {header}")

    report_path.write_text(json.dumps(report))
    print(f"{added} headers added, {len(present)} already present")


if __name__ == "__main__":
    main()
