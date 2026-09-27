#!/usr/bin/env python3
"""Parse Polygon/olymp.sty-generated statement PDFs and send problems to
Competitive Companion-compatible tools (CP Editor, cph, cpbooster, ...).

Usage:
    polygon_pdf_companion.py statements.pdf            # parse + POST to localhost
    polygon_pdf_companion.py statements.pdf --dump     # print JSON instead
    polygon_pdf_companion.py statements.pdf -p A -p C  # only some problems
    polygon_pdf_companion.py statements.pdf --port 27121
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import re
import sys
import urllib.error
import urllib.request
import uuid
from typing import Any

import pdfplumber

# Default ports from competitive-companion's src/hosts/hosts.ts
DEFAULT_PORTS = [
    1327,  # cpbooster
    4244,  # Hightail
    6174,  # Mind Sport
    10042,  # acmX
    10043,  # Caide
    10045,  # CP Editor
    27121,  # Competitive Programming Helper (cph)
]

PROBLEM_RE = re.compile(r"^Problem\s+([A-Z][0-9]?|\d+)[.:]\s+(.*\S)\s*$")
# "Problem A" with the title on the following line (ICPC EC-Final, UCup Finals)
PROBLEM_BARE_RE = re.compile(r"^Problem\s+([A-Z][0-9]?|\d+)\s*$")
# "A \u2013 Three Castles" (EUC style), possibly with a trailing "Page x of y"
PROBLEM_DASH_RE = re.compile(
    r"^([A-Z][0-9]?)\s+[\u2013\u2014\u2212-]\s+(.+?)(?:\s+Page \d+ of \d+)?\s*$"
)
SAMPLE_HEADING_RE = re.compile(
    r"^\s*(Examples?|Sample Input|Sample Output|Sample \d)\b"
)
HEADER_FIELD_RE = re.compile(
    r"^(Input file|Output file|Time limit|Memory limit)\s*:\s*(.*\S)\s*$"
)
TIME_RE = re.compile(
    r"([0-9]+(?:[.,][0-9]+)?)\s*(m(?:illi)?)?\s*(?:second|sec\b|s\b)", re.IGNORECASE
)
MEM_RE = re.compile(
    r"([0-9]+(?:[.,][0-9]+)?)\s*(?:(giga|gibi|mega|mebi|kilo|kibi)bytes?|(G|M|K)i?B)",
    re.IGNORECASE,
)


def classify_header(row: list[str | None]) -> tuple[list[int], list[int]] | None:
    """Return (input column indices, output column indices) if this row looks
    like a sample-table header ("standard input | standard output", possibly
    with extra columns like "explanation"/"Notes", or run-twice communication
    headers like "Alice's Input | Alice's Output | Bob's Input | Bob's Output").
    """
    inputs, outputs = [], []
    for i, cell in enumerate(row):
        c = (cell or "").strip().lower()
        if "input" in c:
            inputs.append(i)
        elif "output" in c:
            outputs.append(i)
    if not inputs or not outputs:
        return None
    return inputs, outputs


@dataclasses.dataclass
class Problem:
    letter: str
    title: str
    group: str
    url: str
    input_file: str = "standard input"
    output_file: str = "standard output"
    time_limit_ms: int = 1000
    memory_limit_mb: int = 256
    interactive: bool = False
    tests: list[dict[str, str]] = dataclasses.field(default_factory=list)
    saw_headed_table: bool = False

    def to_companion(self, batch_id: str, batch_size: int) -> dict[str, Any]:
        def io_spec(value: str, kind: str) -> dict[str, str]:
            if value in ("standard input", "standard output", "stdin", "stdout"):
                return {"type": "stdin" if kind == "input" else "stdout"}
            return {"type": "file", "fileName": value}

        return {
            "name": f"{self.letter}. {self.title}",
            "group": self.group,
            "url": self.url,
            "interactive": self.interactive,
            "memoryLimit": self.memory_limit_mb,
            "timeLimit": self.time_limit_ms,
            "tests": self.tests,
            "testType": "single",
            "input": io_spec(self.input_file, "input"),
            "output": io_spec(self.output_file, "output"),
            "languages": {
                "java": {
                    "mainClass": "Main",
                    "taskClass": self.letter
                    + re.sub(r"[^a-zA-Z0-9]", "", self.title.title()),
                }
            },
            "batch": {"id": batch_id, "size": batch_size},
        }


def parse_time_limit_ms(value: str) -> int | None:
    m = TIME_RE.search(value)
    if not m:
        return None
    amount = float(m.group(1).replace(",", "."))
    return round(amount if m.group(2) else amount * 1000)


def parse_memory_limit_mb(value: str) -> int | None:
    m = MEM_RE.search(value)
    if not m:
        return None
    amount = float(m.group(1).replace(",", "."))
    unit = (m.group(2) or m.group(3) or "mega").lower()
    if unit in ("giga", "gibi", "g"):
        amount *= 1024
    elif unit in ("kilo", "kibi", "k"):
        amount /= 1024
    return round(amount)


def clean_cell(text: str | None) -> str:
    if text is None:
        return ""
    lines = [line.rstrip() for line in text.split("\n")]
    while lines and not lines[-1]:
        lines.pop()
    return "\n".join(lines)


def extract_page_tables(page: Any) -> list[dict[str, Any]]:
    tables = []
    # Generous intersection tolerance: some Polygon PDFs draw table rules with
    # slightly mismatched endpoints, which splits the sample table in half.
    settings = {"snap_tolerance": 3, "join_tolerance": 3, "intersection_tolerance": 10}
    for table in page.find_tables(settings):
        rows = table.extract()
        if not rows or max(len(r) for r in rows) < 2:
            continue
        tables.append({"bbox": table.bbox, "rows": rows, "page_height": page.height})
    return tables


def gather_columns(row: list[str | None], cols: list[int]) -> str:
    parts = [clean_cell(row[i]) for i in cols if i < len(row)]
    return "\n".join(p for p in parts if p)


def add_sample_rows(
    problem: Problem,
    rows: list[list[str | None]],
    input_cols: list[int],
    output_cols: list[int],
    transcript: bool,
) -> None:
    if transcript:
        # Interactive/communication transcript: interleave all rows into a
        # single test, keeping judge/participant order within each column.
        inp = [s for row in rows if (s := gather_columns(row, input_cols))]
        out = [s for row in rows if (s := gather_columns(row, output_cols))]
        if inp or out:
            problem.tests.append(
                {"input": "\n".join(inp) + "\n", "output": "\n".join(out) + "\n"}
            )
        return
    for row in rows:
        inp = gather_columns(row, input_cols)
        out = gather_columns(row, output_cols)
        if not inp and not out:
            continue
        problem.tests.append({"input": inp + "\n", "output": out + "\n"})


def merge_continuation(
    problem: Problem,
    rows: list[list[str | None]],
    input_cols: list[int],
    output_cols: list[int],
    transcript: bool,
) -> None:
    # Table continued on the next page without a repeated header: glue the
    # first continuation row onto the last test, remaining rows are new tests.
    if not problem.tests or not rows or transcript:
        add_sample_rows(problem, rows, input_cols, output_cols, transcript)
        return
    first, rest = rows[0], rows[1:]
    last = problem.tests[-1]
    inp = gather_columns(first, input_cols)
    out = gather_columns(first, output_cols)
    if inp:
        last["input"] = last["input"].rstrip("\n") + "\n" + inp + "\n"
    if out:
        last["output"] = last["output"].rstrip("\n") + "\n" + out + "\n"
    add_sample_rows(problem, rest, input_cols, output_cols, transcript)


def parse_pdf(path: str) -> list[Problem]:
    problems: list[Problem] = []
    current: Problem | None = None
    group = ""
    # Track whether the previous page ended with an open sample table so a
    # header-less table at the top of the next page is treated as its
    # continuation.
    pending_continuation = False
    last_cols: tuple[list[int], list[int], bool] = ([0], [1], False)

    with pdfplumber.open(path) as pdf:
        for page in pdf.pages:
            text = page.extract_text() or ""
            lines = text.split("\n")
            if not group and lines:
                group = lines[0].strip()

            new_problem = None
            for idx, line in enumerate(lines[:8]):
                stripped = line.strip()
                m = PROBLEM_RE.match(stripped)
                if m:
                    new_problem = (m.group(1), m.group(2))
                    break
                m = PROBLEM_BARE_RE.match(stripped)
                if m and idx + 1 < len(lines) and lines[idx + 1].strip():
                    new_problem = (m.group(1), lines[idx + 1].strip())
                    break
                m = PROBLEM_DASH_RE.match(stripped)
                if m and idx < 2:
                    new_problem = (m.group(1), m.group(2))
                    break
            if new_problem and (current is None or new_problem[0] != current.letter):
                letter, title = new_problem
                current = Problem(
                    letter=letter,
                    title=title,
                    group=group,
                    url=f"file://{path}#{letter}",
                )
                problems.append(current)
                pending_continuation = False

            if current is None:
                continue

            if new_problem:
                for line in lines:
                    m = HEADER_FIELD_RE.match(line.strip())
                    if not m:
                        continue
                    field, value = m.group(1), m.group(2)
                    if field == "Input file":
                        current.input_file = value
                    elif field == "Output file":
                        current.output_file = value
                    elif field == "Time limit":
                        current.time_limit_ms = (
                            parse_time_limit_ms(value) or current.time_limit_ms
                        )
                    elif field == "Memory limit":
                        current.memory_limit_mb = (
                            parse_memory_limit_mb(value) or current.memory_limit_mb
                        )

            if re.search(r"\binteractive problem\b|Interaction Protocol", text):
                current.interactive = True

            has_sample_heading = any(SAMPLE_HEADING_RE.match(line) for line in lines)
            page_ends_open = False
            for table in extract_page_tables(page):
                rows = table["rows"]
                cols = classify_header(rows[0])
                if cols is not None:
                    input_cols, output_cols = cols
                    transcript = current.interactive or len(input_cols) > 1
                    add_sample_rows(
                        current, rows[1:], input_cols, output_cols, transcript
                    )
                    last_cols = (input_cols, output_cols, transcript)
                    current.saw_headed_table = True
                elif (
                    pending_continuation
                    and table["bbox"][1] < table["page_height"] * 0.25
                ):
                    merge_continuation(current, rows, *last_cols)
                elif (
                    has_sample_heading
                    and not current.saw_headed_table
                    and all(len(r) == 2 for r in rows)
                ):
                    # Templates without in-table headers (ICPC EC-Final, UCup
                    # Finals, EUC): each bordered 2-column box under a
                    # "Sample Input N"/"Examples" heading is one test.
                    add_sample_rows(current, rows, [0], [1], current.interactive)
                else:
                    continue
                # Table reaching close to the bottom of the page may continue.
                page_ends_open = table["bbox"][3] > table["page_height"] * 0.9
            pending_continuation = page_ends_open

    return problems


def send_to_companion(payload: dict[str, Any], ports: list[int]) -> bool:
    body = json.dumps(payload).encode()
    delivered = False
    for port in ports:
        req = urllib.request.Request(
            f"http://localhost:{port}/",
            data=body,
            headers={"Content-Type": "application/json"},
        )
        try:
            urllib.request.urlopen(req, timeout=2)
            delivered = True
            print(f"  sent to localhost:{port}", file=sys.stderr)
        except (urllib.error.URLError, OSError):
            pass
    return delivered


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pdf", help="Polygon/olymp.sty-generated statements PDF")
    ap.add_argument(
        "-p",
        "--problem",
        action="append",
        metavar="LETTER",
        help="only send these problem letters (repeatable)",
    )
    ap.add_argument(
        "--port",
        type=int,
        action="append",
        help="destination port(s); default: all known Companion ports",
    )
    ap.add_argument("--dump", action="store_true", help="print JSON, don't POST")
    args = ap.parse_args()

    problems = parse_pdf(args.pdf)
    if args.problem:
        wanted = {p.upper() for p in args.problem}
        problems = [p for p in problems if p.letter.upper() in wanted]
    if not problems:
        print("no problems found", file=sys.stderr)
        return 1

    batch_id = str(uuid.uuid4())
    payloads = [p.to_companion(batch_id, len(problems)) for p in problems]

    if args.dump:
        json.dump(payloads, sys.stdout, indent=2)
        print()
        return 0

    ports = args.port or DEFAULT_PORTS
    ok = True
    for problem, payload in zip(problems, payloads):
        print(
            f"{problem.letter}. {problem.title} "
            f"({len(problem.tests)} tests, {problem.time_limit_ms}ms, "
            f"{problem.memory_limit_mb}MB"
            f"{', interactive' if problem.interactive else ''})",
            file=sys.stderr,
        )
        ok &= send_to_companion(payload, ports)
    if not ok:
        print(
            "warning: some problems were not accepted by any listener", file=sys.stderr
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
