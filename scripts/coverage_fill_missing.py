#!/usr/bin/env python3
"""Add zero-coverage entries to a gcovr JSON report for code gcov never saw.

Uninstantiated templates emit no code even when compiled into an instrumented
stub TU, so gcov reports nothing for them: whole template-only headers are
absent from the report, and uninstantiated template functions in otherwise
covered headers are silently missing from the denominator. This fills both in
with count-0 line entries so the docs render them as uncovered.

Function extents and executable lines come from a tree-sitter C++ parse
(pip: tree_sitter + tree_sitter_cpp). A function counts as uninstantiated
when the report has no line entries inside it; its fill lines approximate
gcov's notion of executable lines: the signature line (function entry block),
statements, declarations, constructor field initializers, non-constant loop/if
conditions, and call argument lists. ``assert(false)`` lines are excluded to
match the report's --exclude-lines-by-pattern.

Benchmarked against real gcov data inside instantiated functions, this
selection agrees on >99% of lines (vs ~95% for the previous regex heuristic).
"""

import json
import re
import subprocess
import sys
from pathlib import Path

import tree_sitter_cpp
from tree_sitter import Language, Parser

PARSER = Parser(Language(tree_sitter_cpp.language()))

FUNC_NODES = ("function_definition", "lambda_expression")
STMT_NODES = {
    "expression_statement",
    "return_statement",
    "if_statement",
    "for_statement",
    "while_statement",
    "do_statement",
    "switch_statement",
    "break_statement",
    "continue_statement",
    "goto_statement",
    "co_return_statement",
    "co_yield_statement",
    "throw_statement",
    "for_range_loop",
    "init_declarator",
    "declaration",
    "field_initializer",
    "argument_list",
}
EXCLUDE_RE = re.compile(r".*assert\(false\).*")


def parse_header(path: Path) -> tuple[list[tuple[int, int]], list[set[int]]]:
    """Function line ranges and per-function executable-line sets."""
    source_lines = path.read_text().splitlines()
    tree = PARSER.parse(path.read_bytes())
    ranges: list[tuple[int, int]] = []
    line_sets: list[set[int]] = []

    def add_line(lines: set[int], row: int) -> None:
        if not EXCLUDE_RE.fullmatch(source_lines[row].strip()):
            lines.add(row + 1)

    def walk(node, lines: set[int] | None) -> None:
        for child in node.children:
            kind = child.type
            if kind in FUNC_NODES:
                ranges.append((child.start_point[0] + 1, child.end_point[0] + 1))
                line_sets.append(set())
                # gcov attaches the function entry block to the signature line
                add_line(line_sets[-1], child.start_point[0])
                walk(child, line_sets[-1])
                continue
            if lines is not None:
                if kind in STMT_NODES:
                    add_line(lines, child.start_point[0])
                elif kind == "condition_clause":
                    # constant conditions (while (true)) emit no code
                    if child.text.decode() not in ("(true)", "(false)"):
                        add_line(lines, child.start_point[0])
            walk(child, lines)

    walk(tree.root_node, None)
    return ranges, line_sets


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
    headers = [h for h in headers if not h.endswith(".test.hpp")]

    files_added = lines_added = 0
    for header in headers:
        entry = entries.get(header)
        if entry is None:
            entry = {"file": header, "lines": [], "functions": []}
            report["files"].append(entry)
            files_added += 1
        seen = {line["line_number"] for line in entry["lines"]}
        fill: set[int] = set()
        for (start, end), lines in zip(*parse_header(Path(header))):
            if any(start <= n <= end for n in seen):
                continue  # instantiated: gcov data is authoritative
            fill.update(lines - seen)
        if fill:
            entry["lines"].extend(
                {"line_number": n, "count": 0, "branches": []} for n in sorted(fill)
            )
            entry["lines"].sort(key=lambda line: line["line_number"])
            lines_added += len(fill)
            print(f"added {len(fill)} zero-coverage lines: {header}")

    report_path.write_text(json.dumps(report))
    print(f"{files_added} headers and {lines_added} lines added")


if __name__ == "__main__":
    main()
