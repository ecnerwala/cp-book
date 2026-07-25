#!/usr/bin/env python3
"""Create a new verify/*.test.cpp stub.

Usage: scripts/new_verify.py <problem> [tag]

<problem> may include a folder, e.g. inv_of_formal_power_series or
series/inv_of_formal_power_series; the last component is the yosupo
problem name. Creates verify/<folder>/<problem>[-tag].test.cpp.
"""

import sys
from pathlib import Path

TEMPLATE = """\
// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/{problem}

#include <bits/stdc++.h>
#include <cassert>

int main() {{
\tstd::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

\treturn 0;
}}
"""


def main():
    if len(sys.argv) not in (2, 3):
        sys.exit(f"Usage: {sys.argv[0]} <[folder/]problem> [tag]")
    problem_path = Path(sys.argv[1])
    tag = sys.argv[2] if len(sys.argv) == 3 else None
    problem = problem_path.name
    name = f"{problem}-{tag}" if tag else problem
    path = Path(__file__).resolve().parent.parent / "verify" / problem_path.parent / f"{name}.test.cpp"
    if path.exists():
        sys.exit(f"{path} already exists")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(TEMPLATE.format(problem=problem))
    print(path)


if __name__ == "__main__":
    main()
