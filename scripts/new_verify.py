#!/usr/bin/env python3
"""Create a new verify/*.test.cpp stub.

Usage: scripts/new_verify.py <[folder/]problem[-tag]>

Creates verify/<arg>.test.cpp. The last path component, minus an optional
-tag suffix, is the yosupo problem name, e.g.
series/inv_of_formal_power_series-crt -> problem inv_of_formal_power_series.
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
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <[folder/]problem[-tag]>")
    arg = Path(sys.argv[1])
    problem = arg.name.split("-")[0]
    path = Path(__file__).resolve().parent.parent / "verify" / f"{arg}.test.cpp"
    if path.exists():
        sys.exit(f"{path} already exists")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(TEMPLATE.format(problem=problem))
    print(path)


if __name__ == "__main__":
    main()
