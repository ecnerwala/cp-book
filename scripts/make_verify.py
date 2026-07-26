#!/usr/bin/env python3
"""Create a new verify/*.test.cpp stub.

Usage: scripts/make_verify.py <[folder/]problem[-tag]>

Creates verify/<arg>.test.cpp from contest/template/__PROBLEM_NAME__.cpp
with the competitive-verifier PROBLEM header prepended. The last path
component, minus an optional -tag suffix, is the yosupo problem name, e.g.
series/inv_of_formal_power_series-crt -> problem inv_of_formal_power_series.
"""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <[folder/]problem[-tag]>")
    arg = Path(sys.argv[1])
    problem = arg.name.split("-")[0]
    path = REPO / "verify" / f"{arg}.test.cpp"
    if path.exists():
        sys.exit(f"{path} already exists")
    template = (REPO / "contest" / "template" / "__PROBLEM_NAME__.cpp").read_text()
    header = f"// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/{problem}\n\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(header + template)
    print(path)


if __name__ == "__main__":
    main()
