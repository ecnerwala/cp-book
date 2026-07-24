#!/usr/bin/env python3
"""Add per-file coverage sections to the generated docs site.

For each source file with coverage data, appends a "Coverage" section to its
generated page in the Jekyll source dir, linking to the file's gcovr details
page with line/branch percentages.

Usage: coverage_links.py --coverage-dir coverage --jekyll-dir _jekyll
"""

import argparse
import json
import pathlib
import re
import sys


def details_page_map(index_html: pathlib.Path) -> dict[str, str]:
    """Map source path -> gcovr details page name, parsed from index.html."""
    mapping: dict[str, str] = {}
    for href, text in re.findall(
        r'<a href="([^"]+\.html)"[^>]*>([^<]+)</a>', index_html.read_text()
    ):
        mapping[text.strip()] = href
    return mapping


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverage-dir", type=pathlib.Path, required=True)
    parser.add_argument("--jekyll-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()

    summary = json.loads((args.coverage_dir / "summary.json").read_text())
    pages = details_page_map(args.coverage_dir / "index.html")

    count = 0
    for f in summary["files"]:
        filename = f["filename"]
        page = args.jekyll_dir / (filename + ".md")
        if not page.exists():
            continue
        href = pages.get(filename)
        link = f"{{{{ site.baseurl }}}}/coverage/{href}" if href else None
        parts = []
        for key, label in (("line_percent", "lines"), ("branch_percent", "branches")):
            percent = f.get(key)
            if percent is not None:
                parts.append(f"{percent:.1f}% {label}")
        if not parts:
            continue
        text = ", ".join(parts)
        body = f"[{text}]({link})" if link else text
        with page.open("a") as fp:
            fp.write(f"\n**Coverage:** {body}\n")
        count += 1
    print(f"coverage_links: annotated {count} pages", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
