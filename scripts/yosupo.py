#!/usr/bin/env -S uv run --script
# /// script
# dependencies = ["competitive-verifier"]
# ///
"""Helpers for Library Checker (https://judge.yosupo.jp/).

Usage:
    scripts/yosupo.py login                     # store credentials (~/.config/yosupo/)
    scripts/yosupo.py status                    # local verify/ files vs. judge AC status
    scripts/yosupo.py submissions [--problem P] # list your recent submissions
    scripts/yosupo.py submit verify/foo.test.cpp [--wait]

`submit` bundles the file (inlining src/ headers) and posts it to the judge;
`--wait` polls until judging finishes. Reads (status/submissions) work without
logging in if you pass --user; submitting requires `login` first.
"""

import argparse
import concurrent.futures
import difflib
import getpass
import json
import pathlib
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import bundle as bundle_mod

API_URL = "https://v3.api.judge.yosupo.jp"
FIREBASE_API_KEY = "AIzaSyCmpkoMVbKRDm2H0MJHB0iZ43uQtSqiLV0"
CRED_PATH = pathlib.Path.home() / ".config" / "yosupo" / "credentials.json"

PROBLEM_RE = re.compile(
    r"competitive-verifier:\s*PROBLEM\s+https://judge\.yosupo\.jp/problem/(\S+)"
)


def http_json(
    url: str,
    *,
    data: dict[str, Any] | None = None,
    token: str | None = None,
    method: str | None = None,
) -> dict[str, Any]:
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    body = json.dumps(data).encode() if data is not None else None
    req = urllib.request.Request(url, data=body, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as e:
        detail = e.read().decode(errors="replace")
        raise SystemExit(f"HTTP {e.code} for {url}: {detail}") from None


def api_get(path: str, params: dict[str, Any] | None = None, token: str | None = None) -> dict[str, Any]:
    url = API_URL + path
    if params:
        url += "?" + urllib.parse.urlencode({k: v for k, v in params.items() if v is not None})
    return http_json(url, token=token)


def load_creds() -> dict[str, Any] | None:
    if CRED_PATH.exists():
        return json.loads(CRED_PATH.read_text())
    return None


def save_creds(creds: dict[str, Any]) -> None:
    CRED_PATH.parent.mkdir(parents=True, exist_ok=True)
    CRED_PATH.write_text(json.dumps(creds, indent=2) + "\n")
    CRED_PATH.chmod(0o600)


def id_token() -> str:
    """Return a fresh Firebase ID token, refreshing via the stored refresh token."""
    creds = load_creds()
    if not creds:
        raise SystemExit("Not logged in; run `scripts/yosupo.py login` first.")
    if creds.get("id_token") and creds.get("expires_at", 0) > time.time() + 60:
        return creds["id_token"]
    resp = http_json(
        f"https://securetoken.googleapis.com/v1/token?key={FIREBASE_API_KEY}",
        data={"grant_type": "refresh_token", "refresh_token": creds["refresh_token"]},
    )
    creds["id_token"] = resp["id_token"]
    creds["refresh_token"] = resp["refresh_token"]
    creds["expires_at"] = time.time() + int(resp["expires_in"])
    save_creds(creds)
    return creds["id_token"]


def default_user() -> str | None:
    creds = load_creds()
    return creds.get("user_name") if creds else None


def local_problems() -> dict[str, list[pathlib.Path]]:
    """Map problem slug -> verify/ files, from PROBLEM comments."""
    out: dict[str, list[pathlib.Path]] = {}
    for path in sorted((ROOT / "verify").glob("*.test.cpp")):
        m = PROBLEM_RE.search(path.read_text())
        if m:
            out.setdefault(m.group(1), []).append(path.relative_to(ROOT))
    return out


def normalize_source(source: str) -> list[str]:
    """Source lines modulo #line directives and trailing whitespace."""
    lines = [l.rstrip() for l in source.splitlines()]
    return [l for l in lines if l and not l.startswith("#line ")]


def ac_sources(user: str, slug: str, limit: int) -> list[tuple[int, list[str]]]:
    """(id, normalized source) of the user's most recent AC submissions."""
    resp = api_get(
        "/submissions",
        params={"user": user, "problem": slug, "status": "AC", "limit": limit},
    )
    ids = [s["id"] for s in resp["submissions"]]
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        infos = list(ex.map(lambda i: api_get(f"/submissions/{i}"), ids))
    return [(i, normalize_source(info["source"])) for i, info in zip(ids, infos)]


def source_freshness(
    path: pathlib.Path, user: str, slug: str, limit: int
) -> str:
    """Compare the local bundle against recent AC submissions for slug."""
    local = normalize_source(bundle_mod.bundle(ROOT / path).decode())
    subs = ac_sources(user, slug, limit)
    if not subs:
        return "no AC submissions"
    for sub_id, remote in subs:
        if remote == local:
            return f"up-to-date (matches #{sub_id})"
    sub_id, remote = subs[0]
    added = removed = 0
    for d in difflib.unified_diff(remote, local, n=0):
        added += d.startswith("+") and not d.startswith("+++")
        removed += d.startswith("-") and not d.startswith("---")
    return f"differs from #{sub_id} (+{added}/-{removed} lines)"


def cmd_login(args: argparse.Namespace) -> None:
    email = args.email or input("Email: ")
    password = getpass.getpass("Password: ")
    resp = http_json(
        f"https://identitytoolkit.googleapis.com/v1/accounts:signInWithPassword?key={FIREBASE_API_KEY}",
        data={"email": email, "password": password, "returnSecureToken": True},
    )
    creds = {
        "refresh_token": resp["refreshToken"],
        "id_token": resp["idToken"],
        "expires_at": time.time() + int(resp["expiresIn"]),
    }
    save_creds(creds)
    info = api_get("/auth/current_user", token=creds["id_token"])
    user = (info.get("user") or {}).get("name")
    creds["user_name"] = user
    save_creds(creds)
    print(f"Logged in as {user!r}; credentials saved to {CRED_PATH}")


def cmd_status(args: argparse.Namespace) -> None:
    user = args.user or default_user()
    if not user:
        raise SystemExit("Pass --user or run `scripts/yosupo.py login` first.")
    solved: dict[str, str] = api_get(f"/users/{user}/statistics")["solved_map"]
    local = local_problems()

    nfiles = sum(len(v) for v in local.values())
    print(f"== verify/ files ({nfiles}) vs. judge status for {user} ==")
    for slug, paths in local.items():
        for path in paths:
            line = f"{solved.get(slug, '--'):>9}  {slug}  ({path})"
            if args.compare and slug in solved:
                line += f"  [{source_freshness(path, user, slug, args.depth)}]"
            print(line)

    extra = sorted(set(solved) - set(local))
    if extra and not args.local_only:
        print(f"\n== solved on judge but no verify/ file ({len(extra)}) ==")
        for slug in extra:
            print(f"{solved[slug]:>9}  {slug}")


def cmd_table(args: argparse.Namespace) -> None:
    user = args.user or default_user()
    if not user:
        raise SystemExit("Pass --user or run `scripts/yosupo.py login` first.")
    solved: dict[str, str] = api_get(f"/users/{user}/statistics")["solved_map"]
    categories = api_get("/categories")["categories"]
    local = local_problems()

    for cat in categories:
        rows = []
        for slug in cat["problems"]:
            paths = local.get(slug, [])
            judge = solved.get(slug, "")
            if args.todo and paths and judge:
                continue
            files = ", ".join(str(p) for p in paths)
            rows.append((slug, "x" if paths else "", judge, files))
        if not rows:
            continue
        print(f"\n## {cat['title']}")
        width = max(len(r[0]) for r in rows)
        print(f"{'problem':<{width}}  local  judge      file")
        for slug, has_local, judge, files in rows:
            print(f"{slug:<{width}}  {has_local:<5}  {judge:<9}  {files}")


def cmd_diff(args: argparse.Namespace) -> None:
    path = args.path.resolve()
    m = PROBLEM_RE.search(path.read_text())
    if not m:
        raise SystemExit("No PROBLEM comment found in file.")
    slug = m.group(1)
    local = normalize_source(bundle_mod.bundle(path).decode())
    if args.submission:
        sub_id = args.submission
        remote = normalize_source(api_get(f"/submissions/{sub_id}")["source"])
    else:
        user = args.user or default_user()
        if not user:
            raise SystemExit("Pass --user or run `scripts/yosupo.py login` first.")
        subs = ac_sources(user, slug, 1)
        if not subs:
            raise SystemExit(f"No AC submissions by {user} for {slug}.")
        sub_id, remote = subs[0]
    for line in difflib.unified_diff(
        remote, local, f"submission #{sub_id}", str(args.path), lineterm=""
    ):
        print(line)


def cmd_submissions(args: argparse.Namespace) -> None:
    user = args.user or default_user()
    if not user:
        raise SystemExit("Pass --user or run `scripts/yosupo.py login` first.")
    resp = api_get(
        "/submissions",
        params={"user": user, "problem": args.problem, "limit": args.limit, "skip": args.skip},
    )
    print(f"count: {resp['count']}")
    for s in resp["submissions"]:
        print(
            f"#{s.get('id', 0):>7}  {s.get('status', '?'):>9}  {s.get('time', 0):7.3f}s"
            f"  {s.get('memory', 0) / 2**20:8.1f}MiB  {s.get('submission_time', '-'):>27}"
            f"  {s.get('problem_name', '?')}"
        )


PENDING_RE = re.compile(r"^(WJ|\d+/\d+)$")


def cmd_submit(args: argparse.Namespace) -> None:
    path = args.path.resolve()
    problem = args.problem
    if not problem:
        m = PROBLEM_RE.search(path.read_text())
        if not m:
            raise SystemExit("No PROBLEM comment found; pass --problem <slug>.")
        problem = m.group(1)
    source = bundle_mod.bundle(path).decode()
    resp = http_json(
        f"{API_URL}/submit",
        data={"problem": problem, "source": source, "lang": args.lang},
        token=id_token(),
    )
    sub_id = resp["id"]
    print(f"Submitted: https://judge.yosupo.jp/submission/{sub_id}")
    if not args.wait:
        return
    while True:
        time.sleep(2)
        info = api_get(f"/submissions/{sub_id}")["overview"]
        status = info.get("status", "WJ")
        print(f"  {status}")
        if not PENDING_RE.match(status):
            print(
                f"{status}: {info.get('time', 0):.3f}s,"
                f" {info.get('memory', 0) / 2**20:.1f}MiB"
            )
            sys.exit(0 if status == "AC" else 1)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("login", help="log in and store a refresh token")
    p.add_argument("--email")
    p.set_defaults(func=cmd_login)

    p = sub.add_parser("status", help="compare verify/ files against judge AC status")
    p.add_argument("--user")
    p.add_argument("--local-only", action="store_true", help="skip listing solved-but-missing problems")
    p.add_argument(
        "--compare",
        action="store_true",
        help="also compare each local bundle against your recent AC submission sources",
    )
    p.add_argument("--depth", type=int, default=5, help="how many recent ACs to compare against (default 5)")
    p.set_defaults(func=cmd_status)

    p = sub.add_parser("table", help="all judge problems by category, with local/judge solved status")
    p.add_argument("--user")
    p.add_argument("--todo", action="store_true", help="only show problems missing a local file or a judge AC")
    p.set_defaults(func=cmd_table)

    p = sub.add_parser("diff", help="diff a local bundle against your latest AC submission")
    p.add_argument("path", type=pathlib.Path)
    p.add_argument("--user")
    p.add_argument("--submission", type=int, help="diff against this submission id instead")
    p.set_defaults(func=cmd_diff)

    p = sub.add_parser("submissions", help="list your submissions")
    p.add_argument("--user")
    p.add_argument("--problem")
    p.add_argument("--limit", type=int, default=20)
    p.add_argument("--skip", type=int, default=0)
    p.set_defaults(func=cmd_submissions)

    p = sub.add_parser("submit", help="bundle and submit a file")
    p.add_argument("path", type=pathlib.Path)
    p.add_argument("--problem", help="problem slug (default: from PROBLEM comment)")
    p.add_argument("--lang", default="cpp")
    p.add_argument("--wait", action="store_true", help="poll until judging finishes")
    p.set_defaults(func=cmd_submit)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
