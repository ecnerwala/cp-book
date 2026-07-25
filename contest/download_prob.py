#!/usr/bin/env python3
"""Download and setup problems from Competitive Companion"""

import argparse
import http.server
import json
from pathlib import Path
import re
from typing import Any

import make_prob as make_prob_mod

LISTEN_PORT = 10046

ProblemData = dict[str, Any]

# Returns unmarshalled or None
def listen_once(*, timeout: float | None = None) -> ProblemData | None:
    json_data: ProblemData | None = None

    class CompetitiveCompanionHandler(http.server.BaseHTTPRequestHandler):
        def do_POST(self) -> None:
            nonlocal json_data
            json_data = json.loads(self.rfile.read(int(self.headers['Content-Length'])))
            self.send_response(200)
            self.end_headers()

    with http.server.HTTPServer(('127.0.0.1', LISTEN_PORT), CompetitiveCompanionHandler) as server:
        server.timeout = timeout
        server.handle_request()

    if json_data is not None:
        print(f"Got data {json.dumps(json_data)}")
    else:
        print("Got no data")
    return json_data

def listen_many(
    *,
    num_items: int | None = None,
    num_batches: int | None = None,
    timeout: float | None = None,
) -> list[ProblemData]:
    if num_items is not None:
        res = []
        for _ in range(num_items):
            cur = listen_once(timeout=None)
            assert cur is not None
            res.append(cur)
        return res

    if num_batches is not None:
        res = []

        batches: dict[str, list[int]] = {}
        while len(batches) < num_batches or any(need for need, tot in batches.values()):
            print(f"Waiting for {num_batches} batches:", batches)
            cur = listen_once(timeout=None)
            assert cur is not None
            res.append(cur)

            cur_batch = cur['batch']
            batch_id = cur_batch['id']
            batch_cnt = cur_batch['size']
            if batch_id not in batches:
                batches[batch_id] = [batch_cnt, batch_cnt]
            assert batches[batch_id][0] > 0
            batches[batch_id][0] -= 1

        return res

    first = listen_once(timeout=None)
    assert first is not None
    res = [first]
    while True:
        cnd = listen_once(timeout=timeout)
        if cnd is None:
            break
        res.append(cnd)
    return res

NAME_PATTERN = re.compile(r'^(?:Problem )?([A-Z][0-9]*)\b')

def get_prob_name(data: ProblemData) -> str:
    if 'USACO' in data['group']:
        if 'fileName' in data['input']:
            names: list[str] = [
                str(data['input']['fileName']).removesuffix('.in'),
                str(data['output']['fileName']).removesuffix('.out'),
            ]
            if len(set(names)) == 1:
                return names[0]

    if 'url' in data and data['url'].startswith('https://www.codechef.com'):
        return str(data['url']).rstrip('/').rsplit('/')[-1]

    patternMatch = NAME_PATTERN.search(data['name'])
    if patternMatch is not None:
        return patternMatch.group(1)

    print(f"For data: {json.dumps(data, indent=2)}")
    return input("What name to give? ")

def save_samples(data: ProblemData, prob_dir: Path) -> None:
    with open(prob_dir / 'problem.json', 'w') as f:
        json.dump(data, f)

    for i, t in enumerate(data['tests'], start=1):
        with open(prob_dir / f'sample{i}.in', 'w') as f:
            f.write(t['input'])
        with open(prob_dir / f'sample{i}.out', 'w') as f:
            f.write(t['output'])

# Providing name = '.'
def make_prob(data: ProblemData, name: str | None = None) -> None:
    if name is None:
        name = get_prob_name(data)

    prob_dir = Path('.')/name

    if name == '.':
        print("Using current directory...")
        pass
    elif prob_dir.exists() and prob_dir.is_dir():
        # Skip making it
        print(f"Already created problem {name}...")
    else:
        print(f"Creating problem {name}...")

        if not make_prob_mod.make_prob(name):
            return

    print("Saving samples...")
    save_samples(data, prob_dir)

    print()

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--echo', action='store_true', help='just echo received responses and exit',
    )
    parser.add_argument(
        '--dryrun', action='store_true', help="don't actually create any problems",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        'names', nargs='*', metavar='name', default=[],
        help='problem names (one problem is downloaded per name)',
    )
    group.add_argument('-n', '--number', type=int, help='number of problems')
    group.add_argument(
        '-b', '--batches', type=int, help='number of batches (default: 1 batch)',
    )
    group.add_argument(
        '-t', '--timeout', type=float,
        help='listen until this many seconds pass without a problem',
    )
    args = parser.parse_args()

    if args.echo:
        while True:
            print(listen_once())
    else:
        def run_make_prob(*fn_args: Any, **fn_kwargs: Any) -> None:
            if args.dryrun:
                print(f"make_prob(*args={fn_args}, **kwargs={fn_kwargs})")
                return
            make_prob(*fn_args, **fn_kwargs)

        if args.names:
            datas = listen_many(num_items=len(args.names))
            for data, name in zip(datas, args.names):
                run_make_prob(data, name)
        elif args.number is not None:
            datas = listen_many(num_items=args.number)
            for data in datas:
                run_make_prob(data)
        elif args.batches is not None:
            datas = listen_many(num_batches=args.batches)
            for data in datas:
                run_make_prob(data)
        elif args.timeout is not None:
            datas = listen_many(timeout=args.timeout)
            for data in datas:
                run_make_prob(data)
        else:
            datas = listen_many(num_batches=1)
            for data in datas:
                run_make_prob(data)

if __name__ == '__main__':
    main()
