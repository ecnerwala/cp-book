# Contest tooling

Standalone tooling for setting up and testing contest problems. Nothing here
depends on the library — contest programs are meant to be self-contained
(vendor any needed book code into the file, e.g. with `scripts/bundle.py` or
copy-paste).

## Setup

Copy the base template to the root of your contests tree (copy-once; no live
dependency on this repo):

```sh
cp -r contest/template ~/contests/.template
```

Put `download_prob.py` and `make_prob.py` somewhere on your `PATH` (or invoke
them by path). Both are stdlib-only Python 3.

## Usage

```sh
# Create problem directories from the nearest .template(s)
make_prob.py A B C

# Listen for Competitive Companion (default: one batch, e.g. "parse all")
download_prob.py            # names inferred from the problems
download_prob.py A B C      # explicit names
download_prob.py -n 3       # exactly 3 problems
download_prob.py -t 5       # keep listening until 5s of silence

# Inside a problem directory
make            # build (with sanitizers; set DEBUG := false to disable)
make run        # build and run interactively
make runs       # run against all *.in
make test       # diff *.res against *.out
```

## Layered templates

`make_prob.py` searches upward from the current directory for `.template`
directories. A `.template` containing a file named `__PARENT__` also inherits
from `.template`s further up; templates are applied top-down, so deeper ones
override shallower ones file-by-file. Use this for per-site or per-contest
customization:

```
~/contests/.template/            # base (copied from this repo)
~/contests/codeforces/.template/ # site-specific overrides; contains __PARENT__
```

For sub-file customization of the Makefile, deeper `.template`s can ship two
optional hook files that the base Makefile `-include`s: `config.mk` (included
first; override `LANG`, `DEBUG`, flags — e.g. `EXTRA_CXXFLAGS := -I ...` is
appended to `CXXFLAGS`, and a later `-std=` there overrides the default) and
`local.mk` (included last; add or
override rules without changing the default goal). For custom languages, set
`TARGET`, `EXECUTE`, and `CLEAN_TARGETS` in `config.mk` and provide the build
rule in `local.mk`; to change what `make run` does, point `RUN_TARGET` at your
own target instead of redefining `run`. If a variant diverges
heavily, just ship a full replacement `Makefile` instead.

`__PROBLEM_NAME__` is substituted in file names and contents; nothing else
is touched (`$` has no special meaning), so Makefiles template safely.

If the instantiated directory contains an executable `setup`, it is run once
from inside it (the default one runs `bear -- make all` to generate
`compile_commands.json` for clangd).
