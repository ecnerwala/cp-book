# +--------------------+
# |                    |
# |   GENERAL CONFIG   |
# |                    |
# +--------------------+

PROBLEM_NAME := __PROBLEM_NAME__
DEBUG := true
ALLOW_FAILS := false
LANG := cpp
RUN_TARGET := __run_default

# Local config overrides (e.g. from a layered .template); optional.
# Included before the config is used, so it can override LANG, DEBUG, flags.
-include config.mk

ifeq ($(LANG),cpp)
TARGET := $(PROBLEM_NAME)
EXECUTE := ./$(TARGET)
CLEAN_TARGETS := $(TARGET)
else ifeq ($(LANG),python)
TARGET := $(PROBLEM_NAME).py
EXECUTE := python3 ./$(TARGET)
CLEAN_TARGETS :=
else ifeq ($(TARGET),)
$(error "Unknown LANG '$(LANG)'; set TARGET, EXECUTE, and CLEAN_TARGETS (e.g. in config.mk)")
endif

CXX := g++
CXXFLAGS := -std=c++23 -O2 -g -Wall -Wextra -pedantic -Wshadow -Wformat=2 -Wfloat-equal -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align -Wno-unused-result -Wno-sign-conversion
CXXFLAGS += $(EXTRA_CXXFLAGS)
# No -D_FORTIFY_SOURCE: it rewrites libc calls to *_chk variants that bypass
# ASan's interceptors, silently weakening both checkers.
DEBUG_CXXFLAGS := -fsanitize=address -fsanitize=undefined -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize-recover=all -fstack-protector-all
DEBUG_CXXFLAGS += -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC

#export ASAN_OPTIONS=detect_leaks=0  # Disable leak checking entirely
export LSAN_OPTIONS=exitcode=0

PRECOMPILE_HEADERS := bits/stdc++.h
#PRECOMPILE_HEADERS := bits/extc++.h


# +-------------------+
# |                   |
# |   GENERAL RULES   |
# |                   |
# +-------------------+

all: $(TARGET)
.PHONY: all

clean:
	-rm -rf $(CLEAN_TARGETS)
.PHONY: clean

veryclean:
	-rm -rf $(CLEAN_TARGETS) *.res
.PHONY: veryclean

rebuild: clean all
.PHONY: rebuild

# +---------------------+
# |                     |
# |   C++ COMPILATION   |
# |                     |
# +---------------------+

ifeq ($(DEBUG),true)
CXXFLAGS += $(DEBUG_CXXFLAGS)
endif


PCH_ROOT := $(HOME)/.precompiled_headers

NULL :=
SPACE := $(NULL) $(NULL)
PCH_FILENAME := $(subst $(SPACE),_,$(shell sha256sum $$(which $(CXX)) | cut -d ' ' -f1)-$(shell echo '$(CXXFLAGS)' | sha256sum | cut -d ' ' -f1)).gch
PCH_FILES := $(patsubst %,$(PCH_ROOT)/%.gch/$(PCH_FILENAME),$(PRECOMPILE_HEADERS))

$(PCH_ROOT)/%.gch/$(PCH_FILENAME):
	rm -f $@
	mkdir -p $(dir $@)
	$(LINK.cpp) -x c++-header "$$(echo '#include<$*>' | $(LINK.cpp) -H -E -x c++ - 2>&1 >/dev/null | head -1 | cut -d ' ' -f2)" -o $@
.PRECIOUS: $(PCH_ROOT)/%.gch/$(PCH_FILENAME)


%: %.cpp # Cancel the builtin rule

%: %.cpp $(PCH_FILES)
	mkdir -p archive
	cp $< archive/$*-$$(date +%Y-%m-%d-%H%M%S).cpp
	$(LINK.cpp) -isystem $(PCH_ROOT) $< $(LOADLIBES) $(LDLIBS) -o $@
.PRECIOUS: %


# +-----------------------+
# |                       |
# |   RUNNING / TESTING   |
# |                       |
# +-----------------------+

export TIME=\n  real\t%es\n  user\t%Us\n  sys\t%Ss\n  mem\t%MKB

run: $(RUN_TARGET)
.PHONY: run

__run_default: $(TARGET)
	\time $(EXECUTE)
ifeq ($(DEBUG),true)
	@echo "Built with DEBUG flags enabled, code may be slower than normal"
endif
.PHONY: __run_default

%.res: $(TARGET) %.in
	\time $(EXECUTE) < $*.in > $*.res
ifeq ($(DEBUG),true)
	@echo "Built with DEBUG flags enabled, code may be slower than normal"
endif
.PRECIOUS: %.res

%.out: % # Cancel the builtin rule

__test_%: %.res %.out
ifeq ($(ALLOW_FAILS),true)
	-diff $*.res $*.out
else
	diff $*.res $*.out
endif
.PHONY: __test_%

CASES := $(sort $(basename $(wildcard *.in)))
TESTS := $(sort $(basename $(wildcard *.out)))

runs: $(patsubst %,%.res,$(CASES))
.PHONY: runs

solve: runs
.PHONY: solve

test: $(patsubst %,__test_%,$(TESTS))
.PHONY: test


# Local rule additions/overrides (e.g. from a layered .template); optional.
# Included last so extra rules don't change the default goal.
-include local.mk
