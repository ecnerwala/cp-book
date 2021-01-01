CXX = g++
CXXFLAGS = -Wall -Wextra -pedantic -std=c++14 -O2 -Wshadow -Wformat=2 -Wfloat-equal -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align
DEBUGFLAGS = -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fstack-protector -D_FORTIFY_SOURCE=2

CXXFLAGS += -I..
HEADERS := $(sort $(wildcard ../*.hpp))

DEBUG = false
ifeq ($(DEBUG),true)
	CXXFLAGS += $(DEBUGFLAGS)
endif

MAIN_TARGET := main

CASES := $(sort $(basename $(wildcard *.in)))
TESTS := $(sort $(basename $(wildcard *.out)))
TARGETS := $(CASES)

all: $(TARGET) $(TARGETS)

clean:
	-rm -rf $(CASES) *.res

%.out: %

%: %.cpp $(HEADERS)
	$(LINK.cpp) $< $(LOADLIBES) $(LDLIBS) -o $@

run: $(MAIN_TARGET)
	time ./$(MAIN_TARGET)

%.res: % %.in
	time ./$* < $*.in > $*.res
ifeq ($(DEBUG),true)
	@echo "Built with DEBUG flags enabled, code may be slower than normal"
endif

%.out: % # Cancel the builtin rule

__test_%: %.res %.out
	diff $*.res $*.out

runs: $(patsubst %,%.res,$(CASES))

solve: runs

test: $(patsubst %,__test_%,$(TESTS))

.PHONY: all clean run test __test_% runs solve

.PRECIOUS: %.res

.SUFFIXES:
