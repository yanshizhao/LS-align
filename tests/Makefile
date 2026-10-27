# LS-align test suite -- build & run entry point
#
# Toolchain (verified): MSYS2/MINGW64 g++ 14.2.0 on Windows.
#
# IMPORTANT (environment): dynamic linking of libstdc++ is BROKEN on this box
# (collect2 exits 116). Always link with -static-libstdc++ -static-libgcc, or
# use -static. The README's own example uses -static for the same reason.
#
# Note: the LS-align binary always returns exit code 1 (LSalign.cpp:453) even on
# success, so CLI scripts must NOT rely on the exit status -- they parse stdout.

CXX      := g++
SRC      := ../src
INCLUDE  := -I$(SRC) -I.
STD      := -std=c++98            # shipping binary: fidelity to the original build
STD_UNIT := -std=c++17            # unit harness: required by leak_probe.h (aligned new)
STATIC   := -static-libstdc++ -static-libgcc
WARN     := -Wall -Wextra -Wno-unused-parameter
OPT      := -O2
CXXFLAGS := $(STD) $(INCLUDE) $(WARN) $(OPT) $(STATIC)
UNITFLAGS:= $(STD_UNIT) $(INCLUDE) $(WARN) $(OPT) $(STATIC)
BUILD    := build

# g++ needs a writable temp dir; under make TMP/TEMP/TMPDIR are empty, so g++
# would fall back to C:\Windows\ (permission denied).  Export one explicitly.
TMPDIR_X := $(CURDIR)/build/tmp
export TMP    := $(TMPDIR_X)
export TEMP   := $(TMPDIR_X)
export TMPDIR := $(TMPDIR_X)

# g++ needs the temp dir to exist; create it before any compile.
$(shell mkdir -p $(TMPDIR_X) $(BUILD) 2>/dev/null)

# Windows exec suffix (MSYS make maps this); empty on Linux.
EXE := .exe

.PHONY: all test unit cli asan golden check-golden clean help

all: unit cli

help:
	@echo "make test        - run the default CI suite (unit + cli; no ASan/dataset)"
	@echo "make unit        - build & run white-box unit tests"
	@echo "make cli         - build the binary & run black-box CLI tests"
	@echo "make asan        - build with ASan/UBSan & run leak/UB checks (xfail)"
	@echo "make golden      - (re)generate golden snapshots"
	@echo "make check-golden- diff current output against golden snapshots"
	@echo "make clean       - remove build artifacts"

test: unit cli

# --------------------------------------------------------------------------- #
# White-box unit tests: a single TU that #includes the headers (NOT LSalign.cpp,
# which owns main()).  fixtures/ and example paths are relative to tests/.
# --------------------------------------------------------------------------- #
$(BUILD)/unit$(EXE): unit/test_main.cpp mini_test.h leak_probe.h
	@mkdir -p $(BUILD)
	$(CXX) $(UNITFLAGS) unit/test_main.cpp -lm -o $@

unit: $(BUILD)/unit$(EXE)
	@echo "== running unit tests =="
	@./$(BUILD)/unit$(EXE)

# --------------------------------------------------------------------------- #
# Black-box CLI tests
# --------------------------------------------------------------------------- #
$(BUILD)/LSalign_under_test$(EXE): $(SRC)/LSalign.cpp
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) $(SRC)/LSalign.cpp -lm -o $@

cli: $(BUILD)/LSalign_under_test$(EXE)
	@echo "== running CLI tests =="
	@bash cli/run_cli.sh

# --------------------------------------------------------------------------- #
# Leak check: ASan/UBSan are NOT available on this MinGW toolchain (no libasan
# runtime; -fsanitize=address fails to link).  The dedicated leak checks live in
# unit/test_xfail.cpp using the portable leak_probe.h counter, so they already
# run as part of `make unit`.  This target exists as an explicit entry point and
# makes the environment limitation visible.
# --------------------------------------------------------------------------- #
leakcheck: unit
	@echo "== leak checks (via leak_probe.h; see the [XFAIL] leak_* cases above) =="
	@echo "   NOTE: ASan/UBSan unavailable on this MinGW toolchain -- see report."

# Kept as an alias so scripted invocations of 'make asan' still work.
asan: leakcheck

# --------------------------------------------------------------------------- #
# Golden regression (deterministic: rigid mode only on this fixture set)
# --------------------------------------------------------------------------- #
golden: $(BUILD)/LSalign_under_test$(EXE)
	./$(BUILD)/LSalign_under_test ../example/ace_ligands.mol2 ../example/xtal-lig.mol2 -rf 0 2>/dev/null | grep -vE '^(Using|Totally Using) ' > golden/ace_rigid.txt || true
	./$(BUILD)/LSalign_under_test fixtures/self.mol2 fixtures/copy.mol2 -rf 0 2>/dev/null | grep -vE '^(Using|Totally Using) ' > golden/self_copy.txt || true
	@echo "golden snapshots written (timing lines stripped)."

check-golden: $(BUILD)/LSalign_under_test$(EXE)
	@bash cli/check_golden.sh

clean:
	rm -rf $(BUILD)/unit$(EXE) $(BUILD)/unit_asan$(EXE) $(BUILD)/LSalign_under_test$(EXE) \
	       $(BUILD)/*.o $(BUILD)/sampled $(BUILD)/*.cpp
