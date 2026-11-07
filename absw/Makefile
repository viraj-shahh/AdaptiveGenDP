# Adaptive banded Smith-Waterman: one entry point for everything in absw/.
#
#   make absw       synthesize the Smith-Waterman array in ASAP7 with Cadence Genus
#                   (fetches the ASAP7 libraries on first run) -> flow/genus/<run>/reports
#   make sim        RTL vs golden model (Verilator by default; SIM=vcs on CAEN)
#   make golden     C++ golden model (golden/): selftest + band-policy comparison
#   make baseline   build GenDP simulator + BWA ksw-test, check them against each other
#   make orfs       full RTL-to-GDSII with OpenROAD-flow-scripts (needs ORFS=<path>)
#   make libs       re-fetch the ASAP7 standard-cell libraries only
#   make clean
#
# Knobs:  W=16  MAXLEN=128  CLK_PS=1000  CORNER=TT  VT=RVT   (absw)
#         SIM=verilator|vcs  SIM_WS="8 16 32"                    (sim)

W       ?= 16
MAXLEN  ?= 128
CLK_PS  ?= 1000
SIM     ?= verilator
CORNER  ?= TT
VT      ?= RVT
SIM_WS  ?= 8 16 32

ASAP7_SRC := flow/.asap7_src
ASAP7_LIB := flow/.asap7_lib
RUN       := W$(W)_L$(MAXLEN)_$(CLK_PS)ps_$(VT)_$(CORNER)
GENUS_DIR := flow/genus/$(RUN)

# ASAP7 7.5-track cells (file names as shipped in ORFS flow/platforms/asap7)
LIBS = asap7sc7p5t_AO_$(VT)_$(CORNER)_nldm_211120.lib \
       asap7sc7p5t_INVBUF_$(VT)_$(CORNER)_nldm_220122.lib \
       asap7sc7p5t_OA_$(VT)_$(CORNER)_nldm_211120.lib \
       asap7sc7p5t_SIMPLE_$(VT)_$(CORNER)_nldm_211120.lib \
       asap7sc7p5t_SEQ_$(VT)_$(CORNER)_nldm_220123.lib

.PHONY: absw sim golden baseline orfs libs clean

GOLDEN := build/golden
$(GOLDEN): golden/golden.cpp golden/absw_golden.hpp
	mkdir -p build
	$(CXX) -std=c++17 -O2 -Wall -Wextra -o $@ golden/golden.cpp

golden: $(GOLDEN)
	mkdir -p data
	$(GOLDEN) selftest
	$(GOLDEN) pairs -n 60 --qlen 1000 --tlen-extra 200 --sub-rate 0.04 --indel-rate 0.06 \
	    --max-indel 2 --del-frac 0.65 --h0 30 --max-len 0 --seed 4 --tsv data/long_drift.tsv
	$(GOLDEN) compare data/long_drift.tsv -W 8 16 32

baseline:
	bash run_all.sh

sim:
	SIM=$(SIM) WS="$(SIM_WS)" bash rtl/run_tests.sh

# Sparse checkout: only the ASAP7 liberty and LEF folders of ORFS (~100 MB), no build.
libs:
	rm -f $(ASAP7_LIB)/.done
	$(MAKE) $(ASAP7_LIB)/.done
$(ASAP7_LIB)/.done:
	rm -rf $(ASAP7_SRC) && mkdir -p $(ASAP7_LIB)
	git clone --depth 1 --filter=blob:none --sparse \
	    https://github.com/The-OpenROAD-Project/OpenROAD-flow-scripts.git $(ASAP7_SRC)
	cd $(ASAP7_SRC) && git sparse-checkout set flow/platforms/asap7/lib/NLDM flow/platforms/asap7/lef
	for f in $(ASAP7_SRC)/flow/platforms/asap7/lib/NLDM/asap7sc7p5t_*_nldm_*.lib*; do \
	    case $$f in *FAKE*) continue;; esac; \
	    case $$f in *.gz) gunzip -c $$f > $(ASAP7_LIB)/$$(basename $$f .gz);; *) cp $$f $(ASAP7_LIB)/;; esac; \
	done
	cp $(ASAP7_SRC)/flow/platforms/asap7/lef/asap7_tech_1x_201209.lef \
	   $(ASAP7_SRC)/flow/platforms/asap7/lef/asap7sc7p5t_28_*_1x_220121a.lef $(ASAP7_LIB)/
	touch $@

absw: $(ASAP7_LIB)/.done
	mkdir -p $(GENUS_DIR)
	cd $(GENUS_DIR) && \
	ABSW_RTL="$(CURDIR)/rtl" ABSW_LIBDIR="$(CURDIR)/$(ASAP7_LIB)" ABSW_LIBS="$(LIBS)" \
	ABSW_SDC="$(CURDIR)/flow/asap7/constraint.sdc" ABSW_W=$(W) ABSW_MAXLEN=$(MAXLEN) ABSW_CLK_PS=$(CLK_PS) \
	    genus -batch -no_gui -files "$(CURDIR)/flow/genus/genus.tcl" -log genus
	@echo "reports: $(GENUS_DIR)/reports"
	@grep -E "Total|absw" $(GENUS_DIR)/reports/area.rpt | head -5 || true
	@grep -m1 -E "Slack" $(GENUS_DIR)/reports/timing.rpt || true

orfs:
	@test -n "$(ORFS)" || { echo "set ORFS=<path to OpenROAD-flow-scripts>"; exit 1; }
	$(MAKE) -C $(ORFS)/flow DESIGN_CONFIG="$(CURDIR)/flow/asap7/config.mk"

clean:
	rm -rf build rtl/build flow/genus/W* data
