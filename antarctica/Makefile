# Antarctica timing benchmark: resolution vs. core-count matrix.
#
# Run from antarctica/:
#   make timing
#
# Stages (each idempotent — skips already-built artifacts):
#   meshes        build 10 meshes (5 LC x {10,20}*LC, buffer=20000 m)
#   inversion     single LC=2500 / LC_coarse=25000 inversion (Slurm or mpiexec)
#   redistribute  rewrite inversion checkpoint on 1 core
#   transient     30 short runs (10 mesh combos x 16/32/64 cores)
#   matrix        write TIMING_MATRIX.md from results/timing/*.json

LCS := 500 1000 2000 2500 5000
RATIOS := 20 10
CORES := 16 32 64
BUFFER := 20000
INV_LC := 2500
INV_LC_COARSE := 25000
INV_CHK := mesh/inversion_icepack2_$(INV_LC).h5

.PHONY: timing meshes inversion redistribute transient matrix

timing: matrix

meshes:
	@for lc in $(LCS); do for r in $(RATIOS); do \
	  lcc=$$((lc * r)); fn=mesh/antarctica_$${lcc}_$${lc}_buffered$(BUFFER).msh; \
	  [ -f $$fn ] && { echo "skip $$fn"; continue; }; \
	  python scripts/mesh_antarctica.py --lc $$lc --lc-coarse $$lcc --buffer-m $(BUFFER); \
	done; done

inversion: meshes
	@[ -f $(INV_CHK) ] && echo "skip inversion" || \
	if command -v sbatch >/dev/null 2>&1; then sbatch --wait scripts/batch_runners/timing_inversion.script; \
	else ISMIP7_LC=$(INV_LC) ISMIP7_LC_COARSE=$(INV_LC_COARSE) ISMIP7_BUFFER_M=$(BUFFER) \
	  mpiexec -n 12 python scripts/inversion_icepack2.py; fi

redistribute: inversion
	python scripts/redistribute_checkpoint.py --lc $(INV_LC)

transient: redistribute
	@mkdir -p results/timing
	@for lc in $(LCS); do for r in $(RATIOS); do lcc=$$((lc * r)); for nc in $(CORES); do \
	  out=results/timing/timing_$${lc}_$${lcc}_$${nc}.json; \
	  [ -f $$out ] && { echo "skip $$out"; continue; }; \
	  if command -v sbatch >/dev/null 2>&1; then \
	    sbatch --wait --job-name=timing_$${lc}_$${lcc}_$${nc} --ntasks-per-node=$$nc \
	      --export=ALL,ISMIP7_LC=$$lc,ISMIP7_LC_COARSE=$$lcc,ISMIP7_BUFFER_M=$(BUFFER),ISMIP7_INVERSION=$(CURDIR)/$(INV_CHK),ISMIP7_T_END=2020,ISMIP7_DT=1.0 \
	      scripts/batch_runners/timing_transient.script; \
	  else \
	    ISMIP7_LC=$$lc ISMIP7_LC_COARSE=$$lcc ISMIP7_BUFFER_M=$(BUFFER) ISMIP7_INVERSION=$(CURDIR)/$(INV_CHK) \
	      ISMIP7_T_END=2020 ISMIP7_DT=1.0 mpiexec -n $$nc python scripts/run_timing.py; \
	  fi; \
	done; done; done

matrix: transient
	python scripts/build_timing_matrix.py
