CC      = gcc
CFLAGS  = -O2 -std=c11 -Wall -Wextra
LDLIBS  = -lm
TYPST  ?= typst

FIGS = figures/a1_delta_lambda.svg figures/a1_tomo_rows.svg figures/a1_regr_rows.svg figures/a2_det.svg

all: report/report.pdf

lab1: src/lab1.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

data/results.txt: lab1
	mkdir -p data
	./lab1 data

figures/a1_delta_lambda.svg: data/results.txt plots/a1_delta_lambda.gp plots/common.gp
	mkdir -p figures
	gnuplot plots/a1_delta_lambda.gp

figures/a1_%_rows.svg: data/results.txt plots/a1_rows.gp plots/common.gp
	mkdir -p figures
	gnuplot -e "tag='$*'" plots/a1_rows.gp

figures/a2_det.svg: data/results.txt plots/a2_det.gp plots/common.gp
	mkdir -p figures
	gnuplot plots/a2_det.gp

report/report.pdf: report/report.typ $(FIGS)
	$(TYPST) compile --root . report/report.typ $@

clean:
	rm -rf lab1 data figures report/report.pdf

.PHONY: all clean
