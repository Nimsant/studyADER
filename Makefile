#https://tech.davis-hansson.com/p/make/
SHELL := bash
GCC :=g++
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

MODELS = Burgers.cpp nonConservativeBurgers.cpp Advection.cpp seismic.cpp eqs4testing.cpp gas.cpp

FLUX = flux.cpp

OPTIONS = -O3 -std=c++2a -I fmt/include  -larmadillo
DEFINES = 
					
					
					

ALL: main

GENERATED_FILES = gen_derivatives_l.cpp gen_polynom_at_one.cpp gen_polynom_at_zero.cpp gen_polynomials_l.cpp gen_roots.cpp gen_weights.cpp
gen: $(GENERATED_FILES)
.PHONY: gen

$(GENERATED_FILES): PMpolynoms.py
	python3 $<

main.o: main.cpp $(GENERATED_FILES) $(MODELS) $(FLUX)
	$(GCC) $< $(OPTIONS) $(DEFINES) -c -o $@

main: main.o 
	$(GCC) $< $(OPTIONS) $(DEFINES) -o $@ 

run: main 
	rm -f *.dat
	./$<
	#./$< > run.log
.PHONY: run

show: run
	python3 plot_output.py
.PHONY: show

theory.pdf: theory/a.tex
	cd $(dir $<)
	pdflatex $(notdir $<)
	mv $(notdir $(<:tex=pdf)) ../$@

clean:
	rm -f main.o main
.PHONY: clean

cleanall:
	rm -f main.o main $(GENERATED_FILES) *.dat
.PHONY: cleanall
