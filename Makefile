GCC =g++

MX?=0
MT?=0
N?=8
NT?=1

GENERATED_FILES = gen_derivatives_l.cpp gen_polynom_at_one.cpp gen_polynom_at_zero.cpp gen_polynomials_l.cpp gen_roots.cpp gen_weights.cpp

OPTIONS = -O3 -std=c++2a -I fmt/include  -larmadillo
DEFINES = -DDEFINED_NT=$(NT) \
					-DDEFINED_N=$(N) \
					-DDEFINED_MT=$(MT) \
					-DDEFINED_MX=$(MX) 

ALL: main

gen: $(GENERATED_FILES)

$(GENERATED_FILES): PMpolynoms.py
	python3 $<

main.o: main.cpp $(GENERATED_FILES)
	$(GCC) $< $(OPTIONS) $(DEFINES) -c -o $@

main: main.o 
	$(GCC) $< $(OPTIONS) $(DEFINES) -o $@ 

run: main 
	rm -f *.dat
	./$<  > run.log

show: run
	python3 plot_output.py

clean:
	rm -f main.o main

cleanall:
	rm -f main.o main $(GENERATED_FILES) *.dat
