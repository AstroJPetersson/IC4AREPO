CXX         = /opt/ebsofts/GCCcore/11.2.0/bin/g++
HEADERS     = ./potential.h ./common.h
COMMON_OBJ  = ./potential.o
CXXFLAGS    = -O3


all : potential_check

potential.o : ./potential.cc $(HEADERS) 

potential_check : $(HEADERS) $(COMMON_OBJ) potential_check.o
	$(CXX) $(CXXFLAGS) potential_check.o $(COMMON_OBJ) -lfftw3 -lgsl -lgslcblas -lm -o potential_check

clean:
	@rm -f *.o potential_check ../potential.o


