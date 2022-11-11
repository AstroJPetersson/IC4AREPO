HEADERS     = ./potential.h ./common.h 
COMMON_OBJ  = ./potential.o 
CC          = mpicc
CXXFLAGS    = -O3 -I/home/hd/hd_hd/hd_uw148/lib/eigen -L/opt/bwhpc/common/numlib/gsl/2.1-intel-15.0/lib -I/opt/bwhpc/common/compiler/intel/compxe.2015.3.187/composer_xe_2015.3.187/debugger/gdb/intel64/include -L/home/hd/hd_hd/hd_uw148/lib/lib -I/home/hd/hd_hd/hd_uw148/lib/include -L/opt/bwhpc/common/numlib/gsl/2.1-intel-15.0/lib


all : potential_check

potential.o : ./potential.cc $(HEADERS) 

potential_check : $(HEADERS) $(COMMON_OBJ) potential_check.o
	$(CC) $(CXXFLAGS) potential_check.o $(COMMON_OBJ) -lfftw3 -lgsl -lgslcblas -lm -o potential_check

clean:
	@rm -f *.o potential_check ../potential.o

