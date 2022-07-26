CXX :=g++ -std=c++11 
CXXFLAGS:= -g -O2 -pg  #-fopenmp -W -Wall  
LIBFLAGS:= -lgsl -lgslcblas -lm /software/gsl/2.7.1/lib/libgsl.a            #-lfftw3_omp  -lfftw3
INCFLAG:= -I ./include  -I/software/gsl/2.7.1/include/  -I/usr/include -L/usr/local/lib
OBJDIR = obj
VPATH  = %.cpp src
source = $(wildcard *.cpp) $(notdir $(wildcard src/*.cpp))  
objs = $(source:%.cpp=$(OBJDIR)/%.o)

run: $(objs)  
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)
	@echo Make done
	
$(objs): $(OBJDIR)/%.o : %.cpp
	@mkdir -p $(OBJDIR)  
	$(CXX) $(CXXFLAGS) $(INCFLAG) -w -c $< -o $@ $(LIBFLAGS)

.PHONY: clean
clean:
	-rm -f obj/*.o

.PHONY: re
re:
	make clean
	make
