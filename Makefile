CXX := mpic++
CXXFLAGS:=  -g -fopenmp -O2 -std=c++11 -pg
LIBFLAGS:=   -lgsl -lgslcblas -lm            #-lfftw3_omp  -lfftw3
INCFLAG:= -I ./include  -I/usr/local/include/gsl 
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
