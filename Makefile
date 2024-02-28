CXX:= g++ -std=c++11  -w
INCFLAG:= -I  -I/usr/include -I ./include -I/software/gsl/2.7.1/include   #/usr/local/cuda-10.1/targets/x86_64-linux/include 
LIBFLAGS:= -lgsl -lgslcblas -lm  -lm -lfftw3 -lfftw  /software/gsl/2.7.1/lib/libgsl.a   #  -lcuda -lcufft -lcufftw  -lfftw3_omp  -lfftw3
#CXXFLAGS:= -arch=sm_30 
#CXXFLAGS:=  -arch=sm_30 -std=gnu++11
source = src/*.cpp  #src/*.cu 

OBJDIR = obj
objs = $(source:%.cpp=$(OBJDIR)/%.o)


run: $(source)  
	$(CXX) $(CXXFLAGS) $(INCFLAG) -o $@ $^ $(LIBFLAGS)
	@echo Make done


$(objs): $(OBJDIR)/%.o : %.cpp %.cu
	@mkdir -p $(OBJDIR)  
	$(CXX) $(CXXFLAGS) $(INCFLAG) -w -c $< -o $@ $(LIBFLAGS)


.PHONY: clean
clean:
	-rm -f obj/*.o *.o run

.PHONY: re
re:
	make clean
	make


