
CXX :=nvcc  -std=c++11  -w
INCFLAG:= -I /usr/local/cuda-10.1/targets/x86_64-linux/include  -I/usr/include -I ./include -I/software/gsl/2.7.1/include   
LIBFLAGS:= -lgsl -lgslcblas -lm  -lcuda -lcufft  -lm -lfftw3 -lfftw -lcufftw  #/software/gsl/2.7.1/lib/libgsl.a   #-lfftw3_omp  -lfftw3
#CXXFLAGS:= -arch=sm_30 
#CXXFLAGS:=  -arch=sm_30 -std=gnu++11
source = src/*.cpp src/*.cu 


run: $(source)  
	$(CXX) $(CXXFLAGS) $(INCFLAG) -o $@ $^ $(LIBFLAGS)
	@echo Make done


.PHONY: clean
clean:
	-rm -f obj/*.o *.o run

.PHONY: re
re:
	make clean
	make


