CXX = g++
CMPI = mpicxx
CXXFLAGS = -std=c++11 -Wall -O3
HDRS = io.h init.h oper.h
OBJS = main.o io.o init.o oper.o
LDLIBS = -llapack -lblas -lmpi

default:main

all:	main

%.o:	%.cpp $(HDRS)
	$(CMPI) $(CXXFLAGS) -o $@ -c $<

main:	$(OBJS)
	$(CMPI) -o $@ $^ $(LDLIBS)

#parallell version
mainp: $(OBJS)
	$(CMPI) -o $@ $^ $(LDLIBS)

clean:
	rm -f *.o *.vtk main mainp


c1:	main
	./main 0 1 1 2 0.2 10 5 250 250 0 0 2 10 2500 1

c2:	main
	./main 0 1 1 2 0.2 10 5 250 250 0 3 1 10 2500 2

c3:	main
	./main 0.25 1 1.3 3 0.2 15 8 250 250 0 0 3 -20 -5000 3


#c4p is user defined parameters {a, h1, h2, L, tp, nelem_x, nelem_y,kx, ky, kxy, T_edge, q_edge, T0, q0, casenum}
c4:	main
	./main 0.25 1 1.3 8 0.2 50 8 250 250 0 0 3 -20 -5000 4

#parallell versions
c1p:	main
	mpiexec -np 2	./main 0 1 1 2 0.2 10 5 250 250 0 0 2 10 2500 1

c2p:	main
	mpiexec -np 2	./main 0 1 1 2 0.2 10 5 250 250 0 3 1 10 2500 2

c3p:	main
	mpiexec -np 2	./main 0.25 1 1.3 3 0.2 15 8 250 250 0 0 3 -20 -5000 3

#c4p is user defined parameters {a, h1, h2, L, tp, nelem_x, nelem_y,kx, ky, kxy, T_edge, q_edge, T0, q0, casenum}
c4p:	main
	mpiexec -np 2	./main 0.25 1 1.3 3 0.2 15 8 250 250 0 0 3 -20 -5000 4
