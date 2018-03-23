CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3
OBJS = main.o
LDLIBS = -llapack -lblas

default:main
all:	main clean

%.o:	%.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

main:	$(OBJS)
	$(CXX) -o $@ $^ $(LDLIBS)



clean:
	rm -f *.o *.vtk main

# {a, h1, h2, L, tp, nelem_x, nelem_y,kx, ky, kxy, T_edge, q_edge, T0, q0}
c1:	main
	./main 0 1 1 2 0.2 10 5 250 250 0 0 2 10 2500

c2:	main
	./main 0 1 1 2 0.2 10 5 250 250 0 3 1 10 2500

c3:	main
	./main 0.25 1 1.3 3 0.2 15 8 250 250 0 0 3 -20 -5000
