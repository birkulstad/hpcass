CXX = g++
CXXFLAGS = -Wall
LDLIBS -llapack - lblas -I.

default all:

FE_HT: FE_HT.cpp
	$(CXX) $(CXXFLAGS) -o FE_HT.o FE_HT.cpp $(LSLIBS)

compile: FE_HT.o
	$(CXX) -o all FE_HT.o $(LDLIBS)

all:
	$(CXX) $(CXXFLAGS) -o all FE_HT.o $(LDLIBS)

clean:
	rm -f *.o FE_HT

# {a, h1, h2, L, tp, nelem_x, nelem_y,kx, ky, kxy, T_edge, q_edge, T0, q0}
c1: all
	./all 0 1 1 2 0.2 10 5 250 250 0 0 2 10 2500

c2: all
	./all 0 1 1 2 0.2 10 5 250 250 0 3 1 10 2500

c2: all
	./all 0.25 1 1.3 3 0.2 15 8 250 250 0 0 3 -20 -5000
