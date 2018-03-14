#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iomanip>

#include "cblas.h"

using namespace std;







int main() {
// ----- Defining Materials -----------------
    const double kx = 250.0;  // Thermal conductivity [W/mK]
    const double ky = 250.0;  // Thermal conductivity [W/mK]
    const double kxy = 0.0;  // Thermal conductivity [W/mK]

// ----- Defining Section -----------------
    const double th = .2;  // Thickness [m]

// ----- Integration scheme -----------------
    const int gaussorder = 2;

// ----- Defining Geometry -----------------
    // h = a * x**2 + b * x + h1:    Beam height as a function of x
    const double a = 0.25; // constant in the polynomial describing the beam height
    const double h1 = 1;  // [m] Height at beam left edge
    const double h2 = h1 * 1.3; // [m] Height at beam right edge
    const double L = 3 * h1;  // [m] Beam length

// calculate b constant in the polynomial describing the beam height
    const double b = -a * L + (h2 - h1) / L;

// ----- Meshing Geometry -----------------
    const int nelem_y = 8;  // Number of elements in y-direction
    const int nelem_x = 15;  // Number of elements in the x-direction

    const int nnode_elem = 4;  // number of nodes in each element

    int nNodeDof[4] = {};
    int neDof = 0;

    for (int i; i < nnode_elem; i++){   //OPTIMISE: Change to vector and use cblas to sum
        nNodeDof[i] = 1;    //number of DoF per node (1 = Temperature only)
        neDof += nNodeDof[i];// total number of DoF per element
    }


// Calculatations
    const int nelem = nelem_x * nelem_y;  // Total number of elements
    const int nnode = (nelem_x + 1) * (nelem_y + 1);  // Total number of nodes


    return 0;
}