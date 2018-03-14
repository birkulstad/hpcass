#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>

#include "cblas.h"

using namespace std;

//Declaring function prototypes
void printVector(vector<double>);
void printVector(vector<int>);

void printArray(double*, int);






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
    const double a = 0; // constant in the polynomial describing the beam height
    const double h1 = 1;  // [m] Height at beam left edge
    const double h2 = h1 * 1; // [m] Height at beam right edge
    const double L = 2 * h1;  // [m] Beam length

// calculate b constant in the polynomial describing the beam height
    const double b = -a * L + (h2 - h1) / L;

// ----- Meshing Geometry -----------------
    const int nelem_y = 5;  // Number of elements in y-direction
    const int nelem_x = 10;  // Number of elements in the x-direction

    const int nnode_elem = 4;  // number of nodes in each element

    //OPTIMISE:
    //const double nNodeDof[] = {1.0,1.0,1.0,1.0};
    //const double neDof = cblas_dasum(nnode_elem,nNodeDof,1);

    const vector<double> nNodeDof(nnode_elem, 1);
    const int neDof = accumulate(nNodeDof.begin(),nNodeDof.end(),0);


    //printVector(nNodeDof);
    //cout << neDof << endl;

// Calculatations
    const int nelem = nelem_x * nelem_y;  // Total number of elements
    const int nnode = (nelem_x + 1) * (nelem_y + 1);  // Total number of nodes

// ----- Calculation of Nodal coordinate matrix -> Coord ---------
    // h(x) = a * x**2 + b * x + h1  # Beam height as a function of x
    //OPTIMISE: Create linspace-function and implement where needed
    double x[nelem_x];

    for (double &i : x) {
        i = (L - 0) / (nelem_x + 1);
    }

    //cout << x << endl;
    //printArray(*x,nelem_x);

    //double x2[nelem_x];
    double x2[nelem_x];
    for (int i = 0; i < nelem_x; i++){
        *(x2 + i) = x[i] * x[i];
    }

    //calculating h = a * x^2 + b * x + h1
    double h[nelem_x];

    for (double &i : h) {
        i = h1;
    }
    //printArray(*h,nelem_x);
    cblas_dscal(nelem_x,a,x2,1);    //multiply x^2 by a -> x2
    //printArray(*x2,nelem_x);
    cblas_daxpy(nelem_x, b, x, 1, x2, 1); //multiply x by b and add x^2 + a
    //printArray(*x2,nelem_x);
    cblas_daxpy(nelem_x,1,x2,1,h,1);
    //printArray(*h,nelem_x);

    //double h[nelem_x];
    //double h = a * x2 + b * x + h1

// generating zero array for each element in the mesh
    /*
    double Y[(nelem_x + 1) * (nelem_x + 1)] = {0};
    double* YP = &Y[0];
    Y[2] = 5;
    cout << Y[2] << endl;
     */
    double Y[(nelem_x + 1) * (nelem_y + 1)];
    double* YP = &Y[0];
    //std::fill(Y,(nelem_x + 1) * (nelem_x + 1), 0 );
    printArray(YP,(nelem_x + 1) * (nelem_x + 1));
    cout << YP << endl;
    cout << Y[0] << endl;
    //out << YP << endl;


    return 0;
}


void printVector(vector<double> Vector) {
    const int length = Vector.size();
    cout << "[";
    for (int i = 0; i < length; i++) {
        cout << Vector[i];
        if (i < length - 1){
            cout << ", ";
        }
    }
    cout << "]" << endl;
}

void printVector(vector<int> Vector) {
    const int length = Vector.size();
    cout << "[";
    for (int i = 0; i < length; i++) {
        cout << Vector[i] << ", ";
    }
    cout << "]" << endl;
}

void printArray(double* a, int length){
    cout << "[";
    for (int i = 0; i < length; i++){
        cout << *(a + i);
        if (i < length - 1){
            cout << ", ";
        };
    }
    cout << "]" << endl;
}

