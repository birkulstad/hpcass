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

/* OPTIMISE: stuff that could be optimised
 * HELP: stuff I'm stuck on
 * IMPROVE: stuff that could be implemented in a more elegant/clear way
 * FIX: temporary bugs/problems that need fixing to progress
 *
 */



//Declaring function prototypes
void printVector(vector<double>);
void printVector(vector<int>);

void printArray(double*, int);
void printArray(int*, int);
void printMatrix(double*, int, int);
void printMatrix(int*, int, int);

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
    // nNodeDof[nnode_elem] = {1, 1, 1, 1};
    //double neDof = cblas_dasum(nnode_elem,nNodeDof,1);

    int nNodeDof[nnode_elem] = {1, 1, 1, 1};
    int neDof = 0;
    for (int i : nNodeDof) {
        neDof += i;
    }


    //const vector<int> nNodeDof(nnode_elem, 1);
    //const int neDof = accumulate(nNodeDof.begin(), nNodeDof.end(), 0);


    //printVector(nNodeDof);
    //cout << neDof << endl;

// Calculatations
    const int nelem = nelem_x * nelem_y;  // Total number of elements
    const int nnode = (nelem_x + 1) * (nelem_y + 1);  // Total number of nodes

// ----- Calculation of Nodal coordinate matrix -> Coord ---------
    // h(x) = a * x**2 + b * x + h1  # Beam height as a function of x
    //OPTIMISE: Create linspace-function and implement where needed
    double x[nelem_x + 1];
    //Distribute length of L over number of nodes (linspace)
    for (int i = 0; i < nelem_x + 1; i++) {
        x[i] = (L - 0) / (nelem_x) * i;
    }

    //cout << x << endl;
    //printArray(x, nelem_x + 1);

    //double x2[nelem_x];
    //creating x^2 array
    double x2[nelem_x + 1];
    for (int i = 0; i < nelem_x; i++) {
        *(x2 + i) = x[i] * x[i];
    }

    //calculating h = a * x^2 + b * x + h1
    double h[nelem_x + 1];

    for (double &i : h) {
        i = h1;
    }
    //printArray(*h,nelem_x);
    cblas_dscal(nelem_x + 1, a, x2, 1);    //multiply x^2 by a -> x2
    //printArray(*x2,nelem_x);
    cblas_daxpy(nelem_x + 1, b, x, 1, x2, 1); //multiply x by b and add x^2 + a
    //printArray(*x2,nelem_x);
    cblas_daxpy(nelem_x + 1, 1, x2, 1, h, 1);
    //printArray(h, nelem_x + 1);

    //double h[nelem_x];
    //double h = a * x2 + b * x + h1

    //OPTIMISE: If h is constant (i.e. a and b are zero), Y can be stored as y
    double y[nelem_y + 1];
    //Distribute length of h(x) over number of nodes in y (linspace)
    for (int j = 0; j < nelem_y + 1; j++) {
        y[j] = h[j] * j / (nelem_y) - h[j] / 2;
    }
    //printArray(y, nelem_y + 1);


    //Defining matrix Y of Y-coordinates for each node
    //Defining matrix Coord of [x,y] pairs for each node
    //Calculation of topology matrix NodeTopo
    double X[(nelem_x + 1) * (nelem_y + 1)];
    double Y[(nelem_x + 1) * (nelem_y + 1)];
    double Coord[(nelem_x + 1) * (nelem_y + 1) * 2];
    int NodeTopo[(nelem_x + 1) * (nelem_y + 1)];

    //IMPROVE: Add function/class? that outputs [x,y] coordinates of a [nx,ny] node
    for (int i = 0; i < nelem_y + 1; i++) {
        for (int j = 0; j < nelem_x + 1; j++) {
            X[(nelem_x + 1) * i + j] = x[j];
            Y[(nelem_x + 1) * i + j] = y[i] * h[j];
            NodeTopo[(nelem_x + 1) * i + j] = (nelem_y + 1) * j + i;
            //cout << "i,j = " << i << "," << j << "y,h=" << y[i] << "," << h[j] << "Y=" << Y[nelem_y * i + j] << endl;
            Coord[2 * ((nelem_x + 1) * i + j) + 0] = x[j];
            Coord[2 * ((nelem_x + 1) * i + j) + 1] = y[i] * h[j];
            //cout << "i,j = " << i << "," << j << "    x,y = " << Coord[2 * (nelem_y * i + j) + 0] << "," << Coord[2 * (nelem_y * i + j) + 1] << "    Y = " << Y[nelem_y * i + j] << endl;
            //cout << (nelem_y + 1) * j + i << endl;
        }
    }
    //printMatrix(X,6,11);


    //----- Calculation of topology matrix ElemNode -----------
    int ElemNode[(nelem_x * nelem_y) * (nnode_elem + 1)];
    int elemnr = 0;
    for (int i = 0; i < nelem_x + 1; i++) {
        for (int j = 0; j < nelem_y + 1; j++) {
            if ((i < nelem_x) && (j < nelem_y)) {
                ElemNode[5 * (nelem_y * i + j) + 0] = elemnr;    //Element number ID
                ElemNode[5 * (nelem_y * i + j) + 1] = NodeTopo[((nelem_x + 1) * j + i)]; //Upper left node ID
                ElemNode[5 * (nelem_y * i + j) + 2] = NodeTopo[((nelem_x + 1) * j + i + 1)]; //Upper right node ID
                ElemNode[5 * (nelem_y * i + j) + 3] = NodeTopo[((nelem_x + 1) * (j + 1) + i + 1)]; //Lower right node ID
                ElemNode[5 * (nelem_y * i + j) + 4] = NodeTopo[((nelem_x + 1) * (j + 1) + i)]; //Lower left node ID
                elemnr += 1;
            }
        }
    }

    //printMatrix(ElemNode,50,5);
    //printMatrix(NodeTopo,6,11);
    //printArray(NodeTopo, 6*11);

    double ElemX[(nelem_x * nelem_y) * (nnode_elem)]; // Element nodal x-coordinates
    double ElemY[(nelem_x * nelem_y) * (nnode_elem)]; // Element nodal y-coordinates
    int eNodes[nnode_elem];
    int p, q;   // Nodal indices (for a row-major eg. X or Y representation)
                // calculated for each node associated with an element

    //int s[] = {nelem_y+1,1}; //stride for X stored in column-major format
    //X(i,j) = *X+i*s[0]+j*s[1]

    for (int i = 0; i < nelem; i++){
        for (int k = 0; k < nnode_elem; k++) {
            eNodes[k] = ElemNode[(nnode_elem + 1) * i + k + 1]; // Associated element nodes
            // Calculating indices (i,j)=(p,q) to obtain (x,y)-coordinates from X and Y
            // Note that p and q are calculated bearing in mind NodeTopo numbering of nodes
            p = eNodes[k] % (nelem_y + 1);  // Row number
            q = (eNodes[k] - (eNodes[k] % (nelem_y + 1))) / (nelem_y + 1); // Column number
            //cout << "p,q = " << p << "," << q << "    np+q = " << (nelem_x + 1) * p + q << "    X[] = " << X[(nelem_x + 1) * p + q] << endl;
            ElemX[(nnode_elem) * i + k] = X[(nelem_x + 1) * p + q];
            ElemY[(nnode_elem) * i + k] = Y[(nelem_x + 1) * p + q];
        }
        //printArray(eNodes,4);
    }
    //printMatrix(ElemX,50,4);

    // ----- Generate global dof numbers -----------------

    int globDof[nnode * 2]={0};
    int nNode;

    for (int i = 0; i < nelem; i++){
        for (int k = 0; k < nnode_elem; k++) {
            nNode = ElemNode[(nnode_elem + 1) * i + k + 1];
            // if the already existing ndof of the present node is less than
            // the present elements ndof then replace the ndof for that node]

            if (globDof[2 * nNode] < nNodeDof[k]) {
                globDof[2 * nNode] = nNodeDof[k];
            }
        }
    }
    //printMatrix(globDof,50,2);
    //printArray(nNodeDof,nnode_elem);

    // counting the global dofs and inserting in globDof
    int nDof = 0;
    int eDof = 0;
    for (int i = 0; i < nnode; i++) {
        eDof = globDof[2 * i];
        for (int j = 0; j < eDof; j++) {
            globDof[2 * i + j + 1] = nDof;
            nDof += 1;
        }
    }
    printMatrix(globDof,nnode,2);

    // Assembly of global stiffness matrix K -------------
    // ----- Gauss-points and weights ----------------------------
    // Gauss order is variable gaussorder
    // Points
    double GP[gaussorder] = {-1 / sqrt(3), 1 / sqrt(3)};
    // Weights
    double W[gaussorder] = {1, 1};

    // ----- Conductivity matrix D  -----------
    double D[4] = {kx, kxy, kxy, ky};
    // ----------------------------------------
    double Kp[nDof * nDof];  // Initiation of global stiffness matrix K
    double Ke[nnode_elem * nnode_elem]; // Initialising element stiffness matrix Ke
    double eCoord[nnode_elem * 2]; // For storing the node coordinates
    double eX[nnode_elem]; // For storing the node x-coordinates
    double eY[nnode_elem]; // For storing the node y-coordinates
    double eta, xi; // Natural coordinates
    double N[nnode_elem]; // Shape function matrix
    double GN[2 * nnode_elem]; // Shape function derivative matrix
    double J[gaussorder * gaussorder]; // Jacobian Matrix
    double DetJ[gaussorder * gaussorder]; // For storing the determinants of J
    double detJ; //Determinant of J
    double invJ[gaussorder * gaussorder]; // Inverse Jacobian Matrix
    double B[2 * nnode_elem]; // Some dot product IMPROVE





    for (int i = 0; i < nelem; i++) {
        // - data for element i
        for (int k = 0; k < nnode_elem; k++){
            eNodes[k] = ElemNode[(nnode_elem + 1) * i + k + 1]; // Associated element nodes
            p = eNodes[k] % (nelem_y + 1);                                  // Row number
            q = (eNodes[k] - (eNodes[k] % (nelem_y + 1))) / (nelem_y + 1);  // Column number
            eCoord[2 * k + 0] = X[(nelem_x + 1) * p + q];  // node x-coordinates
            eCoord[2 * k + 1] = Y[(nelem_x + 1) * p + q];  // node y-coordinates
        }


        double gDof[nnode_elem] = {0}; // used to construct scatter matrix

        for (int j = 0; j < nnode_elem; j++) {
            // global dof for node j gathered in element gDof OPTIMISE: gDof == eNodes
            // potential simplification as per above, could have unknown effects in different cases
            gDof[j] = globDof[2 * eNodes[j] + 1];

            eX[j] = eCoord[j + 0]; // Node x-coordinates
            eY[j] = eCoord[j + 1]; // Node y-coordinates
        }


        // ----- Element stiffness matrix, Ke, found by Gauss integration -----------

        for (int j = 0; j < gaussorder; j++){
            for (int k = 0; k < gaussorder; k++){
                eta = GP[i];
                xi = GP[j];

                // FIX elegant way to assign N, GN etc without allocating new memory every time
                // maybe define N1, N2 etc and assign to N
                // shape functions matrix
                N = {0.25 * (1 - xi) * (1 - eta), 0.25 * (1 + xi) * (1 - eta),
                        0.25 * (1 + xi) * (1 + eta), 0.25 * (1 - xi) * (1 + eta)};

                // derivative (gradient) of the shape functions
                GN = {0.25 * -(1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), 0.25 * -(1 + eta),
                0.25 * -(1 - xi), 0.25 * -(1 + xi), 0.25 * (1 + xi), 0.25 * (1 - xi)};

                // Jacobian matrix
                J = np.dot(GN, eCoord)

                // Jacobian determinant
                DetJ = np.linalg.det(J)

                // Inverse Jacobian matrix
                invJ = np.linalg.inv(J)

                B = np.dot(invJ, GN)

                // Ke = Ke + B'*D*B*th*DetJ*Wi*Wj
                Ke = Ke + np.dot(np.dot(B.T, D), B) * th * DetJ * W[i] * W[j]
            }
        }



    }


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
        cout << setw(5) << *(a+i);
        if (i < length - 1){
            cout << ", ";
        };
    }
    cout << "]" << endl;
}

void printArray(int* a, int length){
    cout << "[";
    for (int i = 0; i < length; i++){
        cout << setw(5) << *(a+i);
        if (i < length - 1){
            cout << ", ";
        };
    }
    cout << "]" << endl;
}

void printMatrix(double* a, int M, int N){
    cout << "[";
    for (int i = 0; i < M; ++i){
        for (int j = 0; j < N; ++j){
            if (i + j >= 1){
                cout << ", ";
            };
            if ((j == 0) && (i > 0)){
                cout << endl;
            }
            cout << setw(5) << *(a + i * N + j);

        }

    }
    cout << "]" << endl;
}

void printMatrix(int* a, int M, int N){
    cout << "[";
    for (int i = 0; i < M; ++i){
        for (int j = 0; j < N; ++j){
            if (i + j >= 1){
                cout << ", ";
            };
            if ((j == 0) && (i > 0)){
                cout << endl;
            }
            cout << setw(5) << *(a + i * N + j);

        }

    }
    cout << "]" << endl;
}
