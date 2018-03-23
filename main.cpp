#include <iostream>
#include <cmath>
#include <iomanip>
#include <valarray>
#include <fstream>
#include "cblas.h"
#include <limits>
#include <cstring> // Need this for memset on Linux Remote FIX

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>
#include <ctime>

#define F77NAME(x) x##_
extern "C" {
    // LAPACK routine for solving systems of linear equations
    void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                        const int& lda, int * ipiv, double * B,
                        const int& ldb, int& info);
}

using namespace std;

/*
 * OPTIMISE: stuff that could be optimised
 * HELP: stuff I'm stuck on
 * IMPROVE: stuff that could be implemented in a more elegant/clear way
 * FIX: temporary bugs/problems that need fixing to progress
 */

//Declaring function prototypes
void printMatrix(double*, int, int);
void printMatrix(int*, int, int);

double detMatrix2(double*);
void invMatrix2(double*, double*);

int getInd_i(int, int);
int getInd_j(int, int);

void WriteVtkFile(int, int, int, double* , int* , double* );

// TODO: IMPLEMENT COMMAND LINE INPUT
//  int argc, char* argv[] inside main()
// a = &argv[1]; //starts at 1
// argc checks correct number of inputs

int main(int argc, char* argv[]) {
    // IMPROVE: Diagonalise K-matrix
    // FIX: Setting matrices to 0 at initialisation, so 0 entries are in fact zero if possible


    // --------------- Case Constants -----------------
    // Remove this block once the command line input is enabled

    //const double caseval[] = {a, h1, h2, L, tp, nelem_x, nelem_y,kx, ky, kxy, T_edge, q_edge, T0, q0};
     const double caseval[] = {0, 1, 1, 2, 0.2, 10, 5, 250, 250, 0, 0, 2, 10, 2500}; //case1
    //const double caseval[] = {0, 1, 1, 2, 0.2, 10, 5, 250, 250, 0, 3, 1, 10, 2500}; //case2
    //const double caseval[] = {0.25, 1, 1.3, 3, 0.2, 15, 8, 250, 250, 0, 0, 3, -20, -5000}; //case3

    // ----- Integration scheme -----------------
    const int gaussorder = 2;

    // ----- Defining Geometry -----------------
    const double a = caseval[0]; // constant in the polynomial describing the plate height
    const double h1 = caseval[1];  // [m] Height at plate left edge
    const double h2 = h1 * caseval[2]; // [m] Height at plate right edge
    const double L = h1 * caseval[3];  // [m] Plate length

    // ----- Defining Section -----------------
    const double th = caseval[4];  // Thickness [m]

    // ----- Meshing Geometry -----------------
    const auto nelem_x = int(caseval[5]);  // Number of elements in the x-direction
    const auto nelem_y = int(caseval[6]);  // Number of elements in y-direction

    // ----- Defining Materials -----------------
    const double kx = caseval[7];   // Thermal conductivity [W/mK]
    const double ky = caseval[8];   // Thermal conductivity [W/mK]
    const double kxy = caseval[9];  // Thermal conductivity [W/mK]

    // ----- Boundary Conditions -----------------
    const auto T_edge = int(caseval[10]); // Which edge the constant Temp is applied to 0 = left, 1 = top, 2 = right, 3 = bottom
    const auto q_edge = int(caseval[11]); // Which edge the heat flux is applied to 0 = left, 1 = top, 2 = right, 3 = bottom

    // ----- Boundary Conditions Parameters -----------------
    double T0 = caseval[12];  // Constant temperature at the chosen edge of the plate
    double q0 = caseval[13];  // Constant flux at the chosen edge of the plate

    // ##########################################################################################################################
    /*
    // ----- Defining Geometry -----------------
    const double a = atof(argv[1]); // constant in the polynomial describing the plate height
    const double h1 = atof(argv[2]);  // [m] Height at plate left edge
    const double h2 = h1 * atof(argv[3]); // [m] Height at plate right edge
    const double L = h1 * atof(argv[4]);  // [m] Plate length

    // ----- Defining Section -----------------
    const double th = atof(argv[5]);  // Thickness [m]

    // ----- Meshing Geometry -----------------
    const int nelem_x = atoi(argv[6]);  // Number of elements in the x-direction
    const int nelem_y = atoi(argv[7]);  // Number of elements in y-direction

    // ----- Defining Materials -----------------
    const double kx = atof(argv[8]);   // Thermal conductivity [W/mK]
    const double ky = atof(argv[9]);   // Thermal conductivity [W/mK]
    const double kxy = atof(argv[10]);  // Thermal conductivity [W/mK]

    const int T_edge = atoi(argv[11]); // Which edge the constant Temp is applied to 0 = left, 1 = top, 2 = right, 3 = bottom
    const int q_edge = atoi(argv[12]); // Which edge the heat flux is applied to 0 = left, 1 = top, 2 = right, 3 = bottom

    double T0 = atof(argv[13]);  // Constant temperature at edge
    double q0 = atof(argv[14]);  // Constant flux at right edge of the plate
    */



    // calculate b constant in the polynomial describing the plate height
    const double b = -a * L + (h2 - h1) / L;

    const int nnode_x = nelem_x + 1; // Number of nodes in the x-direction
    const int nnode_y = nelem_y + 1; // Number of nodes in the y-direction

    const int nnode_elem = 4;  // Number of nodes in each element (Quadrilateral)
    int nNodeDof[nnode_elem] = {1, 1, 1, 1}; // Degrees of freedom per node

    // Basic definition calculations
    const int nelem = nelem_x * nelem_y;  // Total number of elements
    const int nnode = nnode_x * nnode_y;  // Total number of nodes

    // ----- Calculation of plate height h(x) ---------
    // h = a * x^2 + b * x + h1
    // OPTIMISE: Create linspace-function and implement where needed
    double x[nnode_x];

    // Distribute length of L over number of nodes (linspace)
    for (int i = 0; i < nnode_x; i++) {
        x[i] = (L - 0) / (nelem_x) * i;
    }
    //printArray(x, nnode_x);

    // Creating x^2 array
    double xsq[nnode_x];
    for (int i = 0; i < nnode_x; i++) {
        *(xsq + i) = x[i] * x[i];
    }

    // Calculating h = a * x^2 + b * x + h1
    double h[nnode_x];

    // Adding h1 to all entries before blas operations
    for (double &i : h) {
        i = h1;
    }
    //printArray(h,nnode_x);

    cblas_dscal(nnode_x, a, xsq, 1);    //multiply x^2 by a -> x2
    //printArray(xsq,nnode_x);

    cblas_daxpy(nnode_x, b, x, 1, xsq, 1); //multiply x by b and add x^2 + a
    //printArray(xsq,nnode_x);

    cblas_daxpy(nnode_x, 1, xsq, 1, h, 1);
    //printArray(h, nnode_x);

    // OPTIMISE: If h is constant (i.e. a and b are zero), Y can be stored as y
    // Coordinate matrices X and Y are mapped in the same way as Node Topology Matrix
    double X[nnode_x * nnode_y];     // Defining matrix X of X-coordinates for each node
    double Y[nnode_x * nnode_y];     // Defining matrix Y of Y-coordinates for each node
    int NodeTopo[nnode_x * nnode_y]; // Calculation of topology matrix NodeTopo

    // Calculating the X, Y and NodeTopo
    for (int i = 0; i < nelem_y + 1; i++) {
        for (int j = 0; j < nnode_x; j++) {
            X[nnode_x * i + j] = x[j];
            Y[nnode_x * i + j] = h[j] / (nelem_y) * i - h[j] / 2;
            NodeTopo[nnode_x * i + j] = nnode_y * j + i;
        }
    }
    //printMatrix(NodeTopo,nnode_y,nnode_x);
    //printMatrix(Y,nnode_y,nnode_x);

    // Calculate list of node [x,y]-coordinates for .vtk-output
    double Coord[nnode * 2];
    for (int i = 0; i < nnode; i++){
        Coord[2 * i + 0] = X[nnode_x * getInd_i(i,nnode_y) + getInd_j(i,nnode_y)];  // node x-coordinate
        Coord[2 * i + 1] = Y[nnode_x * getInd_i(i,nnode_y) + getInd_j(i,nnode_y)];  // node y-coordinate
    }
    //printMatrix(Coord,nnode,2);

    //----- Calculation of topology matrix ElemNode -----------
    int ElemNode[nelem * (nnode_elem + 1)];
    int elemnr = 0;
    for (int i = 0; i < nnode_x; i++) {
        for (int j = 0; j < nnode_y; j++) {
            if ((i < nelem_x) && (j < nelem_y)) {
                ElemNode[5 * (nelem_y * i + j) + 0] = elemnr;    //Element number ID
                ElemNode[5 * (nelem_y * i + j) + 1] = NodeTopo[(nnode_x * j + i)];            //Upper left node ID
                ElemNode[5 * (nelem_y * i + j) + 2] = NodeTopo[(nnode_x * j + i + 1)];        //Upper right node ID
                ElemNode[5 * (nelem_y * i + j) + 3] = NodeTopo[(nnode_x * (j + 1) + i + 1)];  //Lower right node ID
                ElemNode[5 * (nelem_y * i + j) + 4] = NodeTopo[(nnode_x * (j + 1) + i)];      //Lower left node ID
                elemnr += 1;
            }
        }
    }

    //printMatrix(ElemNode,nelem_x * nelem_y,nnode_elem + 1);
    //printMatrix(NodeTopo,nnode_y,nnode_x);

    // ----- Generate global dof numbers -----------------
    int eNodes[nnode_elem];         // Nodes associated with a given element
    int globDof[nnode * 2]={0};     // List of [number of Dof for node, global Dof]

    for (int i = 0; i < nelem; i++){
        for (int k = 0; k < nnode_elem; k++) {
            eNodes[k] = ElemNode[(nnode_elem + 1) * i + k + 1];
            if (globDof[2 * eNodes[k]] < nNodeDof[k]) {
                globDof[2 * eNodes[k]] = nNodeDof[k];
            }
        }
        //printArray(eNodes,4);
    }
    //printMatrix(globDof,nnode,2);

    // Counting the global dofs and inserting in globDof
    int nDof = 0; // Number of global dofs active at a given node
    int eDof = 0; // Number of dofs per node at a given node
    for (int i = 0; i < nnode; i++) {
        eDof = globDof[2 * i];
        for (int j = 0; j < eDof; j++) {
            globDof[2 * i + j + 1] = nDof;
            nDof += 1;
        }
    }
    //printMatrix(eDof,nnode,2);


    // Assembly of global stiffness matrix K -------------
    // ----- Gauss-points and weights ----------------------------

    // Gaussian Points
    const double GP[gaussorder] = {-1 / sqrt(3), 1 / sqrt(3)};

    // Gaussian Weights
    const double W[gaussorder] = {1, 1};

    // Conductivity matrix D
    const double D[2 * 2] = {kx, kxy, kxy, ky};



    // ------ Defining variables used in K-assembly --------------

    // Variables used for storing matrices/related to matrix operations
    auto* Ke = new double[nnode_elem * nnode_elem];   // Initialising element stiffness matrix Ke on heap to allow
                                                        // memset for every new element
    double K[nnode * nnode]={0};  // Initiation of global stiffness matrix K
    double B[2 * nnode_elem]; // Temporary matrix for cblas operations
    double C[nnode_elem * 2]; // Temporary matrix for cblas operations
    double Ke_a;    // alpha-multiplier for cblas calculation of Ke

    // Variables associated with the the natural coordinate system and the Jacobian
    double eta, xi; // Natural coordinates
    double N[nnode_elem]; // Shape function matrix
    double GN[2 * nnode_elem]; // Shape function derivative matrix
    double J[gaussorder * gaussorder]; // Jacobian Matrix
    double invJ[gaussorder * gaussorder]; // Inverse Jacobian Matrix
    double detJ; //Determinant of J

    // Coordinates of nodes associated with the element
    int gDof[nnode_elem] = {0}; // Used to construct scatter matrix
    double eCoord[nnode_elem * 2]; // For storing the node coordinates of an element


    // Using valarray class to initialise N and GN for every element iteration
    // Allows direct assignment of shape function matrix
    valarray <double> vN(nnode_elem);
    valarray <double> vGN(2 * nnode_elem);

    for (int i = 0; i < nelem; i++) {
        for (int k = 0; k < nnode_elem; k++){
            eNodes[k] = ElemNode[(nnode_elem + 1) * i + k + 1]; // Associated element nodes

            // eCoord is defined by obtaining the [i,j]-indices a node using getInd_i and getInd_j
            // These [i,j] indices are then used to extract the relevant element coordinates from X and Y
            eCoord[2 * k + 0] = X[nnode_x * getInd_i(eNodes[k],nnode_y) + getInd_j(eNodes[k],nnode_y)];  // node x-coordinates
            eCoord[2 * k + 1] = Y[nnode_x * getInd_i(eNodes[k],nnode_y) + getInd_j(eNodes[k],nnode_y)];  // node y-coordinates
        }
        //printMatrix(eCoord,nnode_elem,2);

        for (int j = 0; j < nnode_elem; j++) {
            // global dof for node j gathered in element gDof OPTIMISE: gDof == eNodes
            // potential simplification as per above, could have unknown effects in different cases
            gDof[j] = globDof[2 * eNodes[j] + 1];
        }
        //printMatrix(gDof, 1, nnode_elem);

        memset(Ke, 0.0, sizeof(Ke) * nnode_elem * nnode_elem); // Reset Ke to 0 for each new member
        // ----- Element stiffness matrix, Ke, found by Gauss integration -----------

        for (int j = 0; j < gaussorder; j++){
            for (int k = 0; k < gaussorder; k++){
                eta = GP[j];
                xi = GP[k];
                //printArray(GP,2);

                // Shape functions matrix
                vN = {(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                      (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)};
                vN *= 0.25;
                for (int l = 0; l < nnode_elem; l++ ){N[l] = vN[l];};   // Copying values to normal array

                // Derivative (gradient) matrix of the shape functions
                vGN = {-(1 - eta), 1 - eta, 1 + eta, -(1 + eta),
                -(1 - xi), -(1 + xi), 1 + xi, 1 - xi};
                vGN *= 0.25;
                for (int l = 0; l < 2 * nnode_elem; l++ ){GN[l] = vGN[l];}; // Copying values to normal array
                //printMatrix(GN,2, nnode_elem);

                // Jacobian matrix GN * eCoord
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 2, 2, 4, 1, GN, 4, eCoord, 2, 0, J, 2);
                //printMatrix(J, gaussorder, gaussorder);

                // Jacobian determinant (2x2 determinant not worth doing using LAPACK P-L-U)
                detJ = detMatrix2(J);

                // Inverse Jacobian matrix (2x2 inverse not worth doing using LAPACK P-L-U)
                invMatrix2(J, invJ);
                //printMatrix(invJ,2,2);

                // Calculating B = invJ * GN
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 2, 4, 2, 1, invJ, 2, GN, 4, 0, B, 4);
                //printMatrix(B,2,4);

                // Ke = Ke + [B' * D * B] * th * DetJ * Wi * Wj
                // Ke = Ke + np.dot(np.dot(B.T, D), B) * th * detJ * W[j] * W[k];

                //Calculating C = B^T * D
                cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans, 4, 2, 2, 1, B, 4, D, 2, 0, C, 2);
                //printMatrix(C,4,2);

                Ke_a = th * detJ * W[j] * W[k]; // BLAS alpha scalar to multiply with the matrix product

                // Calculating Ke += Ke_a * [C * B]
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 4, 4, 2, Ke_a, C, 2, B, 4, 1, Ke, 4);
            }
        }
        //printMatrix(Ke,nnode_elem, nnode_elem);
        //printMatrix(J,2,2);

        // Adding element stiffness matrices Ke to the global stiffness matrix K
        for (int j = 0; j < nnode_elem; j++){
            for (int k = 0; k < nnode_elem; k++) {
                K[nnode * gDof[j] + gDof[k]] += Ke[nnode_elem * j + k];
                //cout << Ke[nnode_elem * j + k] << endl;
            }
        }



    }

    // Release dynamic memory from heap
    delete[] Ke;

    //printMatrix(K,nnode, nnode);


    // ------------------ Apply boundary conditions -----------------
    // Compute nodal boundary flux vector f --- natural B.C
    int mult_i, mult_j; // Index multiplier defined by which edge is loaded, used to record edge nodes
    int nFluxNodes = 0; // Number of nodes on the chosen edge of the Plate

    // IMPROVE: Add command line input check and error message for default case
    // IMPROVE: Add case with no flux boundary?

    // Assigning number of fluxNodes and index multipliers based on edge chosen
    switch (q_edge){
        // Left edge
        case 0: nFluxNodes = nnode_y;
                mult_i = nnode_x;
                mult_j = 0;
                break;

        // Top edge
        case 1: nFluxNodes = nnode_x;
                mult_i = 1;
                mult_j = nnode_x * nelem_y;
                break;

        // Right edge
        case 2: nFluxNodes = nnode_y;
                mult_i = nnode_x;
                mult_j = nelem_x;
                break;

        // Bottom edge
        case 3: nFluxNodes = nnode_x;
                mult_i = 1;
                mult_j = 0;
                break;

        default: cout << "Error! Choose a value between 0 and 3 corresponding to the edge with heat flux!" << endl;
                 cout << "0 = left, 1 = top, 2 = right, 3 = bottom" << endl;
    }

    // Initialise fluxNodes with size appropriate to the edge chosen
    auto* fluxNodes = new int[nFluxNodes];

    // Recording which nodes are located along the edge chosen
    for (int i = 0; i < nFluxNodes; i++){
        fluxNodes[i] = NodeTopo[mult_i * i + mult_j];
    }
    //printMatrix(fluxNodes,1,nFluxNodes);


    // ----- Defining flux load -----------------------
    int nbe = nFluxNodes - 1;  // Number of elements with flux load

    // Element boundary condition
    auto* n_bc = new double[4 * nbe];
    for (int i = 0; i < nbe; i++) {
        n_bc[nbe * 0 + i] = fluxNodes[i];  // node 1
        n_bc[nbe * 1 + i] = fluxNodes[i + 1];  // node 2
        n_bc[nbe * 2 + i] = q0; // flux value at node 1
        n_bc[nbe * 3 + i] = q0;  //flux value at node 2
    }
    //printMatrix(n_bc, 4, 5);


    // --------------- Calculate Nodal Flux Vector f -------------------
    double f[nnode] = {0};  // initialize nodal flux vector
    int node1; // first node along an element edge
    int node2; // second node along an element edge
    double x1, y1, x2, y2; // Coordinates of the two nodes connected to an element edge
    double flux; // Flux at an element edge
    double leng; // Length of element edge
    double n_bce[2]; // initialising flux values at the two nodes connected to an element edge

    auto* fq = new double[2]; // initialize the nodal source vector

    for (int i = 0; i < nbe; i++) {

        memset(fq, 0.0, sizeof(fq) * 2); // Reset fq to 0 for each new element

        node1 = fluxNodes[i]; // Defining node 1 for the edge element
        node2 = fluxNodes[i + 1]; // Defining node 2 for the edge element
        n_bce[0] = n_bc[nbe * 2 + i];  // Flux value at node 1 of the edge element
        n_bce[1] = n_bc[nbe * 3 + i];  // Flux value at node 2 of the edge element
        //cout << "node1 = " << node1 << ", node2 = " << node2 << endl;

        // Obtaining coordinates for the two nodes of the edge element using indices
        x1 = X[nnode_x * getInd_i(node1,nnode_y) + getInd_j(node1,nnode_y)];  // x coord of the first node
        y1 = Y[nnode_x * getInd_i(node1,nnode_y) + getInd_j(node1,nnode_y)];  // y coord of the first node
        x2 = X[nnode_x * getInd_i(node2,nnode_y) + getInd_j(node2,nnode_y)];  // x coord of the first node
        y2 = Y[nnode_x * getInd_i(node2,nnode_y) + getInd_j(node2,nnode_y)];  // y coord of the second node

        // Calculating length of element edge
        leng = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));  // Element edge length
        detJ = leng / 2;  // 1D Jacobian

        // IMPROVE: Check that sqrt(nnode_elem) is a good way of getting to 2
        for (int j = 0; j < gaussorder; j++){  // integrate in xi direction (1D integration)

            xi = GP[j];

            // 1D  shape functions in parent domain
            vN = {1 - xi, 1 + xi};
            vN *= 0.5;
            for (int l = 0; l < sqrt(nnode_elem); l++ ){N[l] = vN[l];};
            //printMatrix(N,1,2);

            flux = cblas_ddot(2, N, 1, n_bce, 1);
            //cout << flux << endl;

            // CBLAS to perform fq += fq + N.T * flux * detJ * th * W[j];  // nodal flux
            cblas_daxpy(2, flux * detJ * th * W[j], N, 1, fq, 1);
            //printMatrix(fq, 1, 2);
        }
        //printMatrix(fq,1,2);

        //CBLAS to perform fq = -fq (Defining flux as negative integrals)
        cblas_dscal(2,-1,fq,1);
        //printMatrix(fq,1,2);

        // Adding element flux to flux vector f
        f[node1] += fq[0];
        f[node2] += fq[1];

    }
    //printMatrix(f, 1, nnode);

    // Release dynamic memory from the heap
    delete[] fluxNodes;
    delete[] n_bc;
    delete[] fq;


    // ------------------ Apply boundary conditions ----------------- Essential B.C.
    int nTempNodes = 0; // Number of nodes on the edge of the plate

    // IMPROVE: Add check to ensure user input does not select same edge for q and T?
    // IMPROVE: Add case with no flux boundary?

    // Assigning number of tempNodes and index multipliers based on edge chosen
    switch (T_edge){
        // Left edge
        case 0: nTempNodes = nnode_y;
            mult_i = nnode_x;
            mult_j = 0;
            break;

            // Top edge
        case 1: nTempNodes = nnode_x;
            mult_i = 1;
            mult_j = nnode_x * nelem_y;
            break;

            // Right edge
        case 2: nTempNodes = nnode_y;
            mult_i = nnode_x;
            mult_j = nelem_x;
            break;

            // Bottom edge
        case 3: nTempNodes = nnode_x;
            mult_i = 1;
            mult_j = 0;
            break;

        default: cout << "Error! Choose a value between 0 and 3 corresponding to the edge with heat flux!" << endl;
            cout << "0 = left, 1 = top, 2 = right, 3 = bottom" << endl;
    }

    // Initialise fluxNodes with size appropriate to the edge selected
    auto* tempNodes = new int[nTempNodes];

    // Recording which nodes are located along the edge chosen
    for (int i = 0; i < nTempNodes; i++){
        tempNodes[i] = NodeTopo[mult_i * i + mult_j];
    }
    //printMatrix(tempNodes,1,nTempNodes);

    // Initialising matrix containing nodes along chosen edge and their constant temps
    auto* BC = new double[nTempNodes * 2];

    for (int i = 0; i < nTempNodes; i++) {
        BC[2 * i + 0] = tempNodes[i];
        BC[2 * i + 1] = T0;
    }
    //printMatrix(BC,nTempNodes,2);

    // ----------- Assembling global temperature vector ------------
    double T[nnode] = {0}; // Initialize nodal temperature vector

    // Setting temperatures of nodes along chosen edge to T0
    for (int i = 0; i < nTempNodes; i++){
        T[tempNodes[i]] = BC[2 * i + 1];
        //cout << T[i] << endl;
    }


    // --------------------- Partition matrices ----------------------------
    // IMPROVE: Increase efficiency by masking using element-wise operations rather than loops

    const int rDof = nnode - nTempNodes;// Number of nodes not affected by temp BC
    int mask_E[nnode] = {0};            // Known temperature Dof
    double T_E[nTempNodes];             // Known values of T
    double T_F[nnode-nTempNodes];       // Unknown values of T
    double f_E[nTempNodes];             // Unknown values of f
    double f_F[rDof];                   // Known values of f
    double rhs[rDof];                   // RHS of linear solver


    // Creating mask array for known values. If the entry of the array is known -> entry = 1
    for (int i = 0; i < nTempNodes; i++){
        mask_E[tempNodes[i]] = 1;
    }
    //printMatrix(mask_E,1,66);

    // Initialising arbitrary counters
    int counter1 = 0;
    int counter2 = 0;
    int counter3 = 0;

    // Recording known values of f and T to partitions T_E and f_F
    for (int i = 0; i < nnode; i++){
        // If T contains a value not equal to zero record to T_E
        if (T[i] != 0){
            T_E[counter1] = T[i];
            counter1++;
        }
        // Otherwise record f to f_F
        else {
            f_F[counter2] = f[i];
            counter2++;
        }
    }
    //printMatrix(rhs, 1, rDof);
    //printMatrix(T_E,1,nTempNodes);
    //printMatrix(f_F,1,rDof);

    // Initialising partition matrices
    double K_EE[nTempNodes * nTempNodes];
    double K_EF[nTempNodes * rDof];
    double K_FF[rDof * rDof];

    counter1 = 0;
    counter2 = 0;
    counter3 = 0;

    for (int j = 0; j < nnode; j++){
        for (int k = 0; k < nnode; k++) {
            if ((mask_E[j] == 1) && (mask_E[k] == 1)){
                K_EE[counter1] = K[nnode * j + k];
                counter1++;
            }
            else if ((mask_E[j] == 1) && (mask_E[k] != 1)){
                K_EF[counter2] = K[nnode * j + k];
                counter2++;
            }
            else if ((mask_E[j] != 1) && (mask_E[k] != 1)){
                K_FF[counter3] = K[nnode * j + k];
                counter3++;

            }
        }
    }
    //printMatrix(K,nnode,nnode);
    //printMatrix(K_EE,nTempNodes,nTempNodes);
    //printMatrix(K_EF,nTempNodes,rDof);
    //printMatrix(K_FF,rDof,rDof);

    // Pre-setting rhs to f_F, because rhs = f_F + ...cblas routine
    cblas_dcopy(rDof, f_F, 1, rhs, 1);

    // Calculating RHS = f_F - [K_EF^T * T_E]
    //printMatrix(rhs, 1, rDof);
    cblas_dgemv(CblasRowMajor, CblasTrans, nTempNodes, rDof, -1, K_EF, rDof, T_E, 1, 1, rhs, 1);
    //printMatrix(rhs, 1, rDof);

    // Initialising variables required by LAPACK to solve linear system
    int ipiv[rDof];
    int info;

    // Solving the linear system [K_FF * T_F] = RHS  to obtain T_F
    F77NAME(dgesv)(rDof, 1, K_FF, rDof, ipiv, rhs, rDof, info); // Result is saved to rhs

    // Copying output from LAPACK to the correct variable T_F
    cblas_dcopy(rDof, rhs, 1, T_F, 1);
    //printMatrix(T_F, 1, rDof);



    // ---------------------- Calculating the reaction f_E --------------------
    double f_E1[nTempNodes]; // Temporary array to store values for BLAS operations


    // Calculating f_E = [K_EE * T_E] + [K_EF * T_F]

    // Calculating f_E = K_EE * T_E
    cblas_dgemv(CblasRowMajor, CblasNoTrans, nTempNodes, nTempNodes, 1, K_EE, nTempNodes, T_E, 1, 0, f_E, 1);
    //printMatrix(f_E,1,nTempNodes);

    // Calculating f_E1 = K_EF * T_F
    cblas_dgemv(CblasRowMajor, CblasNoTrans, nTempNodes, rDof, 1, K_EF, rDof, T_F, 1, 0, f_E1, 1);
    //printMatrix(T_F, 1, rDof);

    // Calculating f_E = f_E + f_E1
    cblas_daxpy(nTempNodes,1,f_E1,1,f_E,1);
    //printMatrix(f_E, 1, nTempNodes);

    counter1 = 0;
    counter2 = 0;

    for (int i = 0; i < nnode; i++) {
        if (mask_E[i] == 1) {
            T[i] = T_E[counter1];
            f[i] = f_E[counter1];
            counter1++;
        }
        else {
            T[i] = T_F[counter2];
            f[i] = f_F[counter2];
            counter2++;
        }
    }
    printMatrix(T,nnode_x, nnode_y);
    //printMatrix(f,nnode_x, nnode_y);


    // Output result to .vtk-file to plot using ParaView
    WriteVtkFile(nnode_elem, nnode, nelem, Coord, ElemNode, T);



    // Cleaning remaining dynamic arrays from heap
    delete[] tempNodes;
    delete[] BC;

    return 0;
}


int getInd_i(int node, int nnode_y){
    int i = node % (nnode_y); // Nodal Row number
    return i;
}

int getInd_j(int node, int nnode_y){
    int j = (node - (node % (nnode_y))) / (nnode_y);  // Nodal Column number
    return j;
}


double detMatrix2(double* A){
    double detJ;
    detJ = A[3] * A[0] - A[2] * A[1];
    return detJ;
}


void invMatrix2(double* J, double* invJ){
    double detJ = detMatrix2(J);
    invJ[0] = J[3] / detJ;
    invJ[1] = -J[1] / detJ;
    invJ[2] = -J[2] / detJ;
    invJ[3] = J[0] / detJ;
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
            cout << setw(8) << *(a + i * N + j);

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


/// @brief Creates a VTK file of the solutions
void WriteVtkFile(int nnode_elem, int nnode, int nelem, double* Coord, int* ElemNode, double* T) {
/** asdfgsdfg */
    ofstream vOut("Output.vtk", ios::out | ios::trunc); // Output file initialisation

    vOut.precision(10);

    // Title and info
    vOut << "# vtk DataFile Version 4.0" << endl; // VTK Version
    vOut << "vtk output" << endl;   // Stating output
    vOut << "ASCII" << endl;    // Format
    vOut << "DATASET UNSTRUCTURED_GRID" << endl;    // Dataset
    vOut << "POINTS " << nnode << " double" << endl; // Data type

    // Print nodal coordinates
    for (int i = 0; i < nnode ; i++) {
        // [x,y,z] (=0), separated by spaces
        vOut << Coord[2 * i] << " " << Coord[2 * i + 1] << " " << "0 ";
    }
    vOut << endl;

    vOut << "Cells " << nelem << " " << (nnode_elem + 1) * nelem << endl;

    // Print nodes of each element
    for (int i = 0; i < nelem; i++) {
        vOut << nnode_elem << " ";
        for (int j = 1; j < nnode_elem + 1; j++) {
            vOut << ElemNode[(nnode_elem + 1) * i + j] << " ";
        }
        vOut << endl;
    }

    // Print cell types as quadrilateral elements
    vOut << "CELL_TYPES " << nelem << endl;
    for (int i = 0; i < nelem; i++) {
        vOut << "9 " << endl;
    }
    // Adding metadata
    vOut << "POINT_DATA " << nnode << endl;
    vOut << "FIELD FieldData 1" << endl;
    vOut << "disp 1 " << nnode << " double" << endl;

    // Outputting temperature vector T
    for (int i = 0; i < nnode; i++) {
        vOut << T[i] << " ";
    }
    vOut.close();

}