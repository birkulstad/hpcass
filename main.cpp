#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>
#include <valarray>
#include "cblas.h"

#define F77NAME(x) x##_
extern "C" {
    // LAPACK routine for solving systems of linear equations
    void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                        const int& lda, int * ipiv, double * B,
                        const int& ldb, int& info);
}

using namespace std;

/* OPTIMISE: stuff that could be optimised
 * HELP: stuff I'm stuck on
 * IMPROVE: stuff that could be implemented in a more elegant/clear way
 * FIX: temporary bugs/problems that need fixing to progress
 *
 */



//Declaring function prototypes
void printMatrix(double*, int, int);
void printMatrix(int*, int, int);

double detMatrix2(double*);
void invMatrix2(double*, double*);

int getInd_i(int, int);
int getInd_j(int, int);

// TODO: IMPLEMENT COMMAND LINE INPUT
//  int argc, char* argv[] inside main()
// a = &argv[1]; //starts at 1
// argc checks correct number of inputs

int main() {
    // --------------- Case Constants -----------------
    // Remove this block once the command line input is enabled
    //      case[] = {a, h1, h2, L, tp, nelem_x, nelem_y,kx, ky, kxy, T_edge, q_edge, T0, q0};
    //const double caseval[] = {0, 1, 1, 2, 0.2, 10, 5, 250, 250, 0, 0, 2, 10, 2500}; //case1
    //const double caseval[] = {0, 1, 1, 2, 0.2, 10, 5, 250, 250, 0, 3, 1, 10, 2500}; //case2
    //const double caseval[] = {0.25, 1, 1.3, 3, 0.2, 15, 8, 250, 250, 0, 0, 3, -20, -5000}; //case3



    // ----- Defining Materials -----------------
    const double kx = caseval[7];  // Thermal conductivity [W/mK]
    const double ky = caseval[8];  // Thermal conductivity [W/mK]
    const double kxy = caseval[9];  // Thermal conductivity [W/mK]

    // ----- Defining Section -----------------
    const double th = caseval[4];  // Thickness [m]

    // ----- Integration scheme -----------------
    const int gaussorder = 2;

    // ----- Defining Geometry -----------------
    // h = a * x**2 + b * x + h1:    Beam height as a function of x
    const double a = caseval[0]; // constant in the polynomial describing the beam height
    const double h1 = caseval[1];  // [m] Height at beam left edge
    const double h2 = h1 * caseval[2]; // [m] Height at beam right edge
    const double L = h1 * caseval[3];  // [m] Beam length

    // calculate b constant in the polynomial describing the beam height
    const double b = -a * L + (h2 - h1) / L;

    // ----- Meshing Geometry -----------------
    const int nelem_x = caseval[5];  // Number of elements in the x-direction
    const int nelem_y = caseval[6];  // Number of elements in y-direction
    const int nnode_x = nelem_x + 1; // Number of nodes in the x-direction
    const int nnode_y = nelem_y + 1; // Number of nodes in the y-direction

    const int nnode_elem = 4;  // number of nodes in each element

    // IMPROVE: Replace nelem_y + 1 with nnode_y and same with nelem_x

    int nNodeDof[nnode_elem] = {1, 1, 1, 1};
    int neDof = 0;
    for (int i : nNodeDof) {
        neDof += i;
    }

    // Basic definition calculations
    const int nelem = nelem_x * nelem_y;  // Total number of elements
    const int nnode = (nelem_x + 1) * (nelem_y + 1);  // Total number of nodes

    // ----- Calculation of Nodal coordinate matrix -> Coord ---------
    // h(x) = a * x**2 + b * x + h1  // Plate height as a function of x

    // OPTIMISE: Create linspace-function and implement where needed
    double x[nelem_x + 1];

    // Distribute length of L over number of nodes (linspace)
    for (int i = 0; i < nelem_x + 1; i++) {
        x[i] = (L - 0) / (nelem_x) * i;
    }
    //printArray(x, nelem_x + 1);


    // Creating x^2 array
    double xsq[nelem_x + 1];
    for (int i = 0; i < nelem_x + 1; i++) {
        *(xsq + i) = x[i] * x[i];
    }

    // Calculating h = a * x^2 + b * x + h1
    double h[nelem_x + 1];

    for (double &i : h) {
        i = h1;
    }
    //printArray(h,nelem_x + 1);

    cblas_dscal(nelem_x + 1, a, xsq, 1);    //multiply x^2 by a -> x2
    //printArray(xsq,nelem_x + 1);

    cblas_daxpy(nelem_x + 1, b, x, 1, xsq, 1); //multiply x by b and add x^2 + a
    //printArray(xsq,nelem_x + 1);

    cblas_daxpy(nelem_x + 1, 1, xsq, 1, h, 1);
    //printArray(h, nelem_x + 1);

    //OPTIMISE: If h is constant (i.e. a and b are zero), Y can be stored as y

    // Defining matrix Y of Y-coordinates for each node
    // Defining matrix X of X-coordinates for each node
    // Defining matrix Coord of [x,y] pairs for each node
    // Calculation of topology matrix NodeTopo
    double X[(nelem_x + 1) * (nelem_y + 1)];
    double Y[(nelem_x + 1) * (nelem_y + 1)];
    double Coord[(nelem_x + 1) * (nelem_y + 1) * 2];
    int NodeTopo[(nelem_x + 1) * (nelem_y + 1)];

    //IMPROVE: Add function/class? that outputs [x,y] coordinates of a [nx,ny] node
    for (int i = 0; i < nelem_y + 1; i++) {
        for (int j = 0; j < nelem_x + 1; j++) {
            X[(nelem_x + 1) * i + j] = x[j];
            Y[(nelem_x + 1) * i + j] = h[j] / (nelem_y) * i - h[j] / 2;
            NodeTopo[(nelem_x + 1) * i + j] = (nelem_y + 1) * j + i;
        }
    }
    //printMatrix(Coord,(nelem_x + 1) * (nelem_y + 1) , 2);
    //printMatrix(NodeTopo,nelem_y + 1,nelem_x + 1);
    //printMatrix(Y,nelem_y + 1,nelem_x + 1);


    //----- Calculation of topology matrix ElemNode -----------
    int ElemNode[(nelem_x * nelem_y) * (nnode_elem + 1)];
    int elemnr = 0;
    for (int i = 0; i < nelem_x + 1; i++) {
        for (int j = 0; j < nelem_y + 1; j++) {
            if ((i < nelem_x) && (j < nelem_y)) {
                ElemNode[5 * (nelem_y * i + j) + 0] = elemnr;    //Element number ID
                ElemNode[5 * (nelem_y * i + j) + 1] = NodeTopo[((nelem_x + 1) * j + i)];            //Upper left node ID
                ElemNode[5 * (nelem_y * i + j) + 2] = NodeTopo[((nelem_x + 1) * j + i + 1)];        //Upper right node ID
                ElemNode[5 * (nelem_y * i + j) + 3] = NodeTopo[((nelem_x + 1) * (j + 1) + i + 1)];  //Lower right node ID
                ElemNode[5 * (nelem_y * i + j) + 4] = NodeTopo[((nelem_x + 1) * (j + 1) + i)];      //Lower left node ID
                elemnr += 1;
            }
        }
    }

    //printMatrix(ElemNode,nelem_x * nelem_y,nnode_elem + 1);
    //printMatrix(NodeTopo,nelem_y + 1,nelem_x + 1);

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
    int nDof = 0; // Number of dof for a given node
    int eDof = 0; // Number of
    for (int i = 0; i < nnode; i++) {
        eDof = globDof[2 * i];
        for (int j = 0; j < eDof; j++) {
            globDof[2 * i + j + 1] = nDof;
            nDof += 1;
        }
    }
    //printMatrix(globDof,nnode,2);

    // Assembly of global stiffness matrix K -------------
    // ----- Gauss-points and weights ----------------------------
    // Gauss order is variable gaussorder
    // Points
    double GP[gaussorder] = {-1 / sqrt(3), 1 / sqrt(3)};
    // Weights
    double W[gaussorder] = {1, 1};

    // ----- Conductivity matrix D  -----------
    double D[4] = {kx, kxy, kxy, ky};



    // ------ Defining variables used in K-assembly --------------
    double K[nnode * nnode]={0};  // Initiation of global stiffness matrix K
    double* Ke = new double[nnode_elem * nnode_elem];   // Initialising element stiffness matrix Ke on heap to allow
                                                        // memset for every new element
    double eCoord[nnode_elem * 2]; // For storing the node coordinates
    double eX[nnode_elem]; // For storing the node x-coordinates
    double eY[nnode_elem]; // For storing the node y-coordinates
    double eta, xi; // Natural coordinates
    double N[nnode_elem]; // Shape function matrix
    double GN[2 * nnode_elem]; // Shape function derivative matrix
    valarray <double> vN(nnode_elem); //valarray temporary class to allow direct assignement of shape function matrix for each iteration
    valarray <double> vGN(2 * nnode_elem); //valarray temporary class to allow direct assignement of derivative shape function matrix for each iteration
    double J[gaussorder * gaussorder]; // Jacobian Matrix
    double DetJ[gaussorder * gaussorder]; // For storing the determinants of J
    double detJ; //Determinant of J
    double invJ[gaussorder * gaussorder]; // Inverse Jacobian Matrix
    double B[2 * nnode_elem]; // Some dot product IMPROVE
    double C[nnode_elem * 2]; // Some dot product IMPROVE
    double Ke_a;

    int p, q;

    for (int i = 0; i < nelem; i++) {
        // - data for element i
        for (int k = 0; k < nnode_elem; k++){
            eNodes[k] = ElemNode[(nnode_elem + 1) * i + k + 1]; // Associated element nodes
            // TODO: Replace p and q with getInd_i and getInd_j
            p = eNodes[k] % (nelem_y + 1);                                  // Row number
            q = (eNodes[k] - (eNodes[k] % (nelem_y + 1))) / (nelem_y + 1);  // Column number
            eCoord[2 * k + 0] = X[(nelem_x + 1) * p + q];  // node x-coordinates
            eCoord[2 * k + 1] = Y[(nelem_x + 1) * p + q];  // node y-coordinates
        }

        //printMatrix(eCoord,nnode_elem,2);

        int gDof[nnode_elem] = {0}; // used to construct scatter matrix

        for (int j = 0; j < nnode_elem; j++) {
            // global dof for node j gathered in element gDof OPTIMISE: gDof == eNodes
            // potential simplification as per above, could have unknown effects in different cases
            gDof[j] = globDof[2 * eNodes[j] + 1];

            // IMPROVE: Are these variables used for anything?
            eX[j] = eCoord[j + 0]; // Node x-coordinates
            eY[j] = eCoord[j + 1]; // Node y-coordinates
        }
        //printMatrix(gDof, 1, nnode_elem);

        memset(Ke, 0.0, sizeof(Ke) * nnode_elem * nnode_elem); // Reset Ke to 0 for each new member
        // ----- Element stiffness matrix, Ke, found by Gauss integration -----------

        for (int j = 0; j < gaussorder; j++){
            for (int k = 0; k < gaussorder; k++){
                eta = GP[j];
                xi = GP[k];
                //cout << "GP = ";
                //printArray(GP,2);
                //cout << "eta, xi = " << eta << "," << xi << endl;
                // shape functions matrix
                vN = {(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                      (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)};
                vN *= 0.25;
                for (int l = 0; l < nnode_elem; l++ ){N[l] = vN[l];};

                // derivative (gradient) of the shape functions
                vGN = {-(1 - eta), 1 - eta, 1 + eta, -(1 + eta),
                -(1 - xi), -(1 + xi), 1 + xi, 1 - xi};
                vGN *= 0.25;
                for (int l = 0; l < 2 * nnode_elem; l++ ){GN[l] = vGN[l];};
                //printMatrix(GN,2, nnode_elem);

                // Jacobian matrix
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 2, 2, 4, 1, GN, 4, eCoord, 2, 0, J, 2);


                // Jacobian determinant
                // TODO: IMPROVE BY IMPLEMENTING LAPACK DETERMINANT
                detJ = detMatrix2(J);

                // TODO: IMPROVE BY IMPLEMENTING LAPACK INVERSE
                // Inverse Jacobian matrix
                invMatrix2(J, invJ);
                //printMatrix(invJ,2,2);

                // Calculating B = invJ * GN
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 2, 4, 2, 1, invJ, 2, GN, 4, 0, B, 4);

                //B = np.dot(invJ, GN)
                //printMatrix(B,2,4);

                // Ke = Ke + B'*D*B*th*DetJ*Wi*Wj
                //Ke = Ke + np.dot(np.dot(B.T, D), B) * th * detJ * W[j] * W[k];

                //Calculating C = B^T * D
                cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans, 4, 2, 2, 1, B, 4, D, 2, 0, C, 2);
                //printMatrix(C,4,2);

                Ke_a = th * detJ * W[j] * W[k]; // BLAS alpha scalar to multiply with the matrix product

                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 4, 4, 2, Ke_a, C, 2, B, 4, 1, Ke, 4);
            }
        }
        //printMatrix(Ke,4,4);
        //printMatrix(J,2,2);
        //cout << "-----------------------" << endl;
        //printMatrix(gDof,1,4);

        // Adding element stiffness matrices Ke to the global stiffness matrix K
        for (int j = 0; j < nnode_elem; j++){
            for (int k = 0; k < nnode_elem; k++) {
                K[nnode * gDof[j] + gDof[k]] += Ke[nnode_elem * j + k];
                //cout << Ke[nnode_elem * j + k] << endl;
            }
        }



    }
    // Release memory from heap
    delete[] Ke;

    //printMatrix(K,nnode, nnode);

    // Case 1
    // Compute nodal boundary flux vector --- natural B.C
    //  Defined on edges
    int q_edge; // Which edge the heat flux is applied to 0 = left, 1 = top, 2 = right, 3 = bottom
    int mult_i, mult_j;
    int nFluxNodes = 0; // Number of nodes on the edge of the beam

    //IMPROVE: Add command line input check and error message for default case

    q_edge = caseval[11];
    switch (q_edge){
        // Left edge
        case 0: nFluxNodes = nelem_y + 1;
                mult_i = (nelem_x + 1);
                mult_j = 0;
                break;

        // Top edge
        case 1: nFluxNodes = nelem_x + 1;
                mult_i = 1;
                mult_j = (nelem_x + 1) * nelem_y;
                break;

        // Right edge
        case 2: nFluxNodes = nelem_y + 1;
                mult_i = (nelem_x + 1);
                mult_j = nelem_x;
                break;

        // Bottom edge
        case 3: nFluxNodes = nelem_x + 1;
                mult_i = 1;
                mult_j = 0;
                break;

        default: cout << "Error! Choose a value between 0 and 3 corresponding to the edge with heat flux!" << endl;
                 cout << "0 = left, 1 = top, 2 = right, 3 = bottom" << endl;
    }

    // Initialise fluxNodes with size appropriate to the edge selected
    int* fluxNodes = new int[nFluxNodes];

    //printMatrix(NodeTopo,nelem_y + 1,nelem_x + 1);
    //cout << nFluxNodes << endl;

    for (int i = 0; i < nFluxNodes; i++){
        fluxNodes[i] = NodeTopo[mult_i * i + mult_j];
    }

    printMatrix(fluxNodes,1,nFluxNodes);



    // ----- Defining load ----------------------------
    double q_flux = caseval[13];  // Constant flux at right edge of the beam
    int nbe = nFluxNodes - 1;  // Number of elements with flux load

    // Element boundary condition
    double* n_bc = new double[4 * nbe];
    for (int i = 0; i < nbe; i++) {
        n_bc[nbe * 0 + i] = fluxNodes[i];  // node 1
        n_bc[nbe * 1 + i] = fluxNodes[i + 1];  // node 2
        n_bc[nbe * 2 + i] = q_flux; // flux value at node 1
        n_bc[nbe * 3 + i] = q_flux;  //flux value at node 2
    }
    //printMatrix(n_bc, 4, 5);

    //printArray(x,nelem_x+1);
    //printMatrix(X,nelem_y+1,nelem_x+1);
    //printMatrix(Y,nelem_y+1,nelem_x+1);

    // --------------- Calculate Nodal Flux Vector f -------------------
    // FIX: Python wants f to be initialised as nDof, but apparently nDof changes?
    double f[nnode] = {0};  // initialize nodal flux vector
    int node1; // first node
    int node2; // second node
    double x1, y1, x2, y2;
    double flux; // WHAT DOES THIS ACTUALLY DO?
    double leng; // length of edge
    double n_bce[2]; // initialising flux values at the two nodes connected to an edge IMPROVE: check that 2 is constant

    //printMatrix(Coord,(nelem_x + 1) * (nelem_y + 1),2);
    double* fq = new double[2]; // initialize the nodal source vector IMPROVE: check that 2 is constant

    for (int i = 0; i < nbe; i++) {

        memset(fq, 0.0, sizeof(fq) * 2); // Reset fq to 0 for each new element

        node1 = fluxNodes[i]; //int(n_bc[nbe * 0 + i]);
        node2 = fluxNodes[i + 1]; //int(n_bc[nbe * 0 + i]);
        n_bce[0] = n_bc[nbe * 2 + i];  // flux value at an edge node 1
        n_bce[1] = n_bc[nbe * 3 + i];  // flux value at an edge node 2
        //cout << "node1 = " << node1 << ", node2 = " << node2 << endl;

        x1 = X[(nelem_x + 1) * getInd_i(node1,nelem_y) + getInd_j(node1,nelem_y)];  // x coord of the first node
        y1 = Y[(nelem_x + 1) * getInd_i(node1,nelem_y) + getInd_j(node1,nelem_y)];  // y coord of the first node
        x2 = X[(nelem_x + 1) * getInd_i(node2,nelem_y) + getInd_j(node2,nelem_y)];  // x coord of the first node
        y2 = Y[(nelem_x + 1) * getInd_i(node2,nelem_y) + getInd_j(node2,nelem_y)];  // y coord of the second node


        //cout << "x1,y1 = " << x1 << "," << y1 << ", x2,y2 = " << x2 << "," << y2 << endl;

        leng = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));  // edge length
        //cout << leng << endl;
        detJ = leng / 2;  // 1D Jacobian

        // IMPROVE: Check that sqrt(nnode_elem) is a good way of getting to 2
        for (int j = 0; j < gaussorder; j++){  // integrate in xi direction (1D integration)

            xi = GP[j];

            // 1D  shape functions in parent domain
            vN = {1 - xi, 1 + xi};
            vN *= 0.5;
            for (int l = 0; l < sqrt(nnode_elem); l++ ){N[l] = vN[l];};
            //printMatrix(N,1,2);

            flux = cblas_ddot(sqrt(nnode_elem), N, 1, n_bce, 1);
            //cout << flux << endl;
            // CBLAS to perform fq += fq + N.T * flux * detJ * th * W[j];  // nodal flux
            cblas_daxpy(sqrt(nnode_elem), flux * detJ * th * W[j], N, 1, fq, 1);

            //printMatrix(fq, 1, 2);
        }
        //printMatrix(fq,1,2);

        //CBLAS to perform fq = -fq;  // define flux as negative integrals
        cblas_dscal(sqrt(nnode_elem),-1,fq,1);
        //printMatrix(fq,1,2);

        f[node1] += fq[0];
        f[node2] += fq[1];

    }
    // Release memory from the heap
    delete[] fluxNodes;
    delete[] n_bc;
    delete[] fq;

    //printMatrix(f, 1, nnode);



    // ----- Apply boundary conditions ----------------- Essential B.C.
    int T_edge; // Which edge the constant Temp is applied to 0 = left, 1 = top, 2 = right, 3 = bottom
    int nTempNodes = 0; // Number of nodes on the edge of the beam
    double T0 = caseval[12];  // Temperature at boundary

    //IMPROVE: Add check to ensure user input does not select same edge for q and T

    T_edge = caseval[10];
    switch (T_edge){
        // Left edge
        case 0: nTempNodes = nelem_y + 1;
            mult_i = (nelem_x + 1);
            mult_j = 0;
            break;

            // Top edge
        case 1: nTempNodes = nelem_x + 1;
            mult_i = 1;
            mult_j = (nelem_x + 1) * nelem_y;
            break;

            // Right edge
        case 2: nTempNodes = nelem_y + 1;
            mult_i = (nelem_x + 1);
            mult_j = nelem_x;
            break;

            // Bottom edge
        case 3: nTempNodes = nelem_x + 1;
            mult_i = 1;
            mult_j = 0;
            break;

        default: cout << "Error! Choose a value between 0 and 3 corresponding to the edge with heat flux!" << endl;
            cout << "0 = left, 1 = top, 2 = right, 3 = bottom" << endl;
    }

    // Initialise fluxNodes with size appropriate to the edge selected
    int* tempNodes = new int[nTempNodes];


    //printMatrix(NodeTopo,nelem_y + 1,nelem_x + 1);
    //cout << nFluxNodes << endl;

    for (int i = 0; i < nTempNodes; i++){
        tempNodes[i] = NodeTopo[mult_i * i + mult_j];
    }

    //FIX: search for new and ensure all variables initialised on the heap are deleted when appropriate

    double* BC = new double[nTempNodes * 2];

    for (int i = 0; i < nTempNodes; i++) {
        BC[2 * i + 0] = tempNodes[i];
        BC[2 * i + 1] = T0;
    }

    //printMatrix(BC,nTempNodes,2);

    // ----- Assembling global "Force" vector ------------
    // FIX: Python wants OrgDof and T to be initialised as nDof, but apparently nDof changes?
    int OrgDof[nnode] = {0};
    double T[nnode] = {0}; // initialize nodal temperature vector
    int rDof = nnode;  // Reduced number of DOF
    //print(rDof)

    //int ind = BC[:, 0].astype(int) //same as tempNodes array
    for (int i = 0; i < nTempNodes; i++){
        OrgDof[tempNodes[i]] = -1;
        T[tempNodes[i]] = BC[2 * i + 1];
        //cout << T[i] << endl;
    }
    rDof = rDof - nTempNodes;

    //printMatrix(OrgDof,1,nnode);

    //FIX: Is any of this block necessary?
    //OrgDof[ind] = -1
    //T[ind] = BC[:, 1]
    //rDof = rDof - nTempNodes;
    //print(rDof)

    //cout << rDof << endl;

    //FIX: Is any of this block necessary?
    /*
    int RedDof[rDof];
    printMatrix(RedDof,1,rDof);
    int counter1 = 0;
    for (int i = 0; i < nnode; i++) {
        RedDof[i] = 0;
        if (OrgDof[i] == 0) {
            OrgDof[j] = counter1;
            RedDof[counter1] = j;
            counter1 += 1;
        }
    }
    */

    //printMatrix(T,1,nnode);

    // IMPROVE: Could elementwise operations using valarrays be used to do masking more efficiently than loops?
    /*
    valarray <double> asdf = {2, 4, 6, 8};
    valarray <double> resulttt;
    resulttt = asdf * asdf;
    printValarray(resulttt,4);
    */

    // --------------------- Partition matrices ----------------------------

    //mask_E = np.array([(i in TempNodes) for i in range(len(T))])  // known temperature Dof
    // IMPROVE: Replace all nnode - nTempNodes with a single variable

    int mask_E[nnode] = {0}; // known temperature Dof
    double T_E[nTempNodes]; // known values of T
    double T_F[nnode-nTempNodes]; // unknown values of T
    double f_E[nTempNodes]; // unknown values of f
    double f_F[nnode - nTempNodes]; // known f-values
    double rhs[nnode - nTempNodes]; // RHS of linear solver


    // Creating mask array for known values. Known -> entry = 1
    for (int i = 0; i < nTempNodes; i++){
        mask_E[tempNodes[i]] = 1;
    }

    //printMatrix(mask_E,1,66);

    // IMPROVE: Perform this in a different manner?
    int counter1 = 0;
    int counter2 = 0;
    int counter3 = 0;
    for (int i = 0; i < nnode; i++){
        if (T[i] != 0){
            T_E[counter1] = T[i];
            counter1++;
        }
        else {
            f_F[counter2] = f[i];
            rhs[counter2] = f[i]; // Pre-setting rhs to f_F, because rhs = f_F + ...cblas routine
            counter2++;
        }
    }
    //printMatrix(rhs, 1, nnode - nTempNodes);
    //printMatrix(T_E,1,nTempNodes);
    //printMatrix(f_F,1,nnode-nTempNodes);

    // Initialising partition matrices
    double K_EE[nTempNodes * nTempNodes];
    double K_EF[nTempNodes * (nnode - nTempNodes)];
    double K_FF[(nnode - nTempNodes) * (nnode - nTempNodes)];

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
    //printMatrix(K_EF,nTempNodes,(nnode - nTempNodes));
    //printMatrix(K_FF,(nnode - nTempNodes),(nnode - nTempNodes));

    // TODO: Implement LAPACK solver for this block

    // -------------------- Solve for T_F -------------------------
    // Calculate RHS
    //rhs = f_F - np.dot(K_EF.T, T_E)
    //printMatrix(rhs, 1, nnode - nTempNodes);
    cblas_dgemv(CblasRowMajor, CblasTrans, nTempNodes, nnode - nTempNodes, -1, K_EF, nnode - nTempNodes, T_E, 1, 1, rhs, 1);
    //printMatrix(rhs, 1, nnode - nTempNodes);

    //T_F = np.linalg.solve(K_FF, rhs)
    int ipiv[nnode - nTempNodes];
    int info;
    F77NAME(dgesv)(nnode - nTempNodes, 1, K_FF, nnode - nTempNodes, ipiv, rhs, nnode - nTempNodes, info);

    cblas_dcopy(nnode - nTempNodes, rhs, 1, T_F, 1);
    //printMatrix(T_F, 1, nnode - nTempNodes);

    // compute the reaction f_E
    double f_E1[nTempNodes];


    // Calculating f_E = np.dot(K_EE, T_E) + np.dot(K_EF, T_F)

    // Calculating f_E = K_EE * T_E
    cblas_dgemv(CblasRowMajor, CblasNoTrans, nTempNodes, nTempNodes, 1, K_EE, nTempNodes, T_E, 1, 0, f_E, 1);
    //printMatrix(f_E,1,nTempNodes);

    // Calculating f_E1 = K_EF * T_F
    cblas_dgemv(CblasRowMajor, CblasNoTrans, nTempNodes, nnode - nTempNodes, 1, K_EF, nnode - nTempNodes, T_F, 1, 0, f_E1, 1);
    //printMatrix(T_F, 1, nnode - nTempNodes);

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
    printMatrix(T,nelem_x + 1, nelem_y + 1);
    //printMatrix(f,nelem_x + 1, nelem_y + 1);

    // Cleaning remaining dynamic arrays from heap

    delete[] tempNodes;
    delete[] BC;

    return 0;
}

int getInd_i(int node, int nrow){
    int i = node % (nrow + 1); // Row number
    return i;
}

int getInd_j(int node, int nrow){
    int j = (node - (node % (nrow + 1))) / (nrow + 1);  // Column number
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
