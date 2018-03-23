#include <cstring>
#include <iostream>
#include <cmath>
#include <valarray>
#include <mpi.h>

#include "cblas.h"
#include "io.h"
#include "init.h"
#include "oper.h"

using namespace std;

// Written by Birk Andreas Ulstad, Imperial College London, Department of Aeronautics March 2018
// NOTE: All arrays are stored in RowMajor format, dynamic arrays are used sparsely and quickly delete from the heap.


int main(int argc, char* argv[]) {

    // Initialising MPI-variables
    int rank, size, retval_rank, retval_size;
    MPI_Init(&argc, &argv);
    retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    retval_size = MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM){
        cout << "Invalid Communicator!" << endl;
    }
    cout << "I am process " << rank + 1 << " of " << size << endl;

    // ----- Integration scheme -----------------
    const int gaussorder = 2;

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

    // ----- Defining BC Edges -------------------
    const int T_edge = atoi(argv[11]); // Which edge the constant Temp is applied to 0 = left, 1 = top, 2 = right, 3 = bottom
    const int q_edge = atoi(argv[12]); // Which edge the heat flux is applied to 0 = left, 1 = top, 2 = right, 3 = bottom

    // ----- Defining BC Parameters -------------------
    double T0 = atof(argv[13]);  // Constant temperature at edge
    double q0 = atof(argv[14]);  // Constant flux at right edge of the plate

    int casenum = atoi(argv[15]); // Case number for determining .vtk filename

    // ----------------------- Checking inputs, throw errorMessage if invalid ---------------------
    if (h1 <= 0 || h2 <= 0 || L <= 0 || th <= 0){errorMessage(1); return -1;}
    if (nelem_x <= 1 || nelem_y <= 1){errorMessage(2); return -1;}
    if (kx * ky <= 0){errorMessage(3); return -1;}
    if ((T_edge < 0 || T_edge > 3) || (q_edge < 0 || q_edge > 3)){errorMessage(4); return -1;}
    if (T_edge == q_edge){errorMessage(5); return -1;}
    // if (a > 2 * (h1 + h2)/(L * (2 - L))){errorMessage(6); return -1;} // Non-continuous surface condition?

    // Calculate b constant in the polynomial describing the plate height
    const double b = -a * L + (h2 - h1) / L;

    const int nnode_x = nelem_x + 1; // Number of nodes in the x-direction
    const int nnode_y = nelem_y + 1; // Number of nodes in the y-direction

    const int nnode_elem = 4;  // Number of nodes in each element (Quadrilateral)
    int nNodeDof[nnode_elem] = {1, 1, 1, 1}; // Degrees of freedom per node

    // Basic definition calculations
    const int nelem = nelem_x * nelem_y;  // Total number of elements
    const int nnode = nnode_x * nnode_y;  // Total number of nodes

    // Discretising x
    double x[nnode_x];
    // Distribute length of L over number of nodes
    linspace(0, L, nnode_x, x);

    // ----- Calculation of plate height h(x) ---------
    // h = a * x^2 + b * x + h1

    // Creating x^2 array
    double xsq[nnode_x];
    for (int i = 0; i < nnode_x; i++) {
        *(xsq + i) = x[i] * x[i];
    }

    // Calculating h = a * x^2 + b * x + h1
    double h[nnode_x];

    // Adding h1 to all entries before blas operations
    for (int i = 0; i < nnode_x; i++) {
        h[i] = h1;
    }

    // Performing calculations of h using CBLAS Level 1
    cblas_dscal(nnode_x, a, xsq, 1);        // multiply x^2 by a -> x^22
    cblas_daxpy(nnode_x, b, x, 1, xsq, 1);  // multiply x by b and add x^2 + a
    cblas_daxpy(nnode_x, 1, xsq, 1, h, 1);  // add ax^2 + bx to h1 in h

    // Coordinate matrices X and Y are mapped in the same way as Node Topology Matrix (NodeTopo)
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

    // Calculate list of node [x,y]-coordinates for .vtk-output
    double Coord[nnode * 2];
    for (int i = 0; i < nnode; i++){
        Coord[2 * i + 0] = getCoord(i, X, nnode_x, nnode_y);  // node x-coordinate
        Coord[2 * i + 1] = getCoord(i, Y, nnode_x, nnode_y);  // node y-coordinate
    }

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

    // ----- Generate global dof numbers -----------------
    int eNodes[nnode_elem];         // Nodes associated with a given element
    int globDof[nnode * 2];         // List of [number of Dof for node, global Dof]
    populate(globDof, nnode * 2, 0);

    for (int i = 0; i < nelem; i++){
        for (int k = 0; k < nnode_elem; k++) {
            eNodes[k] = ElemNode[(nnode_elem + 1) * i + k + 1];
            if (globDof[2 * eNodes[k]] < nNodeDof[k]) {
                globDof[2 * eNodes[k]] = nNodeDof[k];
            }
        }
    }

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

    // --------------------------------- Assembly of global stiffness matrix K ----------------------------------
    // ----- Gauss-points and weights ----------------------------
    const double GP[gaussorder] = {-1 / sqrt(3), 1 / sqrt(3)};      // Gaussian Points
    const double W[gaussorder] = {1, 1};                            // Gaussian Weights

    // Conductivity matrix D
    const double D[2 * 2] = {kx, kxy, kxy, ky};

    // ----------- Defining variables used in K-assembly --------------

    // Variables used for storing matrices/related to matrix operations
    auto* Ke = new double[nnode_elem * nnode_elem];
    double K[nnode * nnode];  // Initiation of global stiffness matrix K
    double B[2 * nnode_elem]; // Temporary matrix for cblas operations
    double C[nnode_elem * 2]; // Temporary matrix for cblas operations
    double Ke_a;    // alpha-multiplier for cblas calculation of Ke
    populate(K, nnode * nnode, 0);

    // Variables associated with the the natural coordinate system and the Jacobian
    double eta, xi; // Natural coordinates
    double N[nnode_elem]; // Shape function matrix
    double GN[2 * nnode_elem]; // Shape function derivative matrix
    double J[gaussorder * gaussorder]; // Jacobian Matrix
    double invJ[gaussorder * gaussorder]; // Inverse Jacobian Matrix
    double detJ; //Determinant of J

    // Coordinates of nodes associated with the element
    double eCoord[nnode_elem * 2]; // For storing the node coordinates of an element


    // Using valarray class to initialise N and GN for every element iteration
    // Allows direct assignment of shape function matrix
    valarray <double> vN(nnode_elem);
    valarray <double> vGN(2 * nnode_elem);

    for (int i = 0; i < nelem; i++) {
        for (int k = 0; k < nnode_elem; k++){
            eNodes[k] = ElemNode[(nnode_elem + 1) * i + k + 1]; // Associated element nodes
            eCoord[2 * k + 0] = getCoord(eNodes[k], X, nnode_x, nnode_y);  // x-coordinate of associated element nodes
            eCoord[2 * k + 1] = getCoord(eNodes[k], Y, nnode_x, nnode_y);  // y-coordinate of associated element nodes
        }

        memset(Ke, 0, sizeof(Ke) * nnode_elem * nnode_elem); // Reset Ke to 0 for each new member
        // ----- Element stiffness matrix, Ke, found by Gauss integration -----------

        for (int j = 0; j < gaussorder; j++){
            for (int k = 0; k < gaussorder; k++){
                eta = GP[j];
                xi = GP[k];

                // Shape functions matrix
                vN = {(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                      (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)};  // One-operation assignment using valarray
                vN *= 0.25;
                for (int l = 0; l < nnode_elem; l++ ){N[l] = vN[l];};   // Copying values to normal array

                // Derivative (gradient) matrix of the shape functions
                vGN = {-(1 - eta), 1 - eta, 1 + eta, -(1 + eta),
                -(1 - xi), -(1 + xi), 1 + xi, 1 - xi};              // One-operation assignment using valarray
                vGN *= 0.25;
                for (int l = 0; l < 2 * nnode_elem; l++ ){GN[l] = vGN[l];}; // Copying values to normal array

                // Jacobian matrix GN * eCoord
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 2, 2, 4, 1, GN, 4, eCoord, 2, 0, J, 2);
                detJ = detMatrix2(J);   // Jacobian determinant (2x2 determinant not worth doing using LAPACK P-L-U)
                invMatrix2(J, invJ);    // Inverse Jacobian matrix (2x2 inverse not worth doing using LAPACK P-L-U)

                // Calculating B = invJ * GN
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 2, 4, 2, 1, invJ, 2, GN, 4, 0, B, 4);

                // Ke = Ke + [B' * D * B] * th * detJ * W[j] * W[k]
                Ke_a = th * detJ * W[j] * W[k]; // BLAS alpha scalar to multiply with the matrix product

                //Calculating C = B' * D
                cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans, 4, 2, 2, 1, B, 4, D, 2, 0, C, 2);

                // Calculating Ke += Ke_a * [C * B]
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 4, 4, 2, Ke_a, C, 2, B, 4, 1, Ke, 4);
            }
        }
        // -------------------------------- Assembling K-matrix --------------------------------------
        // Adding element stiffness matrices Ke to the global stiffness matrix K
        for (int j = 0; j < nnode_elem; j++){
            for (int k = 0; k < nnode_elem; k++) {
                K[nnode * eNodes[j] + eNodes[k]] += Ke[nnode_elem * j + k];
            }
        }
    }
    // Release dynamic memory from heap
    delete[] Ke;

    // ------------------ Apply boundary conditions -----------------
    // Compute nodal boundary flux vector f --- natural B.C
    auto* fluxNodes = new int[max(nnode_x, nnode_y)];   // Always >= than the largest dimension, to avoid seg fault
    int nFluxNodes = 0;                                 // Number of nodes on the chosen edge of the plate

    // Obtaining node numbers along edge based on edge selected
    getBCNodes(q_edge, nnode_x, nnode_y, nFluxNodes, fluxNodes, NodeTopo);

    // ----- Defining flux load -----------------------
    int nbe = nFluxNodes - 1;  // Number of elements with flux load

    // Element boundary condition
    auto* n_bc = new double[4 * nbe];
    for (int i = 0; i < nbe; i++) {
        n_bc[nbe * 0 + i] = fluxNodes[i];       // node 1
        n_bc[nbe * 1 + i] = fluxNodes[i + 1];   // node 2
        n_bc[nbe * 2 + i] = q0;                 // flux value at node 1
        n_bc[nbe * 3 + i] = q0;                 //flux value at node 2
    }

    // --------------- Calculate Nodal Flux Vector f -------------------
    double f[nnode];            // initialize nodal flux vector
    int node1;                  // first node along an element edge
    int node2;                  // second node along an element edge
    double x1, y1, x2, y2;      // Coordinates of the two nodes connected to an element edge
    double flux;                // Flux at an element edge
    double leng;                // Length of element edge
    double n_bce[2];            // initialising flux values at the two nodes connected to an element edge
    populate(f, nnode, 0);

    auto* fq = new double[2]; // initialize the nodal source vector

    for (int i = 0; i < nbe; i++) {

        memset(fq, 0, sizeof(fq) * 2); // Reset fq to 0 for each new element

        node1 = fluxNodes[i];           // Defining node 1 for the edge element
        node2 = fluxNodes[i + 1];       // Defining node 2 for the edge element
        n_bce[0] = n_bc[nbe * 2 + i];   // Flux value at node 1 of the edge element
        n_bce[1] = n_bc[nbe * 3 + i];   // Flux value at node 2 of the edge element

        // Obtaining coordinates for the two nodes of the edge element using indices
        x1 = getCoord(node1,X,nnode_x, nnode_y);  // x-coord of the first node
        y1 = getCoord(node1,Y,nnode_x, nnode_y);  // y-coord of the first node
        x2 = getCoord(node2,X,nnode_x, nnode_y);  // x-coord of the second node
        y2 = getCoord(node2,Y,nnode_x, nnode_y);  // y-coord of the second node

        // Calculating length of element edge
        leng = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));  // Element edge length
        detJ = leng / 2;  // 1D Jacobian

        for (int j = 0; j < gaussorder; j++){  // Integrate in xi direction (1D integration)

            xi = GP[j];

            // 1D  shape functions in parent domain
            vN = {1 - xi, 1 + xi};
            vN *= 0.5;
            for (int l = 0; l < sqrt(nnode_elem); l++ ){N[l] = vN[l];};

            // Calculating flux = [N * n_bce]
            flux = cblas_ddot(2, N, 1, n_bce, 1);

            // CBLAS to perform fq += fq + N.T * flux * detJ * th * W[j];  // nodal flux
            cblas_daxpy(2, flux * detJ * th * W[j], N, 1, fq, 1);
        }

        //CBLAS to perform fq = -fq (Defining flux as negative integrals)
        cblas_dscal(2,-1,fq,1);

        // Adding element flux to flux vector f
        f[node1] += fq[0];
        f[node2] += fq[1];

    }
    // Release dynamic memory from the heap
    delete[] fluxNodes;
    delete[] n_bc;
    delete[] fq;

    // ------------------ Apply boundary conditions ----------------- Essential B.C.
    // Initialise fluxNodes with size appropriate to the edge selected
    auto* tempNodes = new int[max(nnode_x, nnode_y)];   // Always >= than the largest dimension, to avoid seg fault
    int nTempNodes = 0;                                 // Number of nodes on the edge of the plate

    // Obtaining node numbers along edge based on edge selected
    getBCNodes(T_edge, nnode_x, nnode_y, nTempNodes, tempNodes, NodeTopo);

    // Initialising matrix containing nodes along chosen edge and their constant temps
    auto* BC = new double[nTempNodes * 2];

    for (int i = 0; i < nTempNodes; i++) {
        BC[2 * i + 0] = tempNodes[i];
        BC[2 * i + 1] = T0;
    }

    // ----------- Assembling global temperature vector ------------
    double T[nnode];  // Initialize nodal temperature vector
    populate(T, nnode, 0);

    // Setting temperatures of nodes along chosen edge to T0
    for (int i = 0; i < nTempNodes; i++){
        T[tempNodes[i]] = BC[2 * i + 1];
    }

    // --------------------- Partition matrices ----------------------------
    // IMPROVE: Increase efficiency by masking using element-wise operations rather than loops

    const int rDof = nnode - nTempNodes;// Number of nodes not affected by temp BC
    int mask_E[nnode];                  // Known temperature Dof
    double T_E[nTempNodes];             // Known values of T
    double T_F[nnode-nTempNodes];       // Unknown values of T
    double f_E[nTempNodes];             // Unknown values of f
    double f_F[rDof];                   // Known values of f
    double rhs[rDof];                   // RHS of linear solver
    populate(mask_E, nnode, 0);


    // Creating mask array for known values. If the entry of the array is known -> entry = 1
    for (int i = 0; i < nTempNodes; i++){
        mask_E[tempNodes[i]] = 1;
    }

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

    // Pre-setting rhs to f_F, because rhs = f_F + ...cblas routine
    cblas_dcopy(rDof, f_F, 1, rhs, 1);

    // Calculating RHS = f_F - [K_EF^T * T_E]
    cblas_dgemv(CblasRowMajor, CblasTrans, nTempNodes, rDof, -1, K_EF, rDof, T_E, 1, 1, rhs, 1);

    // Initialising variables required by LAPACK to solve linear system
    int ipiv[rDof];
    int info;

    // Solving the linear system [K_FF * T_F] = RHS  to obtain T_F
    F77NAME(dgesv)(rDof, 1, K_FF, rDof, ipiv, rhs, rDof, info); // Result is saved to rhs

    // Copying output from LAPACK to the correct variable T_F
    cblas_dcopy(rDof, rhs, 1, T_F, 1);

    // ---------------------- Calculating the reaction f_E --------------------
    double f_E1[nTempNodes]; // Temporary array to store values for BLAS operations

    // Calculating f_E = [K_EE * T_E] + [K_EF * T_F]

    // Calculating f_E = K_EE * T_E
    cblas_dgemv(CblasRowMajor, CblasNoTrans, nTempNodes, nTempNodes, 1, K_EE, nTempNodes, T_E, 1, 0, f_E, 1);

    // Calculating f_E1 = K_EF * T_F
    cblas_dgemv(CblasRowMajor, CblasNoTrans, nTempNodes, rDof, 1, K_EF, rDof, T_F, 1, 0, f_E1, 1);

    // Calculating f_E = f_E + f_E1
    cblas_daxpy(nTempNodes,1,f_E1,1,f_E,1);

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

    // Output result to .vtk-file to plot using ParaView
    WriteVtkFile(nnode_elem, nnode, nelem, Coord, ElemNode, T, casenum);

    // Cleaning remaining dynamic arrays from heap
    delete[] tempNodes;
    delete[] BC;

    MPI_Finalize();

    return 0;
}

