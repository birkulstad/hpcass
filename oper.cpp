#include "oper.h"
#include <cmath>

///@brief obtains the i-index of a node, based on the position of the node's number in the [nnode_y x nnode_x] NodeTopo matrix
int getInd_i(int node, int nnode_y){
    /**
     * The i index is based on a row major storage convention in the matrix of which the index will be used to extract from
     * e.g. A[rowlength * i + j]
     * @param node      the number of the node in the NodeTopo(logy) matrix
     * @param nnode_y   number of nodes in y, or equivalently number of nodal rows
     */
    int i = node % (nnode_y);
    return i;
}
///@brief obtains the j-index of a node, based on the position of the node's number in the [nnode_y x nnode_x] NodeTopo matrix
int getInd_j(int node, int nnode_y){
    /**
     * The j index is based on a row major storage convention in the matrix of which the index will be used to extract from
     * e.g. A[rowlength * i + j]
     * @param node      the number of the node in the NodeTopo(logy) matrix
     * @param nnode_x   number of nodes in x, or equivalently number of nodal columns
     */
    int j = (node - (node % (nnode_y))) / (nnode_y);
    return j;
}
///@brief return the coordinate of a numbered node based on the coordinate matrix entered
double getCoord(int nodeID, double* X, int nnode_x, int nnode_y){
    /**
     * @param nodeID        The node number, based on NodeTopo matrix
     * @param X             Coordinate matrix (Usually X or Y)
     * @param nnode_x       Number of nodes in x
     * @param nnode_y       Number of nodes in y
     * @return val          The X-coordinate of the nodeID
     */
    double val = X[nnode_x * getInd_i(nodeID,nnode_y) + getInd_j(nodeID,nnode_y)];
    return val;
}



///@brief calculates the determinant of a 2x2 matrix
double detMatrix2(double* A){
    /**
     * @param A         matrix to calculate det(A) for
     * @return detJ     the determinant of =A
     */
    double detJ;
    detJ = A[3] * A[0] - A[2] * A[1];
    return detJ;
}

///@brief calculates the inverse of a 2x2 matrix
void invMatrix2(double* J, double* invJ){
    /**
     * @param J         the 2x2 matrix stored as an array
     * @param invJ      the location to store the inverse 2x2 matrix in
     */
    double detJ = detMatrix2(J);
    invJ[0] = J[3] / detJ;
    invJ[1] = -J[1] / detJ;
    invJ[2] = -J[2] / detJ;
    invJ[3] = J[0] / detJ;
}


