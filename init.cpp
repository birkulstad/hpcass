#include "init.h"
#include <valarray>


/* Functions that are associated with initialisation and elementwise operations are stored in this file
 */

///@brief Populates an array of doubles with a constant value
void populate(double* a, int length, double val){
    /**
     * @param a         The array to populate
     * @param length    Length of the array to populate
     * @param val       The value to populate with
     */
    for (int i = 0; i < length; i++){
        a[i] = val;
    }
}

///@brief Populates an array of doubles with a constant value
void populate(int* a, int length, int val){
    /**
     * @param a         The array to populate
     * @param length    Length of the array to populate
     * @param val       The value to populate with
     */
    for (int i = 0; i < length; i++){
        a[i] = val;
    }
}

///@brief writes the node numbers of the nodes along an edge selected for a boundary condition into BCNodes
void getBCNodes(int edge, int nnode_x, int nnode_y, int& nBCNodes, int* BCNodes, int* NodeTopo){
    /**
     * @param edge      The edge selected for the boundary condition
     * @param nnode_x   Number of nodes in the x
     * @param nnode_y   Number of nodes in the y
     * @param nBCNodes  Number of nodes along the selected edge
     * @param BCNodes   The nodes along the selected edge
     * @param NodeTopo  The Node Topology matrix
     */
    int mult_i, mult_j; // Defines multipliers
    edgeswitch(edge, nnode_x, nnode_y, nBCNodes, mult_i, mult_j); // Obtains nBCNodes and multipliers
    for (int i = 0; i < nBCNodes; i++){
        BCNodes[i] = NodeTopo[mult_i * i + mult_j];
    }
}

///@brief returns the number of nodes along the selected edge together with multipliers mult_i and mult_j
void edgeswitch(int edge, int nnode_x, int nnode_y,int& nBCNodes, int& mult_i, int& mult_j){
    switch (edge){
        /**
         * @param edge          The edge selected for which to return number of nodes and multipliers
         * @param nnode_x       Number of nodes in x of the mesh
         * @param nnode_y       Number of nodes in y of the mesh
         * @param nBCNodes      Returning the number of nodes along the selected edge
         * @param mult_i        i-multiplier to determine node numbers along the edge
         * @param mult_j        j-multiplier to determine node numbers along the edge
         */
        // Left edge
        case 0: nBCNodes = nnode_y;
            mult_i = nnode_x;
            mult_j = 0;
            break;

            // Top edge
        case 1: nBCNodes = nnode_x;
            mult_i = 1;
            mult_j = nnode_x * (nnode_y - 1);
            break;

            // Right edge
        case 2: nBCNodes = nnode_y;
            mult_i = nnode_x;
            mult_j = (nnode_x - 1);
            break;

            // Bottom edge
        case 3: nBCNodes = nnode_x;
            mult_i = 1;
            mult_j = 0;
            break;

        default:
            break;
    }
}

///@brief returns an array with nnode_x coordinates distributed over the length (end-begin)
void linspace(double begin, double end, int nnode_x, double* arr){
    /**
     * @param begin     Beginning of the linear space
     * @param end       End of the linear space
     * @param nnode_x   Discretisation of the linear space
     * @param arr       Arrayu to store the coordinates in
     */
    for (int i = 0; i < nnode_x; i++){
        arr[i] = (end - begin) / (nnode_x - 1) * i;
    }
}