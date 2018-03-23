#include "io.h"
#include "colormod.h" // namespace Color
#include <fstream>
#include <iomanip>
#include <iostream>


using namespace std;

Color::Modifier def(Color::FG_DEFAULT);
Color::Modifier red(Color::FG_RED);
Color::Modifier green(Color::FG_GREEN);


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

void errorMessage(int e){
    switch (e){
        default:
            cout << red << "ERROR! You probably did something bad" << def << endl;
            break;
        case (1):
            cout << red << "ERROR:  Real lengths must be positive!" << endl;
            cout << green << "Check command line input parameters 2, 3, 4, and 5." << def << endl;
            break;
        case (2):
            cout << red << "ERROR:  The mesh must be at least 2x2!" << endl;
            cout << green << "Check command line input parameters 6 and 7." << def << endl;
            break;
        case(3):
            cout << red << "ERROR:  The conductivity matrix D must be positive definite!" << endl;
            cout << green << "Check command line input parameters 8 and 9." << def << endl;
            cout << green << "D = [kxx, kxy; kxy, kyy]" << def << endl;
            break;
        case(4):
            cout << red << "Error! One or more of the BC edges are invalid" << endl;
            cout << green << "Check command line input parameters 11 and 12." << def << endl;
            cout << green << "0 = left, 1 = top, 2 = right, 3 = bottom" << def << endl;
            break;
        case(5):
            cout << red << "Error! The BCs can not be enforced at the same edge" << endl;
            cout << green << "Check command line input parameters 11 and 12." << def << endl;
            cout << green << "0 = left, 1 = top, 2 = right, 3 = bottom" << def << endl;
            break;
        case(6):
            cout << red << "Error! The dimensions do not produce a continuous object!" << endl;
            cout << green << "Check command line input parameters 1, 2, 3 and 4." << def << endl;
            cout << green << "Should have: a < 2 * (h1 + h2)/(L * (2 - L))" << def << endl;
            break;
    }
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