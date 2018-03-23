#ifndef HPCASS_IO_H
#define HPCASS_IO_H

void WriteVtkFile(int, int, int, double* , int* , double* );
void errorMessage(int);

void printMatrix(double*, int, int);
void printMatrix(int*, int, int);

#endif //HPCASS_IO_H
