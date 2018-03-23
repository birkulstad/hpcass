#ifndef HPCASS_OPER_H
#define HPCASS_OPER_H

#define F77NAME(x) x##_
extern "C" {
// LAPACK routine for solving systems of linear equations
void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                    const int& lda, int * ipiv, double * B,
                    const int& ldb, int& info);
}



double detMatrix2(double*);
void invMatrix2(double*, double*);

int getInd_i(int, int);
int getInd_j(int, int);

double getCoord(int, double*, int, int);


#endif //HPCASS_OPER_H
