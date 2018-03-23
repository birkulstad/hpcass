#ifndef HPCASS_INIT_H
#define HPCASS_INIT_H

void populate(double*, int, double);
void populate(int*, int, int);

void getBCNodes(int, int, int, int&, int*, int*);
void edgeswitch(int, int, int, int&, int&, int&);

void linspace(double, double, int, double*);

#endif //HPCASS_INIT_H
