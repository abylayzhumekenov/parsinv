#ifndef PARSINV_VECTOR_H
#define PARSINV_VECTOR_H


#include <petsc.h>


/**
 * @brief Load a vector
 * 
 * @param comm 
 * @param type 
 * @param n 
 * @param N 
 * @param filename 
 * @param x 
 */
void ParsinvVecLoad(MPI_Comm comm, MatType type, int n, int N, const char filename[], Vec* x);


/**
 * @brief Remove NaN and Inf values from a vector and create a mask for valid elements
 * 
 * @param x 
 * @param y 
 * @param m 
 */
void ParsinvVecRemoveNA(Vec x, Vec* y, int* m);


#endif