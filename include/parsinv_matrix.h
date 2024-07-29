#ifndef PARSINV_MATRIX_H
#define PARSINV_MATRIX_H


#include <petsc.h>
#include <petscmat.h>


/**
 * @brief Read a matrix header
 * 
 * @param filename 
 * @param m 
 * @param n 
 */
void ParsinvMatHeader(const char filename[], int* m, int* n);


/**
 * @brief Load a matrix
 * 
 * @param comm 
 * @param type 
 * @param m 
 * @param n 
 * @param M 
 * @param N 
 * @param filename 
 * @param A 
 */
void ParsinvMatLoad(MPI_Comm comm, MatType type, int m, int n, int M, int N, const char filename[], Mat* A);


/**
 * @brief Create a diagonal matrix from a vector
 * 
 * @param x 
 * @param A 
 */
void ParsinvMatDiagonal(Vec x, Mat* A);


/**
 * @brief Sequential Kronecker product
 * 
 * @param A 
 * @param B 
 * @param C 
 */
void ParsinvMatKroneckerSeq(Mat A, Mat B, Mat* C);


/**
 * @brief Parallel Kronecker product
 * 
 * @param A 
 * @param B 
 * @param C 
 */
void ParsinvMatKroneckerMPI(Mat A, Mat B, Mat* C);


/**
 * @brief A rank-1 update of a matrix
 * 
 * @param A 
 * @param alpha 
 * @param x 
 * @param y 
 */
void ParsinvMatUpdate(Mat A, double alpha, Vec x, Vec y);


/**
 * @brief Hadamard (elementwise) product of two matrices
 * 
 * @param A 
 * @param B 
 * @param C 
 */
void ParsinvMatHadamardSparse(Mat A, Mat B, Mat C);


/**
 * @brief Compute x^T*A*y
 * 
 * @param A 
 * @param x 
 * @param y 
 * @param work 
 * @param val 
 */
void ParsinvVecMatVec(Mat A, Vec x, Vec y, Vec work, double* val);


/**
 * @brief Compute the Schur complement S22 = A22 - A21^T * A11^-1 * A12
 * 
 * @param ksp 
 * @param A11 
 * @param A12 
 * @param A22 
 * @param S12 
 * @param S22 
 */
void ParsinvMatSchur(KSP ksp, Mat A11, Mat A12, Mat A22, Mat S12, Mat* S22);


/**
 * @brief Invert a dense matrix
 * 
 * @param A 
 * @param B 
 */
void ParsinvMatSolveDense(Mat A, Mat* B);


#endif