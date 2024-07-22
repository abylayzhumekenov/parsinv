#ifndef PARSINV_SOLVER_H
#define PARSINV_SOLVER_H


#include <petsc.h>


void ParsinvInverseMUMPS();
// {

    // int n;
    // double *b_array;
    // Mat L;
    // MatInfo info;
    // IS isrow, iscol;

    // MatGetSize(A, &n, NULL);
    // MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);
    // MatGetOrdering(A, MATORDERINGNATURAL, &isrow, &iscol);
    // MatGetFactor(A, MATSOLVERMUMPS, MAT_FACTOR_CHOLESKY, &L);
    // MatCholeskyFactorSymbolic(L, A, isrow, NULL);
    // MatCholeskyFactorNumeric(L, A, NULL);

    // if(reuse == MAT_INITIAL_MATRIX) MatDuplicate(A, MAT_SHARE_NONZERO_PATTERN, B);
    // MatGetInfo(A, MAT_LOCAL, &info);
    // MatSeqAIJGetArray(*B, &b_array);
    // for(int i=0; i<info.nz_allocated; i++) b_array[i] = 1.0;
    // MatSeqAIJRestoreArray(*B, &b_array);
    // MatMumpsGetInverseTranspose(L, *B);

    // MatDestroy(&L);
    // ISDestroy(&isrow);
    // ISDestroy(&iscol);
// }


#endif