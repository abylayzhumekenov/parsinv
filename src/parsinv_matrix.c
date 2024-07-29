#include "parsinv_matrix.h"
#include <petscmat.h>
#include <petscblaslapack.h>


void ParsinvMatHeader(const char filename[], int* m, int* n){

    int count, header[4];
    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_READ, &viewer);
    PetscViewerBinaryRead(viewer, header, 4, &count, PETSC_INT);
    *m = header[1];
    *n = header[2];
    PetscViewerDestroy(&viewer);
}


void ParsinvMatLoad(MPI_Comm comm, MatType type, int m, int n, int M, int N, const char filename[], Mat* A){

    PetscViewer viewer;
    PetscViewerBinaryOpen(comm, filename, FILE_MODE_READ, &viewer);
    MatCreate(comm, A);
    MatSetType(*A, type);
    if((m!=PETSC_DECIDE) || (n!=PETSC_DECIDE) || (M!=PETSC_DECIDE) || (N!=PETSC_DECIDE)){
        MatSetSizes(*A, m, n, M, N);
    }
    MatLoad(*A, viewer);
    PetscViewerDestroy(&viewer);
}


void ParsinvMatDiagonal(Vec x, Mat* A){

    MPI_Comm comm;
    int n, istart, iend, j = 0;
    double* x_array;
    PetscObjectGetComm((PetscObject)x, &comm);
    VecGetLocalSize(x, &n);
    MatCreateAIJ(comm, n, n, PETSC_DETERMINE, PETSC_DETERMINE, 1, NULL, 0, NULL, A);
    MatGetOwnershipRange(*A, &istart, &iend);
    VecGetArray(x, &x_array);
    for(int i=istart; i<iend; i++){
        MatSetValue(*A, i, i, x_array[j], INSERT_VALUES);
        j++;
    }
    VecRestoreArray(x, &x_array);
    MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
}


void ParsinvMatKroneckerSeq(Mat A, Mat B, Mat* C){

    MatSeqAIJKron(A, B, MAT_INITIAL_MATRIX, C);
}


void ParsinvMatKroneckerMPI(Mat A, Mat B, Mat* C){

    int n, m;
    Mat AA, CC;
    
    MatAIJGetLocalMat(A, &AA);
    ParsinvMatKroneckerSeq(AA, B, &CC);
    MatGetSize(AA, &n, NULL);
    MatGetSize(B, NULL, &m);
    MatCreateMPIMatConcatenateSeqMat(PETSC_COMM_WORLD, CC, n*m, MAT_INITIAL_MATRIX, C);

    MatDestroy(&AA);
    MatDestroy(&CC);
}


void ParsinvMatUpdate(Mat A, double alpha, Vec x, Vec y){

    int n;
    const int *ia, *ja;
    double *a_array, *x_array, *y_array;
    PetscBool done;

    MatSeqAIJGetArray(A, &a_array);
    MatGetRowIJ(A, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done);
    VecGetArray(x, &x_array);
    VecGetArray(y, &y_array);

    for(int i=0; i<n; i++){
        for(int ii=ia[i]; ii<ia[i+1]; ii++){
            int j = ja[ii];
            a_array[ii] += alpha * x_array[i] * y_array[j];
        }
    }

    MatRestoreRowIJ(A, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done);
    MatSeqAIJRestoreArray(A, &a_array);
    VecRestoreArray(x, &x_array);
    VecRestoreArray(y, &y_array);
}


void ParsinvMatHadamardSparse(Mat A, Mat B, Mat C){

    MatInfo info;
    double *a_array, *b_array, *c_array;
    MatGetInfo(A, MAT_LOCAL, &info);
    MatSeqAIJGetArray(A, &a_array);
    MatSeqAIJGetArray(B, &b_array);
    MatSeqAIJGetArray(C, &c_array);
    for(int i=0; i<(int)info.nz_allocated; i++) c_array[i] = a_array[i] * b_array[i];
    MatSeqAIJRestoreArray(A, &a_array);
    MatSeqAIJRestoreArray(B, &b_array);
    MatSeqAIJRestoreArray(C, &c_array);
}


void ParsinvVecMatVec(Mat A, Vec x, Vec y, Vec work, double* val){

    MatMult(A, y, work);
    VecDot(x, work, val);
}


void ParsinvMatSchur(KSP ksp, Mat A11, Mat A12, Mat A22, Mat S12, Mat* S22){

    KSPMatSolve(ksp, A12, S12);
    MatTransposeMatMult(A12, S12, MAT_REUSE_MATRIX, PETSC_DEFAULT, S22);
    MatAYPX(*S22, -1.0, A22, SAME_NONZERO_PATTERN);
}


void ParsinvMatSolveDense(Mat A, Mat* B){

    int rank, n;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MatGetLocalSize(A, &n, NULL);
    MatCopy(A, *B, SAME_NONZERO_PATTERN);
    
    int i_array[n], info;
    double *b_array, work[n];
    MatDenseGetArray(*B, &b_array);
    if(!rank){
        dgetrf_(&n, &n, b_array, &n, i_array, &info);
        dgetri_(&n, b_array, &n, i_array, work, &n, &info);
    }
    MatDenseRestoreArray(*B, &b_array);
}