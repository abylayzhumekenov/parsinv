#include "parsinv_vector.h"


void ParsinvVecLoad(MPI_Comm comm, MatType type, int n, int N, const char filename[], Vec* x){

    PetscViewer viewer;
    PetscViewerBinaryOpen(comm, filename, FILE_MODE_READ, &viewer);
    VecCreate(comm, x);
    VecSetType(*x, type);
    if((n!=PETSC_DECIDE) || (N!=PETSC_DECIDE)) VecSetSizes(*x, n, N);
    VecLoad(*x, viewer);
    PetscViewerDestroy(&viewer);
}


void ParsinvVecRemoveNA(Vec x, Vec* y, int* m){

    int n, mm = 0, ok = 0;
    double *x_array, *y_array;
    VecDuplicate(x, y);
    VecSet(*y, 1.0);
    VecGetLocalSize(x, &n);
    VecGetArray(x, &x_array);
    VecGetArray(*y, &y_array);
    for(int i=0; i<n; i++){
        if(isnan(x_array[i]) || isinf(x_array[i])){
            y_array[i] = 0.0;
            x_array[i] = 0.0;
        } else {
            mm++;
        }
    }
    VecRestoreArray(*y, &y_array);
    VecRestoreArray(x, &x_array);
    ok += MPI_Allreduce(&mm, m, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
}
