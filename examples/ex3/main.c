#include "../../include/parsinv.h"


int main(int argc, char** argv){

    int n_over = 1;

    for(int i=0; i<(argc); i++){
        if(!strcmp(argv[i], "-no")) n_over = atoi(argv[i+1]);
        argv[i] = 0;
    }

    int rank, size;
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, NULL, NULL);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    int n_sub;
    Mat Q, *QQ, CC, FF;
    KSP ksp, *ksp_sub;
    PC pc, pc_sub;
    IS *is_sub, *is_inner, *is_outer;
    PetscViewer viewer;

    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "Q", FILE_MODE_READ, &viewer);
    MatCreate(PETSC_COMM_WORLD, &Q);
    MatLoad(Q, viewer);
    // MatZeroEntries(Q);
    MatShift(Q, 1.0);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, Q, Q);
    KSPGetPC(ksp, &pc);
    PCSetOperators(pc, Q, Q);
    PCSetType(pc, PCGASM);
    PCGASMSetOverlap(pc, n_over);
    PCGASMCreateSubdomains(Q, n_over, &n_sub, &is_sub);
    PCSetUp(pc);
    PCGASMGetSubdomains(pc, &n_sub, &is_inner, &is_outer);
    PCGASMGetSubmatrices(pc, &n_sub, &QQ);
    MatDuplicate(*QQ, MAT_COPY_VALUES, &CC);

    PCGASMGetSubKSP(pc, &n_sub, NULL, &ksp_sub);
    KSPGetPC(ksp_sub[0], &pc_sub);
    PCSetType(pc_sub, PCCHOLESKY);
    PCFactorSetMatSolverType(pc_sub, MATSOLVERMUMPS);
    PCFactorSetReuseOrdering(pc_sub, PETSC_TRUE);
    PCFactorGetMatrix(pc_sub, &FF);

    KSPSetUp(ksp);          // needed for KSPSolve?
    PCSetUp(pc);
    KSPSetUp(ksp_sub[0]);   // not needed?
    PCSetUp(pc_sub);

    if(!rank) ISView(*is_inner, PETSC_VIEWER_STDOUT_SELF);
    if(!rank) ISView(*is_outer, PETSC_VIEWER_STDOUT_SELF);

    for(int i=0; i<3; i++){
        MatScale(Q, 0.1);
        
        // /* This does not work!!! */
        // PCSetUp(pc_sub);
        // PCSetUp(pc);

        /* This works!!! */
        PCSetUp(pc);
        PCSetUp(pc_sub);
        
        MatMumpsGetInverseTranspose(FF, CC);
        MatMumpsSetIcntl(FF, 30, 0);    // needed for KSPSolve
        MatView(*QQ, PETSC_VIEWER_STDOUT_WORLD);
        MatView(CC, PETSC_VIEWER_STDOUT_WORLD);
        PetscPrintf(PETSC_COMM_WORLD, "\n");
    }

    PetscFinalize();
    MPI_Finalize();

    return 0;
}