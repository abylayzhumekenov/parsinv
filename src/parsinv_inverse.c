#include "parsinv_inverse.h"
#include "parsinv_matrix.h"
#include "parsinv_random.h"


void ParsinvInverseKSPCreate(MPI_Comm comm, Mat A, int n_over, KSP* ksp){

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int n_sub;
    PC pc, pc_sub;
    KSP *ksp_sub;
    IS *is_sub;

    KSPCreate(comm, ksp);
    KSPSetOperators(*ksp, A, A);
    PCCreate(comm, &pc);
    KSPSetPC(*ksp, pc);
    PCSetOperators(pc, A, A);
    PCSetType(pc, PCGASM);
    PCGASMSetOverlap(pc, n_over);
    PCGASMCreateSubdomains(A, size, &n_sub, &is_sub);
    PCSetUp(pc);
    
    PCGASMGetSubKSP(pc, &n_sub, NULL, &ksp_sub);
    KSPGetPC(ksp_sub[0], &pc_sub);
    PCSetType(pc_sub, PCCHOLESKY);
    PCFactorSetMatSolverType(pc_sub, MATSOLVERMUMPS);
    PCFactorSetReuseOrdering(pc_sub, PETSC_TRUE);
    PCSetUp(pc_sub);
    KSPSetUp(ksp_sub[0]);
    PCSetUp(pc);
    KSPSetUp(*ksp);
}


void ParsinvInverseKSPSetUp(KSP ksp){

    int n_sub;
    KSP* ksp_sub;
    PC pc, pc_sub;
    
    KSPGetPC(ksp, &pc);
    PCGASMGetSubKSP(pc, &n_sub, NULL, &ksp_sub);
    KSPGetPC(ksp_sub[0], &pc_sub);
    KSPSetUp(ksp);
    PCSetUp(pc);
    KSPSetUp(ksp_sub[0]);
    PCSetUp(pc_sub);
}


void ParsinvInverseISCreate(KSP ksp, IS* is_sub, IS* is_over){

    PC pc;
    int n_sub, n_outer;
    const int *is_outer_array;
    IS *is_inner, *is_outer, is_overlap;
    ISLocalToGlobalMapping l2g;

    KSPGetPC(ksp, &pc);
    PCGASMGetSubdomains(pc, &n_sub, &is_inner, &is_outer);
    ISGetSize(*is_outer, &n_outer);
    ISGetIndices(*is_outer, &is_outer_array);
    ISCreateGeneral(PETSC_COMM_WORLD, n_outer, is_outer_array, PETSC_COPY_VALUES, is_sub);
    ISRestoreIndices(*is_outer, &is_outer_array);

    ISLocalToGlobalMappingCreateIS(*is_outer, &l2g);
    ISDifference(*is_outer, *is_inner, &is_overlap);
    ISGlobalToLocalMappingApplyIS(l2g, IS_GTOLM_DROP, is_overlap, is_over);
    ISLocalToGlobalMappingDestroy(&l2g);
    ISDestroy(&is_overlap);
}


void ParsinvInverseMatCreate(KSP ksp, Mat* B){

    int n_sub;
    PC pc;
    Mat *A_sub;

    KSPGetPC(ksp, &pc);
    PCGASMGetSubmatrices(pc, &n_sub, &A_sub);
    MatDuplicate(*A_sub, MAT_COPY_VALUES, B);
}


void ParsinvInverseMatSubmatrix(KSP ksp, Mat* B){

    int n_sub;
    PC pc;
    Mat* C;

    KSPGetPC(ksp, &pc);
    PCGASMGetSubmatrices(pc, &n_sub, &C);
    *B = C[0];
}


void ParsinvInverseMatInvert(KSP ksp, Mat B){

    int n_sub;
    Mat F;
    KSP *ksp_sub;
    PC pc, pc_sub;
    
    KSPGetPC(ksp, &pc);
    PCGASMGetSubKSP(pc, &n_sub, NULL, &ksp_sub);
    KSPGetPC(ksp_sub[0], &pc_sub);
    PCFactorGetMatrix(pc_sub, &F);
    MatMumpsGetInverseTranspose(F, B);
    MatMumpsSetIcntl(F, 30, 0);
}


void ParsinvInverseMatCorrect(KSP ksp, IS is_sub, Mat B, int n_samples, ParsinvRandom* rng){

    int n_sub;
    PC pc, pc_sub;
    KSP *ksp_sub;
    Mat A, F;
    Vec z, z_sub, z_outer, zz_outer;
    Vec x, x_sub, x_outer;
    Vec w, w_sub, w_outer, ww_outer;

    KSPGetOperators(ksp, &A, NULL);
    KSPGetPC(ksp, &pc);
    PCGASMGetSubKSP(pc, &n_sub, NULL, &ksp_sub);
    KSPGetPC(ksp_sub[0], &pc_sub);
    PCFactorGetMatrix(pc_sub, &F);

    MatZeroEntries(B);
    MatCreateVecs(A, &z, &w);
    VecDuplicate(w, &x);

    VecGetSubVector(z, is_sub, &z_sub);
    VecGetSubVector(x, is_sub, &x_sub);
    VecGetSubVector(w, is_sub, &w_sub);
    VecCreateLocalVector(z_sub, &z_outer);
    VecCreateLocalVector(x_sub, &x_outer);
    VecCreateLocalVector(w_sub, &w_outer);
    VecRestoreSubVector(z, is_sub, &z_sub);
    VecRestoreSubVector(x, is_sub, &x_sub);
    VecRestoreSubVector(w, is_sub, &w_sub);

    VecDuplicate(z_outer, &zz_outer);
    VecDuplicate(w_outer, &ww_outer);

    for(int i=0; i<n_samples; i++){

        ParsinvRandomNormalIID(rng, z);
        MatMult(A, z, w);
        KSPSolve(ksp, z, x);
            
        VecGetSubVector(z, is_sub, &z_sub);
        VecGetSubVector(x, is_sub, &x_sub);
        VecGetSubVector(w, is_sub, &w_sub);
        VecGetLocalVector(z_sub, z_outer);
        VecGetLocalVector(x_sub, x_outer);
        VecGetLocalVector(w_sub, w_outer);
        
        MatSolve(F, w_outer, ww_outer);
        MatSolve(F, z_outer, zz_outer);
        VecAXPY(ww_outer, -1.0, z_outer);
        VecAXPY(zz_outer, -1.0, x_outer);
        ParsinvMatUpdate(B, 1.0/n_samples, ww_outer, zz_outer);

        VecRestoreLocalVector(z_sub, z_outer);
        VecRestoreLocalVector(x_sub, x_outer);
        VecRestoreLocalVector(w_sub, w_outer);
        VecRestoreSubVector(z, is_sub, &z_sub);
        VecRestoreSubVector(x, is_sub, &x_sub);
        VecRestoreSubVector(w, is_sub, &w_sub);
    }

    VecDestroy(&z_outer);
    VecDestroy(&x_outer);
    VecDestroy(&w_outer);
    VecDestroy(&zz_outer);
    VecDestroy(&ww_outer);
    VecDestroy(&z);
    VecDestroy(&x);
    VecDestroy(&w);

    // ---------------------------------------------------------------

    /* Only the correction part!... */

    /* RBMC estimator */
    // C_aa = Q_aa^-1 + (Q_aa^-1 Q_as x_s) * (...)^T
    //      = Q_aa^-1 + (Q_aa^-1 (Qx)_a - x_a) * (...)^T
    //      = Q_aa^-1 + (Q_aa^-1 w_a - x_a) * (...)^T                           where w = Qx

    /* RBCC (RBH) estimator */
    // C_aa = Q_aa^-1 + (Q_aa^-1 Q_as z_s) * (Q_aa^-1 Q_as (Q^-1z)_s)^T
    //      = Q_aa^-1 + (Q_aa^-1 (Qz)_a - z_a) * (Q_aa^-1 z_a - (Q^-1z)_a)^T
    //      = Q_aa^-1 + (Q_aa^-1 w_a - z_a) * (Q_aa^-1 z_a - y_a)^T             where w = Qz, y = Q^-1z

    // ---------------------------------------------------------------
}


void ParsinvInverseMatMatTrace(Mat A, Mat B, Mat C, IS is_over, double* trace){

    Vec c;
    int ok = 0;
    double sum = 0;

    ParsinvMatHadamardSparse(A, B, C);
    MatCreateVecs(C, &c, NULL);
    MatGetRowSum(C, c);
    VecISSet(c, is_over, 0.0);
    VecSum(c, &sum);
    VecDestroy(&c);
    ok += MPI_Allreduce(&sum, trace, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
}


void ParsinvInverseMatSolve(KSP ksp_postr, 
                            Mat Qyy, Mat Qub_postr, Mat Auy, Mat Ayb, Mat Cbb_postr, 
                            Vec y, Vec wy, Vec wu, Vec wu2, Vec wb, Vec wb2, Vec xu, Vec xb){

    MatMult(Qyy, y, wy);
    MatMult(Auy, wy, wu);
    MatMultTranspose(Ayb, wy, wb);
    KSPSolve(ksp_postr, wu, wu2);
    MatMultTranspose(Qub_postr, wu2, wb2);
    VecAXPY(wb, -1.0, wb2);
    MatMult(Cbb_postr, wb, xb);

    MatMult(Qub_postr, xb, wu);
    KSPSolve(ksp_postr, wu, xu);
    VecScale(xu, -1.0);
    VecAXPY(xu, 1.0, wu2);
}