#include "../../include/parsinv.h"
#include <math.h>


/**
 * @brief Run example 1
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */
int main(int argc, char** argv){

    ParsinvRandom rng;
    PetscViewer viewer;
    Mat J0, J1, J2;
    Mat G0, G1, G2, G3;
    Mat At, As, Auy, Ayb;
    Mat Muu[10], Myy, Mub, Mbb;
    Vec y, obs;

    Mat Qyy, Qyy_, Cyy;
    Mat Quu_prior, Qub_prior, Qbb_prior, Sub_prior, Sbb_prior, Cbb_prior, Quu_prior_sub, Cuu_prior_sub;
    Mat Quu_postr, Qub_postr, Qbb_postr, Sub_postr, Sbb_postr, Cbb_postr, Quu_postr_sub, Cuu_postr_sub;
    Mat Quu_prior_, Qub_prior_, Qbb_prior_, Sub_prior_, Sbb_prior_, Cbb_prior_, Quu_prior_sub_, Cuu_prior_sub_;
    Mat Quu_postr_, Qub_postr_, Qbb_postr_, Sub_postr_, Sbb_postr_, Cbb_postr_, Quu_postr_sub_, Cuu_postr_sub_;
    Mat Wuu, Wub, Wbb, Wuu_prior_sub, Wuu_postr_sub, Wuu_prior_sub_, Wuu_postr_sub_;
    Vec xu, xb, xu_, xb_;
    Vec wu, wu2, wb, wb2, wy;
    KSP ksp_prior, ksp_postr, ksp_prior_, ksp_postr_;
    IS  is_sub, is_over;

    int nu_t, nu_s, nu, nb;
    int ny_t, ny_s, ny, ny_o;
    int nu_t_local, nu_local;
    int ny_t_local, ny_local;

    ParsinvManifold manifold = PARSINV_MANIFOLD_S2;
    int n_over = 1;
    int n_iter = 100;
    int n_samples = 10;
    int n_restart = 1000;
    int use_gd = 0;
    double lrate0   = 0.5;
    double lrate    = lrate0;
    double drate    = 1.0;
    double dpoly    = 0.0;
    double epsilon  = 5e-3;
    double norm0    = 1.0;
    double normg    = 0.0;
    double rtol     = 1e-3;
    double atol     = 1e-3;
    double theta[4] = {0, 0, 0, 0};
    double grad [4] = {0, 0, 0, 0};
    double hess [4] = {0, 0, 0, 0};
    double work[16];
    double time = 0, memory = 0;

    /* Set parameters from the options */
    for(int i=0; i<(argc); i++){
        if(!strcmp(argv[i], "-ni"))  n_iter     = atoi(argv[i+1]);
        if(!strcmp(argv[i], "-ns"))  n_samples  = atoi(argv[i+1]);
        if(!strcmp(argv[i], "-no"))  n_over     = atoi(argv[i+1]);
        if(!strcmp(argv[i], "-nr"))  n_restart  = atoi(argv[i+1]);
        if(!strcmp(argv[i], "-gd"))  use_gd     = atoi(argv[i+1]);
        if(!strcmp(argv[i], "-lr"))  lrate0     = atof(argv[i+1]);
        if(!strcmp(argv[i], "-dr"))  drate      = atof(argv[i+1]);
        if(!strcmp(argv[i], "-dp"))  dpoly      = atof(argv[i+1]);
        if(!strcmp(argv[i], "-ee"))  epsilon    = atof(argv[i+1]);
        if(!strcmp(argv[i], "-rt"))  rtol       = atof(argv[i+1]);
        if(!strcmp(argv[i], "-at"))  atol       = atof(argv[i+1]);
        if(!strcmp(argv[i], "-hh")){
            for(int j=0; j<4; j++) theta[j] = atof(argv[i+1+j]);
        }
        argv[i] = 0;
    }

    int rank, size;
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, NULL, NULL);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    ParsinvRandomCreate(PARSINV_RANDOM_UNIQUE, 0, 0, &rng);

    // -----------------------------------------------------------------------------------------------------------------

    ParsinvCheckpoint(PETSC_COMM_WORLD, &time, &memory);

    ParsinvMatHeader("data/At", &nu_t, &ny_t);
    ParsinvMatHeader("data/As", &nu_s, &ny_s);
    ParsinvMatHeader("data/Ab", &ny, &nb);
    ParsinvMatLoad(PETSC_COMM_WORLD, MATMPIAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, "data/J0", &J0);
    ParsinvMatLoad(PETSC_COMM_WORLD, MATMPIAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, "data/J1", &J1);
    ParsinvMatLoad(PETSC_COMM_WORLD, MATMPIAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, "data/J2", &J2);
    ParsinvMatLoad(PETSC_COMM_SELF,  MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, "data/G0", &G0);
    ParsinvMatLoad(PETSC_COMM_SELF,  MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, "data/G1", &G1);
    ParsinvMatLoad(PETSC_COMM_SELF,  MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, "data/G2", &G2);
    ParsinvMatLoad(PETSC_COMM_SELF,  MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, "data/G3", &G3);
    ParsinvMatLoad(PETSC_COMM_WORLD, MATMPIAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, "data/At", &At);
    ParsinvMatLoad(PETSC_COMM_SELF,  MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, "data/As", &As);

    MatGetLocalSize(At, &nu_t_local, &ny_t_local);  // At^T
    nu       = nu_t       * nu_s;
    nu_local = nu_t_local * nu_s;
    ny_local = ny_t_local * ny_s;
    ParsinvMatLoad(PETSC_COMM_WORLD, MATMPIDENSE, ny_local, nb*(!rank), ny, nb, "data/Ab", &Ayb);
    ParsinvVecLoad(PETSC_COMM_WORLD, VECMPI,      ny_local,             ny,     "data/y",  &y);
    ParsinvVecRemoveNA(y, &obs, &ny_o);

    ParsinvLog(PETSC_COMM_WORLD, "Global:\tnt=%i\tns=%i\tnu=%i\tny=%i\n", nu_t, nu_s, nu, ny);
    ParsinvLog(PETSC_COMM_WORLD, "Local: \tnt=%i\tns=%i\tnu=%i\tny=%i\n", nu_t_local, nu_s, nu_local, ny_local);

    ParsinvMatKroneckerMPI(J0, G0, &Muu[0]);
    ParsinvMatKroneckerMPI(J0, G1, &Muu[1]);
    ParsinvMatKroneckerMPI(J0, G2, &Muu[2]);
    ParsinvMatKroneckerMPI(J0, G3, &Muu[3]);
    ParsinvMatKroneckerMPI(J1, G0, &Muu[4]);
    ParsinvMatKroneckerMPI(J1, G1, &Muu[5]);
    ParsinvMatKroneckerMPI(J1, G2, &Muu[6]);
    ParsinvMatKroneckerMPI(J2, G0, &Muu[7]);
    ParsinvMatKroneckerMPI(J2, G1, &Muu[8]);
    ParsinvMatKroneckerMPI(At, As, &Auy);

    ParsinvMatDiagonal(obs, &Myy);
    MatDiagonalScale(Auy, NULL, obs);
    MatDiagonalScale(Ayb, obs, NULL);
    MatDiagonalScale(Myy, obs, NULL);
    MatDuplicate(Myy, MAT_COPY_VALUES, &Qyy);
    MatDuplicate(Myy, MAT_COPY_VALUES, &Qyy_);
    MatDuplicate(Myy, MAT_COPY_VALUES, &Cyy);
    
    MatMatTransposeMult(Auy, Auy, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Muu[9]);
    MatMatTransposeMult(Auy, Auy, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Wuu);
    MatMatMult         (Auy, Ayb, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Mub);
    MatMatMult         (Auy, Ayb, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Wub);
    MatTransposeMatMult(Ayb, Ayb, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Mbb);
    MatTransposeMatMult(Ayb, Ayb, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Wbb);
    
    MatDuplicate(Muu[9], MAT_COPY_VALUES, &Quu_prior);
    for(int i=0; i<9; i++) MatAXPY(Quu_prior, 1.0, Muu[i], DIFFERENT_NONZERO_PATTERN);
    MatDuplicate(Quu_prior, MAT_COPY_VALUES, &Quu_postr);
    MatDuplicate(Quu_prior, MAT_COPY_VALUES, &Quu_prior_);
    MatDuplicate(Quu_prior, MAT_COPY_VALUES, &Quu_postr_);
    
    MatMatMult(Auy, Ayb, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Qub_prior);
    MatDuplicate(Qub_prior, MAT_COPY_VALUES, &Qub_postr);
    MatDuplicate(Qub_prior, MAT_COPY_VALUES, &Qub_prior_);
    MatDuplicate(Qub_prior, MAT_COPY_VALUES, &Qub_postr_);
    MatDuplicate(Qub_prior, MAT_COPY_VALUES, &Sub_prior);
    MatDuplicate(Qub_prior, MAT_COPY_VALUES, &Sub_postr);
    MatDuplicate(Qub_prior, MAT_COPY_VALUES, &Sub_prior_);
    MatDuplicate(Qub_prior, MAT_COPY_VALUES, &Sub_postr_);

    MatTransposeMatMult(Ayb, Ayb, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Qbb_prior);
    MatDuplicate(Qbb_prior, MAT_COPY_VALUES, &Qbb_postr);
    MatDuplicate(Qbb_prior, MAT_COPY_VALUES, &Qbb_prior_);
    MatDuplicate(Qbb_prior, MAT_COPY_VALUES, &Qbb_postr_);
    MatDuplicate(Qbb_prior, MAT_COPY_VALUES, &Cbb_prior);
    MatDuplicate(Qbb_prior, MAT_COPY_VALUES, &Cbb_postr);
    MatDuplicate(Qbb_prior, MAT_COPY_VALUES, &Cbb_prior_);
    MatDuplicate(Qbb_prior, MAT_COPY_VALUES, &Cbb_postr_);
    MatDuplicate(Qbb_prior, MAT_COPY_VALUES, &Sbb_prior);
    MatDuplicate(Qbb_prior, MAT_COPY_VALUES, &Sbb_postr);
    MatDuplicate(Qbb_prior, MAT_COPY_VALUES, &Sbb_prior_);
    MatDuplicate(Qbb_prior, MAT_COPY_VALUES, &Sbb_postr_);
    
    ParsinvInverseKSPCreate(PETSC_COMM_WORLD, Quu_prior, n_over, &ksp_prior);
    ParsinvInverseKSPCreate(PETSC_COMM_WORLD, Quu_postr, n_over, &ksp_postr);
    ParsinvInverseKSPCreate(PETSC_COMM_WORLD, Quu_prior_, n_over, &ksp_prior_);
    ParsinvInverseKSPCreate(PETSC_COMM_WORLD, Quu_postr_, n_over, &ksp_postr_);
    ParsinvInverseKSPSetUp(ksp_prior);
    ParsinvInverseKSPSetUp(ksp_postr);
    ParsinvInverseKSPSetUp(ksp_prior_);
    ParsinvInverseKSPSetUp(ksp_postr_);
    ParsinvInverseISCreate(ksp_prior, &is_sub, &is_over);
    ParsinvInverseMatSubmatrix(ksp_prior, &Quu_prior_sub);
    ParsinvInverseMatSubmatrix(ksp_postr, &Quu_postr_sub);
    ParsinvInverseMatSubmatrix(ksp_prior_, &Quu_prior_sub_);
    ParsinvInverseMatSubmatrix(ksp_postr_, &Quu_postr_sub_);
    ParsinvInverseMatCreate(ksp_prior, &Cuu_prior_sub);
    ParsinvInverseMatCreate(ksp_postr, &Cuu_postr_sub);
    ParsinvInverseMatCreate(ksp_prior_, &Cuu_prior_sub_);
    ParsinvInverseMatCreate(ksp_postr_, &Cuu_postr_sub_);
    ParsinvInverseMatCreate(ksp_prior, &Wuu_prior_sub);
    ParsinvInverseMatCreate(ksp_postr, &Wuu_postr_sub);
    ParsinvInverseMatCreate(ksp_prior_, &Wuu_prior_sub_);
    ParsinvInverseMatCreate(ksp_postr_, &Wuu_postr_sub_);

    MatCreateVecs(Quu_prior, &xu, &xu_);
    MatCreateVecs(Quu_prior, &wu, &wu2);
    MatCreateVecs(Qbb_prior, &xb, &xb_);
    MatCreateVecs(Qbb_prior, &wb, &wb2);
    MatCreateVecs(Qyy, &wy, NULL);

    ParsinvCheckpoint(PETSC_COMM_WORLD, &time, &memory);
    ParsinvLog(PETSC_COMM_WORLD, "Starting optimizer with options:\n"
                                "-ni %i\n"
                                "-ns %i\n"
                                "-no %i\n"
                                "-nr %i\n"
                                "-gd %i\n"
                                "-lr %f\n"
                                "-dr %f\n"
                                "-dp %f\n"
                                "-ee %f\n"
                                "-rt %f\n"
                                "-at %f\n"
                                "-hh %f %f %f %f\n\n", 
        n_iter, n_samples, n_over, n_restart, use_gd, lrate0, drate, dpoly, epsilon, rtol, atol, theta[0], theta[1], theta[2], theta[3]);

    // ---------------------------------------------------------------------------------------

    for(int iter=0; iter<n_iter; iter++){

        ParsinvLog(PETSC_COMM_WORLD, "Iteration: %i\n", iter);

        /* LIKELIHOOD */
        ParsinvAssembleQyy(Myy, theta, Qyy);
        ParsinvVecMatVec(Qyy, y, y, wy, &work[0]);

        /* PRIOR */
        ParsinvAssembleQuu_prior(Muu, theta, manifold, Quu_prior);
        ParsinvAssembleQub_prior(Qub_prior);
        ParsinvAssembleQbb_prior(Qbb_prior);
        ParsinvInverseKSPSetUp(ksp_prior);
        ParsinvInverseMatInvert(ksp_prior, Cuu_prior_sub);
        ParsinvInverseMatCorrect(ksp_prior, is_sub, Wuu_prior_sub, n_samples, &rng);

        /* POSTERIOR */
        ParsinvAssembleQuu_postr(Muu, theta, manifold, Quu_postr);
        ParsinvAssembleQub_postr(Mub, theta, Qub_postr);
        ParsinvAssembleQbb_postr(Mbb, theta, Qbb_postr);
        ParsinvInverseKSPSetUp(ksp_postr);
        ParsinvInverseMatInvert(ksp_postr, Cuu_postr_sub);
        ParsinvInverseMatCorrect(ksp_postr, is_sub, Wuu_postr_sub, n_samples, &rng);
        ParsinvMatSchur(ksp_postr, Quu_postr, Qub_postr, Qbb_postr, Sub_postr, &Sbb_postr);
        ParsinvMatSolveDense(Sbb_postr, &Cbb_postr);
        ParsinvInverseMatSolve(ksp_postr, Qyy, Qub_postr, Auy, Ayb, Cbb_postr, 
                               y, wy, wu, wu2, wb, wb2, xu, xb);
        ParsinvVecMatVec(Quu_postr, xu, xu, wu, &work[1]);
        ParsinvVecMatVec(Qub_postr, xu, xb, wu, &work[2]);
        ParsinvVecMatVec(Qbb_postr, xb, xb, wb, &work[3]);

        /* HYPERPRIOR */
        ParsinvHyperparPrior(theta, manifold, &work[4]);
        
        for(int k=0; k<4; k++){

            theta[k] += epsilon;

            /* LIKELIHOOD */
            ParsinvAssembleQyy(Myy, theta, Qyy_);
            ParsinvVecMatVec(Qyy_, y, y, wy, &work[5]);

            /* PRIOR */
            ParsinvAssembleQuu_prior(Muu, theta, manifold, Quu_prior_);
            ParsinvAssembleQub_prior(Qub_prior_);
            ParsinvAssembleQbb_prior(Qbb_prior_);
            ParsinvInverseKSPSetUp(ksp_prior_);
            ParsinvInverseMatMatTraceSparseSeq(Cuu_prior_sub, Quu_prior_sub_, Cuu_prior_sub_, is_over, &work[6]);
            ParsinvInverseMatMatTraceSparseSeq(Wuu_prior_sub, Quu_prior_sub_, Cuu_prior_sub_, is_over, &work[14]);

            /* POSTERIOR */
            ParsinvAssembleQuu_postr(Muu, theta, manifold, Quu_postr_);
            ParsinvAssembleQub_postr(Mub, theta, Qub_postr_);
            ParsinvAssembleQbb_postr(Mbb, theta, Qbb_postr_);
            ParsinvInverseKSPSetUp(ksp_postr_);
            ParsinvMatSchur(ksp_postr_, Quu_postr_, Qub_postr_, Qbb_postr_, Sub_postr_, &Sbb_postr_);
            ParsinvMatSolveDense(Sbb_postr_, &Cbb_postr_);
            ParsinvInverseMatSolve(ksp_postr_, Qyy_, Qub_postr_, Auy, Ayb, Cbb_postr_, 
                                   y, wy, wu, wu2, wb, wb2, xu_, xb_);
            ParsinvVecMatVec(Quu_postr_, xu_, xu_, wu, &work[7]);
            ParsinvVecMatVec(Qub_postr_, xu_, xb_, wu, &work[8]);
            ParsinvVecMatVec(Qbb_postr_, xb_, xb_, wb, &work[9]);
            ParsinvInverseMatMatTraceSparseSeq(Cuu_postr_sub, Quu_postr_sub_, Cuu_postr_sub_, is_over, &work[10]);
            ParsinvInverseMatMatTraceSparseSeq(Wuu_postr_sub, Quu_postr_sub_, Cuu_postr_sub_, is_over, &work[15]);
            
            /* HYPERPRIOR */
            ParsinvHyperparPrior(theta, manifold, &work[11]);

            grad[k] = (((ny*epsilon*(k==3) - work[5] + work[0]) +                                           // likelihood
                        (work[6] - nu) -                                                                    // prior
                        (work[10] - nu - work[7]-2*work[8]-work[9] + work[1]+2*work[2]+work[3])) / 2.0 +    // posterior
                        (work[11] - work[4])) / epsilon;                                                    // hypeprior
            
            // ----------------------------------------------------------------------------------------------------------

            theta[k] -= 2*epsilon;

            /* LIKELIHOOD */
            ParsinvAssembleQyy(Myy, theta, Qyy_);
            ParsinvVecMatVec(Qyy_, y, y, wy, &work[5]);

            /* PRIOR */
            ParsinvAssembleQuu_prior(Muu, theta, manifold, Quu_prior_);
            ParsinvAssembleQub_prior(Qub_prior_);
            ParsinvAssembleQbb_prior(Qbb_prior_);
            ParsinvInverseKSPSetUp(ksp_prior_);
            ParsinvInverseMatInvert(ksp_prior_, Cuu_prior_sub_);
            ParsinvInverseMatMatTraceSparseSeq(Cuu_prior_sub_, Quu_prior_sub, Wuu_prior_sub_, is_over, &work[6]);

            /* POSTERIOR */
            ParsinvAssembleQuu_postr(Muu, theta, manifold, Quu_postr_);
            ParsinvAssembleQub_postr(Mub, theta, Qub_postr_);
            ParsinvAssembleQbb_postr(Mbb, theta, Qbb_postr_);
            ParsinvInverseKSPSetUp(ksp_postr_);
            ParsinvInverseMatInvert(ksp_postr_, Cuu_postr_sub_);
            ParsinvMatSchur(ksp_postr_, Quu_postr_, Qub_postr_, Qbb_postr_, Sub_postr_, &Sbb_postr_);
            ParsinvMatSolveDense(Sbb_postr_, &Cbb_postr_);
            ParsinvInverseMatSolve(ksp_postr_, Qyy_, Qub_postr_, Auy, Ayb, Cbb_postr_, 
                                   y, wy, wu, wu2, wb, wb2, xu_, xb_);
            ParsinvVecMatVec(Quu_postr_, xu_, xu_, wu, &work[7]);
            ParsinvVecMatVec(Qub_postr_, xu_, xb_, wu, &work[8]);
            ParsinvVecMatVec(Qbb_postr_, xb_, xb_, wb, &work[9]);
            ParsinvInverseMatMatTraceSparseSeq(Cuu_postr_sub_, Quu_postr_sub, Wuu_postr_sub_, is_over, &work[10]);
            
            /* HYPERPRIOR */
            ParsinvHyperparPrior(theta, manifold, &work[11]);

            hess[k] = (((ny*epsilon*(k==3) - work[0] + work[5]) +                                           // likelihood
                        (work[6] - nu) -                                                                    // prior
                        (work[10] - nu - work[1]-2*work[2]-work[3] + work[7]+2*work[8]+work[9])) / 2.0 +    // posterior
                        (work[4] - work[11])) / epsilon;                                                    // hypeprior
            hess[k] = (grad[k] - hess[k]) / epsilon;                        // form hessian from two gradients
            hess[k] = (hess[k] < -1.0) ? hess[k] : -1.0;                    // use gradient descent if negative hess is not PSD
            hess[k] = (!use_gd)        ? hess[k] : -1.0;                    // use gradient ascent if gd = 1
            grad[k] += (work[14] + work[15]) / 2.0 / epsilon;               // correction part
            hess[k] = (fabs(grad[k] / hess[k]) < 1.0) ? hess[k] : -ceil(fabs(grad[k]));   // limit the update magnitude
            
            theta[k] += epsilon;
        }

        normg = 0.0;
        for(int k=0; k<4; k++)              normg += grad[k] * grad[k];
        if(iter == 0)                       norm0 = normg;
        if(dpoly > 0.5)                     lrate = pow(iter%n_restart+1, -dpoly) * lrate0;
        else                                lrate = pow(drate, iter%n_restart) * lrate0;

        ParsinvLog(PETSC_COMM_WORLD, "theta\t\tgrad\t\thess\n");
        for(int k=0; k<4; k++)              ParsinvLog(PETSC_COMM_WORLD, "%f\t%f\t%f\n", theta[k], grad[k], hess[k]);
        ParsinvLog(PETSC_COMM_WORLD, "Abs |g|:\t%f\n", pow(normg, 0.5));
        ParsinvLog(PETSC_COMM_WORLD, "Rel |g|:\t%f\n", pow(normg / norm0, 0.5));
        ParsinvLog(PETSC_COMM_WORLD, "lrate:  \t%f\n", lrate);

        if(normg         < atol*atol){      ParsinvLog(PETSC_COMM_WORLD, "Converged!\n"); break; }
        if(normg / norm0 < rtol*rtol){      ParsinvLog(PETSC_COMM_WORLD, "Converged!\n"); break; }
        if(norm0 / normg < rtol*rtol){      ParsinvLog(PETSC_COMM_WORLD, "Diverged! \n"); break; }
        for(int k=0; k<4; k++)              theta[k] -= lrate * grad[k] / hess[k];

        ParsinvLog(PETSC_COMM_WORLD, "\n");
        ParsinvCheckpoint(PETSC_COMM_WORLD, &time, &memory);
    }


    /* Save results */
    Vec theta_vec, su, sb;
    MatCreateVecs(Quu_postr, &su, NULL);
    MatCreateVecs(Cbb_postr, &sb, NULL);
    
    if(!rank){
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, 4, theta, &theta_vec);
        PetscViewerBinaryOpen(PETSC_COMM_SELF, "R/data/theta", FILE_MODE_WRITE, &viewer);
        VecView(theta_vec, viewer);
    }

    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "R/data/muu", FILE_MODE_WRITE, &viewer);
    VecView(xu, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "R/data/mub", FILE_MODE_WRITE, &viewer);
    VecView(xb, viewer);
    
    MatGetDiagonal(Cbb_postr, sb);
    VecPow(sb, 0.5);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "R/data/sdb", FILE_MODE_WRITE, &viewer);
    VecView(sb, viewer);
    
    MatMatMult(Sub_postr, Cbb_postr, MAT_REUSE_MATRIX, PETSC_DEFAULT, &Wub);
    ParsinvMatHadamardDenseMPI(Sub_postr, Wub, Sub_postr_);
    ParsinvInverseMatMatDiagonalDenseMPI(Sub_postr, Wub, Sub_postr_, wu);
    MatAXPY(Cuu_postr_sub, 1.0, Wuu_postr_sub, SAME_NONZERO_PATTERN);
    ParsinvInverseMatDiagonalSparseSeq(Cuu_postr_sub, is_over, su);
    VecAXPY(su, 1.0, wu);
    VecPow(su, 0.5);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "R/data/sdu", FILE_MODE_WRITE, &viewer);
    VecView(su, viewer);

    /* Finalize */
    PetscFinalize();
    MPI_Finalize();

    return 0;
}