#include "parsinv_stmatrix.h"
#include "parsinv_hyperpar.h"
#include <math.h>


void ParsinvAssembleQyy(Mat Myy, double* theta, Mat Qyy){
    
    double tauy = exp(theta[3]);
    MatZeroEntries(Qyy);
    MatAXPY(Qyy, tauy, Myy, SAME_NONZERO_PATTERN);
}


void ParsinvAssembleQuu_prior(Mat* Muu, double* theta, int manifold_s, Mat Quu_prior){

    double gamma[4];
    ParsinvHyperparSet(theta, gamma, manifold_s);

    MatZeroEntries(Quu_prior);

    MatAXPY(Quu_prior, pow(gamma[2], 2) * pow(gamma[1], 0) * pow(gamma[0], 6) * 1, Muu[0], SUBSET_NONZERO_PATTERN);
    MatAXPY(Quu_prior, pow(gamma[2], 2) * pow(gamma[1], 0) * pow(gamma[0], 4) * 3, Muu[1], SUBSET_NONZERO_PATTERN);
    MatAXPY(Quu_prior, pow(gamma[2], 2) * pow(gamma[1], 0) * pow(gamma[0], 2) * 3, Muu[2], SUBSET_NONZERO_PATTERN);
    MatAXPY(Quu_prior, pow(gamma[2], 2) * pow(gamma[1], 0) * pow(gamma[0], 0) * 1, Muu[3], SUBSET_NONZERO_PATTERN);

    MatAXPY(Quu_prior, pow(gamma[2], 2) * pow(gamma[1], 1) * pow(gamma[0], 4) * 1, Muu[4], SUBSET_NONZERO_PATTERN);
    MatAXPY(Quu_prior, pow(gamma[2], 2) * pow(gamma[1], 1) * pow(gamma[0], 2) * 2, Muu[5], SUBSET_NONZERO_PATTERN);
    MatAXPY(Quu_prior, pow(gamma[2], 2) * pow(gamma[1], 1) * pow(gamma[0], 0) * 1, Muu[6], SUBSET_NONZERO_PATTERN);

    MatAXPY(Quu_prior, pow(gamma[2], 2) * pow(gamma[1], 2) * pow(gamma[0], 2) * 1, Muu[7], SUBSET_NONZERO_PATTERN);
    MatAXPY(Quu_prior, pow(gamma[2], 2) * pow(gamma[1], 2) * pow(gamma[0], 0) * 1, Muu[8], SUBSET_NONZERO_PATTERN);
}


void ParsinvAssembleQuu_postr(Mat* Muu, double* theta, int manifold_s, Mat Quu_postr){

    double gamma[4];
    ParsinvHyperparSet(theta, gamma, manifold_s);

    ParsinvAssembleQuu_prior(Muu, theta, manifold_s, Quu_postr);
    MatAXPY(Quu_postr, gamma[3], Muu[9], SUBSET_NONZERO_PATTERN);
}


void ParsinvAssembleQub_prior(Mat Qub_prior){

    MatZeroEntries(Qub_prior);
}


void ParsinvAssembleQub_postr(Mat Mub, double* theta, Mat Qub_postr){

    double tauy = exp(theta[3]);
    ParsinvAssembleQub_prior(Qub_postr);
    MatAXPY(Qub_postr, tauy, Mub, SAME_NONZERO_PATTERN);
}


void ParsinvAssembleQbb_prior(Mat Qbb_prior){

    MatZeroEntries(Qbb_prior);
    MatShift(Qbb_prior, 1e-3);
}


void ParsinvAssembleQbb_postr(Mat Mbb, double* theta, Mat Qbb_postr){

    double tauy = exp(theta[3]);
    ParsinvAssembleQbb_prior(Qbb_postr);
    MatAXPY(Qbb_postr, tauy, Mbb, SAME_NONZERO_PATTERN);
}
