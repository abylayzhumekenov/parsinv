#ifndef PARSINV_STMATRIX_H
#define PARSINV_STMATRIX_H


#include <petsc.h>


void ParsinvAssembleQyy(Mat Myy, double* theta, Mat Qyy);


void ParsinvAssembleQuu_prior(Mat* M, double* theta, int manifold_s, Mat Q);


void ParsinvAssembleQuu_postr(Mat* M, double* theta, int manifold_s, Mat Q);


void ParsinvAssembleQub_prior(Mat Qub_prior);


void ParsinvAssembleQub_postr(Mat Mub, double* theta, Mat Qub_postr);


void ParsinvAssembleQbb_prior(Mat Qbb_prior);


void ParsinvAssembleQbb_postr(Mat Mbb, double* theta, Mat Qbb_postr);


#endif