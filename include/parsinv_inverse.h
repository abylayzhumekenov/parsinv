#ifndef PARSINV_INVERSE_H
#define PARSINV_INVERSE_H


#include <petsc.h>
#include "parsinv_random.h"


void ParsinvInverseKSPCreate(MPI_Comm comm, Mat A, int n_over, KSP* ksp);


void ParsinvInverseKSPSetUp(KSP ksp);


void ParsinvInverseISCreate(KSP ksp, IS* is_sub, IS* is_over);


void ParsinvInverseMatCreate(KSP ksp, Mat* B);


void ParsinvInverseMatSubmatrix(KSP ksp, Mat* B);


void ParsinvInverseMatInvert(KSP ksp, Mat B);


void ParsinvInverseMatCorrect(KSP ksp, IS is_sub, Mat B, int n_samples, ParsinvRandom* rng);


void ParsinvInverseMatMatTrace(Mat A, Mat B, Mat C, IS is_over, double* trace);


void ParsinvInverseMatSolve(KSP ksp_postr, 
                            Mat Qyy, Mat Qub_postr, Mat Auy, Mat Ayb, Mat Cbb_postr, 
                            Vec y, Vec wy, Vec wu, Vec wu2, Vec wb, Vec wb2, Vec xu, Vec xb);


#endif