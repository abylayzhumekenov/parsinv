#include "parsinv_hyperpar.h"
#include <math.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


double ParsinvHyperparC1(int d, double alpha){

    return tgamma(alpha - d/2.0) / tgamma(alpha) / pow(4*M_PI, d/2.0);
}


double ParsinvHyperparC2(int d, double alpha, double gamma_s){

    return ParsinvHyperparC1(d, alpha) / pow(gamma_s, 2*alpha-d);
}


double ParsinvHyperparC2Sphere(double alpha, double gamma_s){

    if(gamma_s < 0.5){
        return 1 / (4*M_PI * pow(gamma_s, 2*alpha));
    } else if(gamma_s > 2.0){
        return 1 / (4*M_PI * (alpha-1) * pow(gamma_s, 2*alpha-2));
    } else {
        double c_sphere = 0;
        for(int k=0; k<50; k++){
            c_sphere += (2*k+1) / (4*M_PI * pow(gamma_s*gamma_s + k*(k+1), alpha));
        }
        return c_sphere;
    }
}


void ParsinvHyperparSet(double* theta, double* gamma, ParsinvManifold manifold_s){
    double c1, c2;
    double alphas[3] = {1, 2, 1};
    double smooth[3];
    smooth[0] = alphas[2] + alphas[1] * (alphas[0] - 0.5);
    smooth[2] = smooth[0] - 1.0;
    smooth[1] = fmin(alphas[0] - 0.5, smooth[2] / alphas[1]);

    c1 = ParsinvHyperparC1(1, 1);
    gamma[0] = sqrt(8 * smooth[2]) / exp(theta[0]);
    gamma[1] = exp(theta[1]) * pow(gamma[0], alphas[1]) / sqrt(8 * (alphas[0]-0.5));
    if(manifold_s == PARSINV_MANIFOLD_R2) c2 = ParsinvHyperparC2(2,    smooth[0], gamma[0]);      // if R2
    else                                  c2 = ParsinvHyperparC2Sphere(smooth[0], gamma[0]);      // if S2
    gamma[2] = sqrt(c1 * c2 / gamma[1] / pow(exp(theta[2]), 2));
    gamma[3] = exp(theta[3]);
}


void ParsinvHyperparPrior(double* theta, ParsinvManifold manifold, double* val){

    *val = 0.0;
    double gamma [4];
    double lambda[4];
    double alphas[3] = {1, 2, 1};
    double smooth[3];
    double prob  [4] = {0.01, 0.01, 0.01, 0.01};
    double quant [4] = {1.00, 1.00, 1.00, 1.00};

    ParsinvHyperparSet(theta, gamma, manifold);
    smooth[0] = alphas[2] + alphas[1] * (alphas[0] - 0.5);
    smooth[2] = smooth[0] - 1.0;
    smooth[1] = fmin(alphas[0] - 0.5, smooth[2] / alphas[1]);

    lambda[0] = -log(prob[0]) * pow(quant[0] / sqrt(8.0 * smooth[1]), 2.0 / 2.0);
    lambda[1] = -log(prob[1]) * pow(quant[1] / sqrt(8.0 * smooth[2]), 1.0 / 2.0);
    lambda[2] = -log(prob[2]) * pow(gamma[0], -smooth[1]) 
              * sqrt(tgamma(smooth[1]) / tgamma(smooth[1] + 2.0 / 2.0) / pow(4 * M_PI, 2.0 / 2.0));
    lambda[3] = -log(prob[3]) / quant[3];

    *val += log(2.0 / 2.0 * lambda[0]) - 2.0 / 2.0 * theta[0] - lambda[0] * exp(-2.0 / 2.0 * theta[0]);
    *val += log(1.0 / 2.0 * lambda[1]) - 1.0 / 2.0 * theta[1] - lambda[1] * exp(-1.0 / 2.0 * theta[1]);
    *val += log(lambda[2]) + theta[2] - lambda[2] * exp(theta[2]);
    *val += log(1.0 / 2.0 * lambda[3]) - 1.0 / 2.0 * theta[3] - lambda[3] * exp(-1.0 / 2.0 * theta[3]);
}