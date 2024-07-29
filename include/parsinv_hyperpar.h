#ifndef PARSINV_HYPERPAR_H
#define PARSINV_HYPERPAR_H


/**
 * @brief PARSINV_MANIFOLD_R2 = 0, PARSINV_MANIFOLD_S2 = 1
 * 
 */
typedef enum {
    PARSINV_MANIFOLD_R2 = 0,
    PARSINV_MANIFOLD_S2 = 1
} ParsinvManifold;


/**
 * @brief Compute C1 constant
 * 
 * @param d 
 * @param alpha 
 * @return double 
 */
double ParsinvHyperparC1(int d, double alpha);


/**
 * @brief Compute C2 constant when space is R2
 * 
 * @param d 
 * @param alpha 
 * @param gamma_s 
 * @return double 
 */
double ParsinvHyperparC2(int d, double alpha, double gamma_s);


/**
 * @brief Compute C2 constant when space is S2
 * 
 * @param alpha 
 * @param gamma_s 
 * @return double 
 */
double ParsinvHyperparC2Sphere(double alpha, double gamma_s);


/**
 * @brief Convert internal theta hyperparameters to gamma
 * 
 * @param theta [lrange_s, lrange_t, lsigma, lprec]
 * @param gamma [gamma_s, gamma_t, gamma_e, prec]
 * @param manifold_s R2 or S2
 */
void ParsinvHyperparSet(double* theta, double* gamma, ParsinvManifold manifold_s);


/**
 * @brief Evaluate the PC prior for the hyperparameters
 * 
 * @param theta 
 * @param manifold 
 * @param val 
 */
void ParsinvHyperparPrior(double* theta, ParsinvManifold manifold, double* val);


#endif