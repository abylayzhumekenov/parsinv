#ifndef PARSINV_RANDOM_H
#define PARSINV_RANDOM_H


#include <petsc.h>
#include "parsinv_pcg.h"


/**
 * @brief Random number generator (RNG)
 * 
 */
typedef pcg64_random_t ParsinvRandom;


/**
 * @brief Type of RNG. 
 * 
 */
typedef enum {
    PARSINV_RANDOM_DEFAULT = 0,
    PARSINV_RANDOM_UNIQUE  = 1,
    PARSINV_RANDOM_USER    = 2
} ParsinvRandomType;


/**
 * @brief Create and seed an RNG
 * 
 * DEFAULT and UNIQUE types are similar and ignore initial state and sequence.
 * For reproducibility choose USER and fix state and sequence. 
 * Remember to set different sequences to different MPI processes, otherwise random numbers are correlated.
 * 
 * @param type 
 * @param initstate 
 * @param initseq 
 * @return ParsinvRandom 
 */
ParsinvRandom ParsinvRandomCreate(ParsinvRandomType type, pcg128_t initstate, pcg128_t initseq);


/**
 * @brief Generate a standard uniform random number
 * 
 * @param rng 
 * @return double 
 */
double ParsinvRandomUniform(ParsinvRandom* rng);


/**
 * @brief Generate a standard normal random number
 * 
 * @param rng 
 * @return double 
 */
double ParsinvRandomNormal(ParsinvRandom* rng);


/**
 * @brief Generate a standard normal random vector
 * 
 * @param rng 
 * @param z 
 */
void ParsinvRandomNormalIID(ParsinvRandom* rng, Vec z);


#endif