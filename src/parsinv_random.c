#include "parsinv_random.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


void ParsinvRandomCreate(ParsinvRandomType type, pcg128_t initstate, pcg128_t initseq, ParsinvRandom* rng){

    if(type == PARSINV_RANDOM_DEFAULT)  pcg64_srandom_r(rng, 0,          (intptr_t)&rng);
    if(type == PARSINV_RANDOM_UNIQUE)   pcg64_srandom_r(rng, time(NULL), (intptr_t)&rng);
    if(type == PARSINV_RANDOM_USER)     pcg64_srandom_r(rng, initstate,  initseq);
}


double ParsinvRandomUniform(ParsinvRandom* rng){

    return ldexp(pcg64_random_r(rng), -64);
}


double ParsinvRandomNormal(ParsinvRandom* rng){

    double u1 = ParsinvRandomUniform(rng);
    double u2 = ParsinvRandomUniform(rng);

    return sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
}


void ParsinvRandomNormalIID(ParsinvRandom* rng, Vec z){

    int istart, iend;
    VecGetOwnershipRange(z, &istart, &iend);
    for(int i=istart; i<iend; i++) VecSetValue(z, i, ParsinvRandomNormal(rng), INSERT_VALUES);
    VecAssemblyBegin(z);
    VecAssemblyEnd(z);
}


