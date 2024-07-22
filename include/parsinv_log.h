#ifndef PARSINV_LOG_H
#define PARSINV_LOG_H


#include <stdlib.h>
#include <stdio.h>
#include <petsc.h>


/**
 * @brief Macro for error message
 * 
 */
#define PARSINV_ERROR_MSG(code, msg) \
        printf("%s:%s:%i: %s\n", __FILE__, __FUNCTION__, __LINE__, msg); \
        exit(code)


/**
 * @brief Print formatted message (debug only)
 * 
 * @param comm 
 * @param msg 
 * @param ... 
 */
void ParsinvLog(MPI_Comm comm, const char msg[], ...);


/**
 * @brief Profile checkpoint (blocking, debug only)
 * 
 * @param comm 
 * @param time 
 * @param memory 
 */
void ParsinvCheckpoint(MPI_Comm comm, double* time, double* memory);


#endif