#include "parsinv_log.h"


void ParsinvLog(MPI_Comm comm, const char msg[], ...){

    #ifdef DEBUG
    int rank;
    MPI_Comm_rank(comm, &rank);

    va_list arg_list;
    va_start(arg_list, msg);
    if(!rank) {
        vprintf(msg, arg_list);
    }
    va_end(arg_list);
    #endif
}


void ParsinvCheckpoint(MPI_Comm comm, double* time, double* memory){

    #ifdef DEBUG
    double time_old = *time, memory_old = *memory;

    MPI_Barrier(comm);
    PetscTime(time);
    PetscMallocGetCurrentUsage(memory);

    ParsinvLog(comm, "Time:\t%f sec\n\tMemory:\t%i bytes\n\n", *time-time_old, (int)(*memory-memory_old));
    #endif
}