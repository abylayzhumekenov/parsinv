CC = mpicc
CFLAGS = -g -Wall -fopenmp

INC = ../../include
LIB = ../../lib

IFLAGS = -I$(INC) \
		 -I$(PETSC_DIR)/include \
		 -I$(PETSC_DIR)/$(PETSC_ARCH)/include 
LFLAGS = -L$(LIB) -lparsinv -lpetsc

all: main.c
	$(CC) $(CFLAGS) -o main main.c $(IFLAGS) -Wl,-rpath,$(LIB) $(LFLAGS) 