CC = mpicc
CFLAGS = -g -Wall -fopenmp

INC = ../../include
LIB = ../../lib

IFLAGS = -I$(INC) \
		 -I$(PETSC_DIR)/include \
		 -I$(PETSC_DIR)/$(PETSC_ARCH)/include 
LFLAGS = -L$(LIB) -lparsinv \
		 -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc

.PHONY: all clean

all: main.c
	$(CC) $(CFLAGS) -o main main.c $(IFLAGS) -Wl,-rpath,$(LIB) $(LFLAGS)

clean:
	rm -f main