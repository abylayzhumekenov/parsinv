CC = mpicc
CFLAGS = -g -Wall -fPIC -O1 -DDEBUG
TFLAGS = -shared

INC = include
SRC = src
OBJ = obj
LIB = lib

INCS = $(wildcard $(INC)/*.h)
SRCS = $(wildcard $(SRC)/*.c)
OBJS = $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SRCS))
DEPS = $(patsubst $(SRC)/%.c, $(OBJ)/%.d, $(SRCS))
TRGT = $(LIB)/libparsinv.so

IFLAGS = -I$(INC) \
		 -I$(PETSC_DIR)/include \
		 -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LFLAGS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
		 -lpetsc

.PHONY: all debug release clean

all: debug

debug: $(TRGT)

release: CFLAGS := $(filter-out -O1 -DDEBUG, $(CFLAGS)) -O3
release: clean $(TRGT)

$(TRGT): $(OBJS)
	$(CC) $(TFLAGS) -o $(TRGT) $(OBJS) $(LFLAGS)

$(OBJ)/%.o: $(SRC)/%.c Makefile
	$(CC) $(CFLAGS) $(IFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f $(OBJS) $(DEPS) $(TRGT)

-include $(DEPS)