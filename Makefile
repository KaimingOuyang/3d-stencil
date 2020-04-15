
DEBUG=0
CC=mpicc
CFLAGS=-g -Wall
BIN=3d-stencil

ifeq ($(DEBUG), 0)
	CFLAGS+= -O2 
else
	CFLAGS+= -O
endif

ifeq ($(CC), ampicc)
	CFLAGS += -balancer GreedyRefineLB -DAMPI_LOAD_BALANCE
else
	CFLAGS += -fpie -pie -rdynamic -pthread
endif

FILES=3d-stencil.c

all:$(BIN)

$(BIN):$(FILES)
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: clean

clean:
	@rm 3d-stencil
