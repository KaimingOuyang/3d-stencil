DEBUG=0
CC=mpicc
CFLAGS=-g -Wall
ifeq ($(DEBUG), 0)
	CFLAGS+=-O2
else
	CFLAGS+=-O
endif

FILES=3d-stencil.c

all:3d-stencil

3d-stencil:$(FILES)
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: clean

clean:
	@rm 3d-stencil
