
CC=mpicc
CFLAGS=-O2 -g -Wall
FILES=3d-stencil.c

all:3d-stencil

3d-stencil:$(FILES)
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: clean

clean:
	@rm 3d-stencil
