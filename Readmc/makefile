# Makefile for dqa - library maintenance
# Author: M. Jarrell

.SUFFIXES: .a .o .for .h

RM = rm -f

OBJS = readmc_linux.o

readmc: $(OBJS)
	f77 -O -o readmc_linux *.o -lblas 

readmc_linux.o: readmc_linux.for
	f77 -O -c readmc_linux.for

clean:
	$(RM) *.o readmc_linux

