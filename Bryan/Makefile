# Makefile for bryan - library maintenance
# Author: M. Jarrell
# Please not that the compiler flags for optimization will
# have to be changed for machines other than SUNS.
.SUFFIXES: .a .o .for .h
#
# Define a few things
#
F77=g77
LINK=g77
FLAGS=-O
BLAS=-L/usr/local/Linux_PIIISSE1/lib -lf77blas -latlas


RM = rm -f

OBJS = control.o error.o init.o main.o\
	output.o  update.o  

bryan: $(OBJS)
	$(LINK) -o bryan_linux *.o $(BLAS)

control.o: defs.h control.for
	$(F77) $(FLAGS) -c control.for

error.o: defs.h error.for
	$(F77) $(FLAGS) -c error.for

init.o: defs.h init.for
	$(F77) $(FLAGS) -c init.for

main.o: defs.h main.for
	$(F77) $(FLAGS) -c main.for

output.o: defs.h output.for
	$(F77) $(FLAGS) -c output.for

update.o: defs.h update.for
	$(F77) $(FLAGS) -c update.for

clean:
	$(RM) *.o bryan

