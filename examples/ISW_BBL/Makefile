SUNTANSHOME=../../main
include $(SUNTANSHOME)/Makefile.in

ifneq ($(MPIHOME),)
  CC = $(MPIHOME)/bin/mpicc
  MPIDEF = 
  MPIINC = -I$(MPIHOME)
else
  CC = gcc
  MPIDEF = -DNOMPI
  MPIINC = 
endif

ifneq ($(PARMETISHOME),)
  PARMETISINC = -I$(PARMETISHOME)/ParMETISLib
endif

LD = $(CC) 
CFLAGS = -O2
MATHLIB = -lm

EXEC = bl_test
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC)
DEFS = $(MPIDEF)
NUMPROCS = 8

all:	data

test:	data
	sh $(EXEC).sh $(NUMPROCS)

data: $(SUN)

.c.o:	
	$(LD) -c $(CFLAGS) $(INCLUDES) $(DEFS) $*.c

$(SUN):	initialization.o boundaries.o state.o uservertcoordinate.o sources.o 
	cp initialization.o boundaries.o state.o uservertcoordinate.o sources.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

clean:
	rm -f *.o

clobber: clean
	rm -rf *~ */*~ \#*\# PI* $(EXEC) gmon.out data rundata/*~
