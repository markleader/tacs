include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = plate.o

default: ${OBJS}
	${CXX} -o plate plate.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o plate

test: default
	mpirun -np 1 plate
	mpirun -np 2 plate
	mpirun -np 4 plate