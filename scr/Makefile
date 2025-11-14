FC=gfortran
MKLROOT=/home/cerda/.soft/mkl_2024/mkl/2024.0
FFLAGS=-O2 -fopenmp -fdefault-integer-8 -I"${MKLROOT}/include"
LINK=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group ${MKLROOT}/../../compiler/2024.0/lib/libiomp5.so -lpthread -lm -ldl
#LINK=-llapack -lblas 
OUT=DIAB_2D

OBJ=method.o block.o fchk_reader.o mat_reader.o general.o set_var.o

all: ${OBJ} main.f90  
	${FC} ${FFLAGS} ${OBJ} main.f90 ${LINK} -o ${OUT}
	rm *.o
fchk_reader.o: fchk_reader.f90  
	${FC} ${FFLAGS} -c fchk_reader.f90
mat_reader.o: mat_reader.f90  
	${FC} ${FFLAGS} -c mat_reader.f90
set_var.o: set_var.f90  
	${FC} ${FFLAGS} -c set_var.f90
general.o: general.f90  
	${FC} ${FFLAGS} -c general.f90
block.o: block.f90  
	${FC} ${FFLAGS} -c block.f90
method.o: method.f90  
	${FC} ${FFLAGS} -c method.f90
