.SUFFIXES: .C

CC=cc
C_FLAGS=-g

CPP=/usr/bin/g++
#CPP_FLAGS=-g -Wall -pedantic
CPP_FLAGS=-O3 -Wall -pedantic

.c.o:	$< 
	$(CC) $(C_FLAGS)  -c  $< 

.C.o:	$< 
	$(CPP) $(CPP_FLAGS) -c  $< 

OBJ = CrackFittingModule.o MakeCrackModule.o Standard.o EuropaCycloids.o ReadInTable.o \
	StressCalcModuleVisco.o Parser.o BFGSOptimization.o MatrixAlgebra.o CycloidParameters.o

ec:	$(OBJ)
	$(CPP) $(CPP_FLAGS) $(OBJ) -o $@

clean:
	- rm -f ec *.o a.out

depend :
	perl ./listdepend.scr dependencies.txt

include dependencies.txt     
