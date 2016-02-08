CC=mpic++
INCPATH=../include/
LIBPATH=../lib/
LIBS=-llapack -lblas -lf2c -lm 
EXE_MAIN=main
SRC=main.cpp matlap.cpp

all:
	mpic++ -I $(INCPATH) -L $(LIBPATH) -o $(EXE_MAIN) $(SRC) $(LIBS) -std=c++11
	mv $(EXE_MAIN) ../bin

