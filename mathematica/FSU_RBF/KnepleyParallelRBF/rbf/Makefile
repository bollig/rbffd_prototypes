CFLAGS           = -g
FFLAGS           = 
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = 
MANSEC           = 
EXAMPLESC        = main.c 
EXAMPLESF        = 
EXAMPLESCH       = 
EXAMPLESFH       = 
EXAMPLESMATLAB   = 
DIRS             = . 

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

all: main

main: main.o
	${CLINKER} -o $@ $< ${PETSC_SNES_LIB}

test: test.o
	${CLINKER} -o $@ $< ${PETSC_SNES_LIB}

