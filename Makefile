EXECUTABLE=twofluid
SOURCES=init.c alloc.c main.c readinputs.c utils.c output.c boundary.c restart.c algo_driver.c implicit.c
HEADER=edisk.h rk45.h

LDFLAGS=-llapack -lblas -lm -lgomp

CFLAGS=-c -fopenmp -Wall -O3 


BIN=bin/
SRC=src/
IN=inputs/
PY=src/pyutils/

CC=gcc-4.9

#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
OBJECTS=$(SOURCES:.c=.o)
CSOURCES=$(addprefix $(SRC),$(SOURCES))
COBJECTS=$(addprefix $(BIN),$(OBJECTS))
CHEADER=$(addprefix $(SRC),$(HEADER))




all: $(CSOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(COBJECTS) 
	$(CC)  $(COBJECTS) $(LDFLAGS) -o $@

$(BIN)%.o: $(SRC)%.c $(CHEADER) 
	$(CC) $(CFLAGS) $< -o $@

	
clean:
	rm $(COBJECTS) $(EXECUTABLE) 
