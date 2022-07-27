NESTDIR=/home/unimelb.edu.au/jnewstead/Codes/nest
CC = g++
FLAGS = -g -O2 -std=c++11 -Wall
LIBS = -lgsl -lgslcblas -lstdc++ -lm -lNESTCore
INCLUDE = -I./inc -I${NESTDIR}/include/NEST -I${NESTDIR}/include/Detectors -I./detectors -I../gcem/include
OBJECTS = src/TestSpectra.o src/Target.o src/PhysicalConstants.o

default: thinNEST

thinNEST: $(OBJECTS) src/thinNEST.o
	$(CC) $(FLAGS) $^ -o $@ $(INCLUDE) $(LIBS)

src/%.o: src/%.cpp
	$(CC) $< $(FLAGS) $(INCLUDE) -c -o $@

clean:
	-rm src/*.o
	-rm -f ./thinNEST
