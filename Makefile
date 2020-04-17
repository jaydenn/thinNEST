CC = g++
FLAGS = -g -O2 -std=c++11 -Wall -fopenmp
LIBS = -L../nest/build -lgsl -lgslcblas -lstdc++ -lm -lNEST_core
INCLUDE = -I./inc -I../nest/include/NEST -I../nest/include/Detectors -I./detectors
OBJECTS = src/TestSpectra.o src/migdalRate.o

default: thinNEST

thinNEST: $(OBJECTS) src/thinNEST.o
	$(CC) $(FLAGS) $^ -o $@ $(INCLUDE) $(LIBS)

src/%.o: src/%.cpp
	$(CC) $< $(FLAGS) $(INCLUDE) -c -o $@

clean:
	-rm src/*.o
	-rm -f ./thinNEST
