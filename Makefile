CC = g++
FLAGS = -g -O2 -std=c++11 -Wall
LIBS = -L../nest_v2/build -lgsl -lgslcblas -lstdc++ -lm -lNESTcore
INCLUDE = -I./inc -I../nest_v2/include/NEST -I../nest_v2/include/Detectors -I./detectors
OBJECTS = src/TestSpectra.o src/migdalRate.o

default: thinNEST

thinNEST: $(OBJECTS) src/thinNEST.o
	$(CC) $(FLAGS) $^ -o $@ $(INCLUDE) $(LIBS)

src/%.o: src/%.cpp
	$(CC) $< $(FLAGS) $(INCLUDE) -c -o $@

clean:
	-rm src/*.o
	-rm -f ./thinNEST
