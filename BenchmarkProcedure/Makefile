# This makefile builds a generic tester for all heaps
CC 				 = g++
CPPFLAGS        += -O3 -Wall -g -DFOO=BAR

# paths to other files (starting from current directory)
BH				 = ../BinaryHeap
POLY			 = ../Polynomial
IO				 = ../IOHandling


all: benchmark


# build your application
benchmark: Options.o Polynomial.o IOHandling.o BinaryHeap.o BenchmarkProcedure.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

BenchmarkProcedure.o: BenchmarkProcedure.cpp $(POLY)/Polynomial.h $(IO)/IOHandling.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c BenchmarkProcedure.cpp

Options.o: ../Options/Options.h ../Options/Options.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c ../Options/Options.cpp

IOHandling.o: $(IO)/IOHandling.h $(IO)/IOHandling.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(IO)/IOHandling.cpp

Polynomial.o: $(POLY)/Polynomial.h $(POLY)/Polynomial.cpp $(POLY)/Monomial.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(POLY)/Polynomial.cpp

BinaryHeap.o: $(BH)/BinaryHeap.cpp $(BH)/BinaryHeap.h $(POLY)/Polynomial.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(BH)/BinaryHeap.cpp


clean:
	rm -f *.o benchmark

