CC=g++
CREATEEXE=$(CC) -o $@ $^


all: main 

main: Options.o Chain.o BinaryHeapWithChaining.o PolynomialArithmetic.o main.cpp 
	$(CREATEEXE)


Chain.o: ../BinaryHeapWithChaining/Chain.cpp ../BinaryHeapWithChaining/Chain.h
	$(CC) $(CFLAGS) -c ../BinaryHeapWithChaining/Chain.cpp

BinaryHeapWithChaining.o: ../BinaryHeapWithChaining/BinaryHeapWithChaining.cpp ../BinaryHeapWithChaining/BinaryHeapWithChaining.h ../../Options/Options.h
	$(CC) $(CFLAGS) -c ../BinaryHeapWithChaining/BinaryHeapWithChaining.cpp
	
PolynomialArithmetic.o: ../PolynomialArithmetic/PolynomialArithmetic.cpp ../PolynomialArithmetic/PolynomialArithmetic.h ../../Options/Options.h
	$(CC) $(CFLAGS) -c ../PolynomialArithmetic/PolynomialArithmetic.cpp
	
Options.o: ../../Options/Options.h ../../Options/Options.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c ../../Options/Options.cpp

clean:
	rm -f *.o main


