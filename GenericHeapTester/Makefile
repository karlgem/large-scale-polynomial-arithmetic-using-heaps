# This makefile builds a generic tester for all heaps


all: main
# path to stxxl.mk file
# from your stxxl installation
include  $(HOME)/stxxl/trunk/stxxl.mk

# use the variables from stxxl.mk
CXX              = $(STXXL_CXX)
CPPFLAGS        += $(STXXL_CPPFLAGS)
LDLIBS          += $(STXXL_LDLIBS)

# add your own optimization, warning, debug, ... flags
# (these are *not* set in stxxl.mk)
CPPFLAGS        += -O3 -Wall -g -DFOO=BAR -fno-inline

# paths to other files (starting from current directory)
FH				 = ../FunnelHeap
FHM				 = ../FunnelHeapWithMerging
FHC				 = ../FunnelHeapComponents

BH				 = ../BinaryHeap
BHC				 = ../BinaryHeapWithChaining

POLY			 = ../Polynomial

# build your application
# (my_example.o is generated from my_example.cpp automatically)
main: Options.o Buffer.o Link.o KMerger.o FunnelHeap.o FunnelHeapWithMerging.o BinaryHeap.o BinaryHeapWithChaining.o GenericHeapTester.o 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

GenericHeapTester.o: GenericHeapTester.cpp $(POLY)/Monomial.h $(FHM)/FunnelHeapWithMerging.h $(FH)/FunnelHeap.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c GenericHeapTester.cpp

Options.o: ../Options/Options.h ../Options/Options.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c ../Options/Options.cpp

BinaryHeap.o: $(BH)/BinaryHeap.cpp $(BH)/BinaryHeap.h $(POLY)/Polynomial.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(BH)/BinaryHeap.cpp

BinaryHeapWithChaining.o: $(BHC)/BinaryHeapWithChaining.cpp $(BHC)/BinaryHeapWithChaining.h $(POLY)/Polynomial.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(BHC)/BinaryHeapWithChaining.cpp

# New Funnel heap with merging classes

Buffer.o: $(FHC)/Buffer.cpp $(FHC)/Buffer.h $(POLY)/Monomial.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(FHC)/Buffer.cpp

KMerger.o: $(FHC)/KMerger.cpp $(FHC)/KMerger.h 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(FHC)/KMerger.cpp

Link.o: $(FHC)/Link.cpp $(FHC)/Link.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(FHC)/Link.cpp

FunnelHeap.o: $(FH)/FunnelHeap.cpp $(FH)/FunnelHeap.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(FH)/FunnelHeap.cpp

FunnelHeapWithMerging.o: $(FHM)/FunnelHeapWithMerging.cpp $(FHM)/FunnelHeapWithMerging.h $(POLY)/Polynomial.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(FHM)/FunnelHeapWithMerging.cpp




clean:
	rm -f *.o main

