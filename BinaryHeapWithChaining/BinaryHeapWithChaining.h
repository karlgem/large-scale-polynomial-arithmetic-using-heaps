/*
 *  BinaryHeapWithChaining.h
 *  PAL v1.0
 *
 *  Created by Karl Gemayel on 11/06/2012.
 *  Copyright 2012 American University of Beirut. All rights reserved.
 *
 */

/** 
 * \class BinaryHeapWithChaining
 * \brief A Max Heap implementation using the Binary Heap design with a tweak; elements
 * having equal degress are 'chained' together.
 *
 * The Binary Heap is a tree structure where each node has a maximum of two children
 * nodes. Moreover, in a max heap, the value stored in the parent node is greater than
 * or equal to the values stored in its children. This recursive definition results
 * in a partially ordered set of elements.
 */

#ifndef BinaryHeapWithChaining_h
#define BinaryHeapWithChaining_h

#include "../Heap/Heap.h"
#include "../Polynomial/Polynomial.h"
#include "../Options/Options.h"


using namespace std;


class BinaryHeapWithChaining : public Heap {
private:
	heap_size_t heapSize;		/**< the size of the heap */
	heap_size_t heapCapacity;	/**< the capacity of the heap */
	heap_t heap;			/**< the actual heap itself */

	/**
	 * Chaining Implementation:
	 * 	The chaining mechanism used is based on an indexing structure that
	 *	acts as a buffer between the heap elements and the elements in the
	 *	chains. It works in the following way:
	 *		After reading an monomial from the heap, its degree is found and
	 *		using this degree, we go to the appropriate indexing structure
	 *		that contains the locations of the chains. Every element in that
	 *		structure will contain a number that represents the actual 'physical'
	 *		distance between the END of the indexing structure, to the start of
	 *		the requested chain.
	 *		Each memory block representing a chain contains, as its first element,
	 *		a value that represents the chain's capacity and size (i.e. number of 
	 *		elements it contains). After that, space (of size equal to the chain's
	 *		capacity) is reserved to hold the chain's elements, where the first 'size'
	 *		elements are valid.
	 */
	heap_t::iterator chains;	/**< points to the start of the block of chains */
	heap_t::iterator chainDistances; /**< points to the element holding the parameter of the first chain */
	size_t numberOfChains;
	
	
	// heap statistics
	heap_size_t largestHeapSize;
	
// private functions
private:
	/** 
	 * A function that maintains the max-heap property. It is called on a certain
	 * which is suspected to violate that property
	 *
	 * \param i the index of the node we want to address
	 */
	void maxHeapify(size_t i);
	
	/** 
	 * Increases the value of the element in node (i) to the value keyDegree
	 *
	 * \param i the index of the node whose key we want to increase
	 * \param keyDegree the value of the new key
	 */
	void heapIncreaseKey(size_t i, monom_t keyDegree);
	
	
	// function prototypes related to chaining
	size_t calculateMemorySpaceForChaining();
	size_t getDistanceToChainOfDegree(deg_t degree);
	void addMonomialToChain (monom_t monomial);
	monom_t extractMonomialFromChainOfDegree(deg_t degree);
	void setupChainParameters ();
	bool getMonomial (deg_t degree, monom_t &monomial);
	bool degreeIsInHeap(deg_t degree);
	void setDegreePresent(deg_t degree, bool present);
	
public:
	// functions
	
	/**	
	 * A constructor that initializes the heap with the specified size.
	 *
	 * \param sizeOfHeap the size of the heap
	 */
	BinaryHeapWithChaining (heap_size_t sizeOfHeap);
	
	~BinaryHeapWithChaining() {
		// print some statistics of the heap
		if (verboseLevel(VERBOSE_LOW)) {
			cout << "Heap Statistics:" << endl;
			cout << "\tLargest Capacity allocated: " << heap.capacity() << endl;
			cout << "\tLargest Size reached      : " << largestHeapSize << endl;
			cout << "\tNumber of chains needed   : " << numberOfChains << endl;

			cout << endl;
		}
	}
	/**	
	 * Inserts an element that has a degree, coefficient, and an id into the heap
	 * in its correct position so as to not disturb the max-heap property.
	 *
	 * \param degree the degree of the element
	 * \param coef the coefficient of the element
	 * \param the id of the element
	 */
	void insert (deg_t degree, coef_t coef, ID_t f_id);
	void insert (deg_t degree, coef_t coef);
	
	/** 
	 * Returns the value of the max element of the heap without actually extracting it.
	 *
	 * @return The value of the max element 
	 */
	monom_t peek();
	
	/** 
	 * Extracts the max element of the heap and returns it.
	 *
	 * \return The value of the max element
	 */
	monom_t poll();
	
	/**
	 * Returns the size of the heap.
	 *
	 * \return The size of the heap
	 */
	heap_size_t size();
	
	/** 
	 * Returns true if the heap is empty, and false otherwise
	 *
	 * \return if heap is empty - true. Otherwise - false
	 */
	bool isEmpty();
};

#endif

