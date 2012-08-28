/*
 *  BinaryHeap.h
 *  PAL v1.0
 *
 *  Created by Karl Gemayel on 16/05/2012.
 *  Copyright 2012 American University of Beirut. All rights reserved.
 *
 */

/** 
 * \class BinaryHeap
 * \brief A Max Heap implementation using the standard Binary Heap design.
 *
 * The Binary Heap is a tree structure where each node has a maximum of two children
 * nodes. Moreover, in a max heap, the value stored in the parent node is greater than
 * or equal to the values stored in its children. This recursive definition results
 * in a partially ordered set of elements.
 */

#ifndef BinaryHeap_h
#define BinaryHeap_h

#include "../Heap/Heap.h"
#include "../Polynomial/Polynomial.h"
#include "../Options/Options.h"

using namespace std;


class BinaryHeap : public Heap {
private:
	heap_size_t heapSize;		/**< the size of the heap */
	heap_size_t heapCapacity;	/**< the capacity of the heap */
	heap_t heap;			/**< the actual heap itself */
	
	// heap statistics
	heap_size_t largestHeapSize;		// the largest size reached by the heap
	
	
// private functions
private:
	/** 
	 * A function that maintains the max-heap property. It is called on a certain
	 * which is suspected to violate that property
	 *
	 * \param i the index of the node we want to address
	 */
	void maxHeapify(heap_size_t i);
	
	/** 
	 * Increases the value of the element in node (i) to the value keyDegree
	 *
	 * \param i the index of the node whose key we want to increase
	 * \param keyDegree the value of the new key
	 */
	void heapIncreaseKey(heap_size_t i, monom_t keyDegree);
	
public:
	// functions
	
	/**	
	 * A constructor that initializes the heap with the specified size.
	 *
	 * \param sizeOfHeap the size of the heap
	 */
	BinaryHeap (heap_size_t sizeOfHeap);
	
	~BinaryHeap() {
		if (verboseLevel(VERBOSE_LOW)) {
			cout << "Heap Statistics:" << endl;
			cout << "\tLargest Capacity allocated: " << heap.capacity() << endl;
			cout << "\tLargest Size reached      : " << largestHeapSize << endl;
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

