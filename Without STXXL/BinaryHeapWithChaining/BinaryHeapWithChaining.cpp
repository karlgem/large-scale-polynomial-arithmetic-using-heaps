/*
 *  Polynomial_Heap.cpp
 *  Multiple Pairs With Chaining (2)
 *
 *  Created by Karl Gemayel on 7/21/11.
 *  Copyright 2011 American University of Beirut. All rights reserved.
 *
 */

#include "BinaryHeapWithChaining.h"
#include "Chain.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <fstream>
#include <string>


// FUNCTION PROTOTYPES

inline size_t parent(size_t i);
inline size_t left(size_t i);
inline size_t right(size_t i);


/*
 initializeHeap: sets the heap capacity to n, and returns the array
 */
BinaryHeapWithChaining::BinaryHeapWithChaining (size_t n, size_t numOfChains) {
//	heapSize = 0;
//  heapCapacity = n;
    numberOfChains = numOfChains;
	largestHeapSize = 0;
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Heap Capacity = " << n << endl;

	// reserve space in vector
	heap.reserve(n);
	
	// allocate space for chains
	chains = new Chain *[numberOfChains];
	
	// intiially, set all chains to null...chains are allocated as necessary
	for (size_t i = 0; i < numberOfChains; i++) {
		chains[i] = NULL;
	}
}


/*
 * heapIncreaseKey: increases the degree of the i'th element to key
 *
 * @param A: the current monomial A
 * @param i: index of element to be increased
 * @param key: value of the current degree
 */
void BinaryHeapWithChaining::heapIncreaseKey(size_t i, deg_t keyDegree) {
	
	if (keyDegree < GET_DEGREE(heap[i])) {			// if key < degree of i'th monomial
		printf("ERROR: new key is smaller than current key\n");
		return;
	}
	
	setDegree(heap[i], keyDegree);    // set the new degree
	
	// loop while i > 0 and Degree(parent(i)) < Degree(i)
	while (i > 0 && (GET_DEGREE(heap[parent(i)]) < (GET_DEGREE(heap[i])))) {		
		// swap the two elements in the heap
		monom_t temp = heap[i];
		heap[i] = heap[parent(i)];
		heap[parent(i)] = temp;

		i = parent(i);      // set i to current position
	}	
}


void BinaryHeapWithChaining::insert(deg_t degree, coef_t coef) {
	insert(degree, coef, 0);
}
/*
 insert: inserts a word into its correct position
 */

void BinaryHeapWithChaining::insert(deg_t degree, coef_t coef, ID_t f_id) {
	
	size_t cap = heap.capacity();
	bool full = heap.size() == cap;
	
	
	Chain *currChain = chains[degree];
	
	// allocate space for a chain of that degree if it doesn't exit
	if (chains[degree] == NULL) {
		if (verboseLevel(VERBOSE_HIGH)) cout << "Allocating chain for degree " << degree << endl;
		chains[degree] = new Chain();
		currChain = chains[degree];
	}
	
	monom_t key;
	
	if (heap.size() == 0) {		// heap is originally empty
		key = createMonomial(degree, coef, 0);
		//HERE heap[heapSize++] = key;
		heap.push_back(key);
		return;
	}

	
	// if a monomial of the same degree already exists in the heap, 
	// add this element to the chain of that degree
	if (currChain != NULL && currChain->isValid()) {		
		currChain->push(createMonomial(degree, coef, f_id));
		if (verboseLevel(VERBOSE_HIGH)) {
			cout << "\tElement Chained" << endl;
		}
	}
	// otherwise, put the monomial at the end of the heap and perform heapIncreaseKey
	// and set the chain of that degree to be valid
	else {
		w_type keyDegree = degree;		// get current degree
		
		// create an 0 degree poly
		key = createMonomial(0, coef, f_id);
	
		//HERE heap[heapSize++] = key ;
		heap.push_back(key);
		heapIncreaseKey(heap.size()-1, keyDegree);
		
		// indicate that a monomial of that degree is now present in
		// the heap by setting the 'present' boolean to true
		currChain->setValid(true);
		
		// update the variable containing the largest heap size reached
		if (largestHeapSize < heap.size()) largestHeapSize = heap.size();
		
		// if heap was full before insert (i.e. if we increased the capacity
		if (full && verboseLevel(VERBOSE_LOW)) {
			cout << "Increase heap capacity from " << cap << " to " << heap.capacity() << endl;
		}
	}
}

/*
 peek: peeks at the maximum value without extracting it
 */
monom_t BinaryHeapWithChaining::peek() {
	if (heap.empty()) {
		// should throw exception
		return 0;
	}
	
	return heap[0];
}

/*
 poll: extracts the maximum value from the heap
 */
monom_t BinaryHeapWithChaining::poll() {
	
	if (heap.empty()) {
		// should throw exception
		return 0;
	}

	monom_t max = heap[0];
	heap[0] = heap[heap.size() - 1];
//	heapSize = heapSize-1;		// decrement heap size
	heap.pop_back();
	maxHeapify(0);

	// now get the elements from the chain relevant to the monomial extracted,
	// and add their coefficients to that monomial
	deg_t currDegree = GET_DEGREE(max);
	coef_t currCoef = GET_COEF(max);
	
	Chain *currChain = chains[currDegree];
	monom_t tempMonomial = 0;
	while (!currChain->isEmpty()) {
		tempMonomial = currChain->pop();
		currCoef += GET_COEF(tempMonomial);
	}
	
	// set the (new) coefficient of the current monomial (max)
	setCoef(max, currCoef);
	
	// indicate that all monomial of that degree have been removed from
	// the heap by setting the 'present' boolean to false
	currChain->setValid(false);
	
	
	// Now that the all monomials of the current degree have been removed from the heap and chains,
	// and given that we know that no new monomials of that degree will be inserted in the heap again,
	// we can safely de-allocate space for the chain the holds monomials of that degree
	delete currChain;
	chains[currDegree] = NULL;
	
	return max;
}

/*
 maxHeapify: reorders the heap
 */

void BinaryHeapWithChaining::maxHeapify(size_t i) {
	size_t l = left(i);
	size_t r = right(i);
	
	size_t largest = -1;
    if (l <= heap.size() && GET_DEGREE(heap[l]) > GET_DEGREE(heap[i])) {
        largest = l;
    }
    else {
        largest = i;
    }
    
    if (r < heap.size() && GET_DEGREE(heap[r]) > GET_DEGREE(heap[largest])) {
        largest = r;
    }
    
    if (largest != i) {
        // exchange heap[i] with heap[largest]
        monom_t temp = heap[i];
        heap[i] = heap[largest];
        heap[largest] = temp;
        
        maxHeapify(largest);
    }	
}




/*
 parent: returns the index of the parent of i
 */
inline size_t parent(size_t i) {	
	if (i == 0)
        return 0;
    
    if (i % 2 == 0) {
        return (size_t) (i/2) - 1;
    }
    
	return (size_t) (i/2);
}

/*
 left: returns the index of the left child of i
 */
inline size_t left(size_t i) {
    
	return 2*i + 1;
}

/*
 right: returns the index of the right child of i
 */
inline size_t right(size_t i) {
	return 2*i + 2;
}



bool BinaryHeapWithChaining::isEmpty() {
	return (heap.size() == 0 ? true : false);
}


size_t BinaryHeapWithChaining::size() {
//HERE	return heapSize;
	return heap.size();
}



