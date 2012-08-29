/*
 *  BinaryHeap.cpp
 *  Binary Heap With ID
 *
 *  Created by Karl Gemayel on 16/05/2012.
 *  Copyright 2012 American University of Beirut. All rights reserved.
 *
 */
#include <algorithm>
#include <vector>
#include <iostream>

#include "BinaryHeap.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

using namespace std;

// FUNCTION PROTOTYPES


// accessing children and parent
heap_size_t parent(heap_size_t i);
heap_size_t left(heap_size_t i);
heap_size_t right(heap_size_t i);

/*
 * Function Implementation
 */

/*
 * initializeHeap: sets the heap capacity to n, and returns the array
 *
 * @param n: the heap's capacity
 */
BinaryHeap::BinaryHeap(heap_size_t n) {
	heapSize = 0;
    heapCapacity = n;
	
	// reserve space for the heap
	heap.reserve(n);
	
	// initialize variables that gather statistics
	largestHeapSize = 0;
}


/*
 * heapIncreaseKey: increases the degree of the i'th element to key
 *
 * @param i: the index of the mononmial
 * @param keyDegree: the degree to which the mononmial's degree will be set
 */

void BinaryHeap::heapIncreaseKey(heap_size_t i, deg_t keyDegree) {
	
	if (keyDegree < GET_DEGREE(heap[i])) {			// if keyDegree < degree of i'th monomial
		printf("ERROR: new key is smaller than current key\n");
		return;
	}
	

    setDegree(heap[i], keyDegree);    // set the new degree
	
	// loop while i > 0 and Degree(parent(i)) < Degree(i)
	while (i > 0 && (GET_DEGREE(heap[parent(i)])) < (GET_DEGREE(heap[i]))) {
        
        // swap the two elements in the heap
		monom_t temp = heap[i];
		heap[i] = heap[parent(i)];
		heap[parent(i)] = temp;
		
		i = parent(i);      // set i to current position
	}	
}

/*
 * insert: inserts a word into its correct position
 *
 *  @param heap:    the max heap
 *  @param degree:  the degree to be inserted, where degree >= 0
 *  @param coef:    the coefficient to be inserted, where coef >= 0
 *  @param f_id:    the ID of the current f monomial, where f_id >= 0
 */

void BinaryHeap::insert(deg_t degree, coef_t coef, ID_t f_id) {
    
    if (heapSize == heapCapacity) {
		printWarning("Heap is full, can't insert element");
        return;
    }
	
	monom_t key;
	
	if (heapSize == 0) {		// heap is originally empty
		key = createMonomial(degree, coef, f_id);         // pack the elements into a word
		heap[heapSize++] = key;                            // insert new element into the heap
		
		updateStatistics();			// update heap statistics used for performance evaulation
		
		return;
	}
	
	
	deg_t keyDegree = degree;		// get current degree
	
	// create a 0 degree poly
	key = createMonomial(0, coef, f_id);
	
	heap[heapSize++] = key;
	heapIncreaseKey(heapSize-1, keyDegree);
	
	updateStatistics();			// update heap statistics used for performance evaulation
	
}


void BinaryHeap::insert (deg_t degree, coef_t coef) {
	insert (degree, coef, 0);			// insert an element into the heap with ID=0
}

/*
 peek: peeks at the maximum value without extracting it
 */
monom_t BinaryHeap::peek() {
	if (heapSize == 0) {
		// should throw exception
		return 0;
	}
	
	return heap[0];
}

/*
 poll: extracts the maximum value from the heap
 */
monom_t BinaryHeap::poll() {
	if (heapSize < 0) {
		printf("ERROR: heap underflow");
		return 0;
	}
	
	if (heapSize == 0) {			// heap is empty
		return 0;
	}
	
	monom_t max = heap[0];
	heap[0] = heap[heapSize - 1];
	heapSize = heapSize-1;		// decrement heap size
	maxHeapify(0);
	
	return max;
}

/*
 maxHeapify: reorders the heap
 */

void BinaryHeap::maxHeapify(heap_size_t i) {
	heap_size_t l = left(i);
	heap_size_t r = right(i);

    heap_size_t largest = -1;
    if (l <= heapSize && GET_DEGREE(heap[l]) > GET_DEGREE(heap[i])) {
        largest = l;
    }
    else {
        largest = i;
    }
    
    if (r < heapSize && GET_DEGREE(heap[r]) > GET_DEGREE(heap[largest])) {
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
inline heap_size_t parent(heap_size_t i) {
    if (i == 0)
        return 0;
    
    if (i % 2 == 0) {
        return (heap_size_t) (i/2) - 1;
    }
    
	return (heap_size_t) (i/2);
}

/*
 left: returns the index of the left child of i
 */
inline heap_size_t left(heap_size_t i) {
    
	return 2*i + 1;
}

/*
 right: returns the index of the right child of i
 */
inline heap_size_t right(heap_size_t i) {
	return 2*i + 2;
}


bool BinaryHeap::isEmpty() {
	return (heapSize == 0 ? true : false);
}


heap_size_t BinaryHeap::size() {
	return heapSize;
}



/**
 *	Updates the statistics of the heap
 */
inline void BinaryHeap::updateStatistics() {
	largestHeapSize = (heapSize > largestHeapSize) ? heapSize : largestHeapSize;		// update largest heap size reached
}