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

#include "../Options/Options.h"




// FUNCTION PROTOTYPES

inline size_t parent(size_t i);
inline size_t left(size_t i);
inline size_t right(size_t i);


size_t createDistancePresentPair(size_t distance, bool present) ;


/*
 initializeHeap: sets the heap capacity to n, and returns the array
 */
BinaryHeapWithChaining::BinaryHeapWithChaining (heap_size_t n) {
	heapSize = 0;
    heapCapacity = n;

	// calculate the space needed to hold the heap elements
	// as well as the overhead needed for chaining
	heap_size_t spaceToReserve = n; 			// space for heap elements
	spaceToReserve += calculateMemorySpaceForChaining();		// space for chaining elements and overhead
	
	if (verboseLevel(VERBOSE_LOW)) {
		cout << "Space to Reserve = " << spaceToReserve << endl;
		cout << "Heap Capacity = " << heapCapacity << endl;
	}


	// reserve space in vector
	heap.reserve(spaceToReserve);

	// setup default parameters for all chains
	setupChainParameters();
	
	// initialize variables that gather statistics
	largestHeapSize = 0;
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
//	cout << "Inserting - Size = " << heapSize << endl;
	if (heapSize == heapCapacity) {
		printWarning("Heap is full, can't insert element");
		return;
	}
	
	monom_t key;
	
	if (heapSize == 0) {		// heap is originally empty
		key = createMonomial(degree, coef, 0);
		heap[heapSize++] = key;
		
		updateStatistics();
		
		return;
	}
	
	
	
	// if a monomial of the same degree already exists in the heap, 
	// add this element to the chain of that degree
	if (degreeIsInHeap(degree)) {
		addMonomialToChain(createMonomial(degree, coef, f_id));
		if (verboseLevel(VERBOSE_HIGH)) {
			cout << "\tElement Chained" << endl;
		}
	}
	// otherwise, put the monomial at the end of the heap and perform heapIncreaseKey
	else {
		w_type keyDegree = degree;		// get current degree
		
		// create an 0 degree poly
		key = createMonomial(0, coef, f_id);
	
		heap[heapSize++] = key ;
		heapIncreaseKey(heapSize-1, keyDegree);
		
		// indicate that a monomial of that degree is now present in
		// the heap by setting the 'present' boolean to true
		setDegreePresent(keyDegree, true);
	}
	
	updateStatistics();
}

/*
 peek: peeks at the maximum value without extracting it
 */
monom_t BinaryHeapWithChaining::peek() {
	if (heapSize == 0) {
		// should throw exception
		return 0;
	}
	
	return heap[0];
}

/*
 poll: extracts the maximum value from the heap
 */
monom_t BinaryHeapWithChaining::poll() {
//	cout << "Polling - Size = " << heapSize << endl;
	if (heapSize < 0) {
		printf("ERROR: heap underflow");
		return NULL;
	}
	
	if (heapSize == 0) {
		return 0;
	}

	monom_t max = heap[0];
	heap[0] = heap[heapSize - 1];
	heapSize = heapSize-1;		// decrement heap size
	maxHeapify(0);

	// now get the elements from the chain relevant to the monomial extracted,
	// and add their coefficients to that monomial
	deg_t currDegree = GET_DEGREE(max);
	coef_t currCoef = GET_COEF(max);
	monom_t tempMonomial = 0;
	while (!chainIsEmpty(chains+getDistanceToChainOfDegree(currDegree))) {
		tempMonomial = extractMonomialFromChainOfDegree (currDegree);
		currCoef += GET_COEF(tempMonomial);
	}
	
	// set the (new) coefficient of the current monomial (max)
	setCoef(max, currCoef);
	
	// indicate that all monomial of that degree have been removed from
	// the heap by setting the 'present' boolean to false
	setDegreePresent(currDegree, false);
	
	return max;
}

/*
 maxHeapify: reorders the heap
 */

void BinaryHeapWithChaining::maxHeapify(size_t i) {
	size_t l = left(i);
	size_t r = right(i);
	
	size_t largest = -1;
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
	return (heapSize == 0 ? true : false);
}


heap_size_t BinaryHeapWithChaining::size() {
	return heapSize;
}



bool BinaryHeapWithChaining::getMonomial (deg_t degree, monom_t &monomial) {
	for (heap_size_t i = 0; i < heapSize; i++) {
		if (degree == GET_DEGREE(heap[i])) {
			monomial = heap[i];
			return true;
		}
	}
	
	return false;
}


// BELOW SECTION IS RELATED TO CHAINING

/* Packing Distance to chain and Element in Heap 
 * 	This packing packs two values together in each element of the indexing structure that
 *	holds the positions of all the chains. For the i'th element of the indexing structure
 *	the first value represents the distance to the chain of the i'th degree, and the second
 * 	value acts as a boolean that indicates whether an element of degree i is already present
 *	in the heap
 */
#define DISTANCE_MASK 0xFFFFFFFFFFFFFFFEULL	/**< mask representing distance */
#define PRESENT_MASK  0x0000000000000001ULL	/**< mask representing boolean indicating if element is in the heap	*/


#define DISTANCE_SHIFT 1          /**< displacement of distance from right */


// the below macros automatically unpack the distance and 'present' values
// based on the 'masks' and 'shifts' above
#define GET_DISTANCE(W) ((size_t) ((W & DISTANCE_MASK) >> DISTANCE_SHIFT))			/**< gets the ID of the monomial */
#define GET_PRESENT(W)  ((size_t) ((W & PRESENT_MASK)))		/**< gets the degreeof the monomial */


bool BinaryHeapWithChaining::degreeIsInHeap(deg_t degree) {
	return (GET_PRESENT(chainDistances[degree]) == PRESENT_MASK);
}

void BinaryHeapWithChaining::setDegreePresent(deg_t degree, bool present) {
	chainDistances[degree] = chainDistances[degree] & DISTANCE_MASK;		// set present to false
	
	if (present) 
		chainDistances[degree] = chainDistances[degree] | PRESENT_MASK;		// set present to true 
	
}



size_t createDistancePresentPair(size_t distance, bool present) {
	size_t elem = distance << DISTANCE_SHIFT;		// set distance value, present = false
	
	if (present) 
		elem = elem | PRESENT_MASK;
	
	return elem;
	
}



void BinaryHeapWithChaining::setupChainParameters () {
	chainDistances = heap.begin() + heapCapacity + 1;		// 'chainsDistances' now points to the start of the memory block containing the chains distances

	// get the filename of the file that contains the chaining info
	string infoFilename;
	getBHWithChainingFilename(infoFilename);
	
	// open the file containing chaining info
	ifstream cInfo (infoFilename.c_str());
	
	// The info file contains lines of numbers. The first line contains the number of chains
	// that we will need in the whole multiplication process. The remaining lines each contain 
	// two numbers seperated by the space character. The first of the two numbers represents 
	// the degree for which the chain will be created, and the second number represents the 
	// capacity of the chain of that degree.
	// NOTE: the lines containing degree-capacity pairs are MUST be sorted in increasing order of
	//		 degrees, and must contain EVERY degree, starting from 0 up to the largest degree used
	
	
	cInfo >> numberOfChains;		// set the number of chains
	chains = chainDistances + numberOfChains;		// 'chains' now points to the start of the memory block containing the chains

	
	//cParameter_t param = 0;
	size_t distanceIndex = 0;		// the index of the element (in the indexing structure) containing the distance to the chain
	deg_t currDegree = 0;
	const chain_size_t currChainSize = 0;
	chain_size_t currChainCapacity = 0;
	
	chain_t currChain = chains;		// point to the first chain
	
	// 'distanceToChain' represents the distance from the END of the indexing structure
	// to the start of the next chain
	size_t distanceToChain = 0;		
	
	int counter = 0;
	while (!cInfo.eof()) {
		counter++;
		// cout << "Counter = " << counter << endl;
		// read the degree and chain size
		cInfo >> currDegree;
		cInfo >> currChainCapacity;
		
		// set the distance/present values to the chain
		*(chainDistances + distanceIndex) = createDistancePresentPair(distanceToChain, false);
		
		// set the chain parameters having size 0
		setChainParameters(currChain, currChainSize, currChainCapacity);
		
//		cout << "cSize = " << getChainSize(currChain) << ", cCap = " << getChainCapacity(currChain) << endl;
		// update the "distance" to the next chain
		distanceToChain = distanceToChain + currChainCapacity + 1;		// the extra '1' accounts for the chain parameters

		distanceIndex++;		// increment the index for the next chain distance element
		if (!cInfo.eof()) currChain = currChain + currChainCapacity + 1;
	}
	
	
	// close the info file
	cInfo.close();
	
	// print for info
	if (verboseLevel(VERBOSE_HIGH)) {
		cout << "-----Chains Information-----" << endl;
		
		for (size_t i = 0; i < numberOfChains; i++) {
			cout << "Chain #" << i << ": Cap = " << getChainCapacity(chains + getDistanceToChainOfDegree(i)) << endl;
		}
	}
}

size_t BinaryHeapWithChaining::getDistanceToChainOfDegree(deg_t degree) {
	return GET_DISTANCE(*(chainDistances + degree)); 
}


void BinaryHeapWithChaining::addMonomialToChain (monom_t monomial) {
	// get the degree of the monomial
	deg_t degree = GET_DEGREE(monomial);
	
	
	// get a reference to the chain of that degree
	chain_t degChain = chains + getDistanceToChainOfDegree(degree);
	
	// make sure the chain (of predefined size) can hold more elements
	chain_size_t chainSize = getChainSize (degChain);
	
 	if (verboseLevel(VERBOSE_HIGH))
		cout << "Chain size = " << chainSize << "/" << getChainCapacity(degChain) << ". Adding degree = " << degree << ", Coef = " << GET_COEF(monomial) << endl;
	
	if (chainSize == getChainCapacity(degChain)) {
		printError("Chain is full. Possible Internal error");
		return;
	}
	
	// otherwise, add the new monomial to the chain and increment its size
	degChain[CHAIN_OVERHEAD + chainSize] = monomial;
	incrementChainSize(degChain);
}

monom_t BinaryHeapWithChaining::extractMonomialFromChainOfDegree(deg_t degree) {
	// get a reference to the chain of that degree
	chain_t degChain = chains + getDistanceToChainOfDegree(degree);
	
	chain_size_t chainSize = getChainSize(degChain);
	
	if (chainSize == 0) {
		// should throw exception
		return 0;
	}
	
	// otherwise, remove an element from the chain and decrement its size
	monom_t element = degChain[CHAIN_OVERHEAD + chainSize-1];			// get the last element in the chain
	decrementChainSize(degChain);
	

	
	return element;
}


/*
 * calculateMemorySpaceForChaining: computes the space needed to hold the chains,
 *	as well as any overhead needed by the chaining mechanism to perform the chaining
 *	process. This information is retrieved from the preprocessing information file.
 */
size_t BinaryHeapWithChaining::calculateMemorySpaceForChaining() {
	// open chaining information file
	
	string cInfoFilename;
	getBHWithChainingFilename(cInfoFilename);
	
	ifstream cInfo (cInfoFilename.c_str());
	
	size_t total = 0;
	
	size_t numOfChains;
	cInfo >> numOfChains;		// get the maxDegree
	
	total += numOfChains;		// we need 'maxDegree' elements to index into every chain
	
	deg_t currDegree;
	size_t currChainCapacity;
	// now add the capacity of each chain, along with the chain's overhead
	while (!cInfo.eof()) {
		
		// read the degree and chain size
		cInfo >> currDegree;
		cInfo >> currChainCapacity;
		
		total += currChainCapacity + CHAIN_OVERHEAD;		// add the chain capacity plus chain overhead
	}
	
	
	return total;
	
}




/**
 *	Updates the statistics of the heap
 */
inline void BinaryHeapWithChaining::updateStatistics() {
	largestHeapSize = (heapSize > largestHeapSize) ? heapSize : largestHeapSize;		// update largest heap size reached
}
