//
//  FunnelHeap.cpp
//
//  Created by Karl Gemayel on 23/08/12.
//  Copyright 2012 American University of Beirut. All rights reserved.


#include <iostream>
#include <cstdlib>
#include <math.h>
#include <stack>
#include <queue>
#include "FunnelHeap.h"

#include <assert.h>

using namespace std;


// Enable assertions to ensure that the functions are operating correctly on a technical level
#define ASSERTIONS true


#define LINK_OVERHEAD (int) 0       // number of elements in beginning of a link or
#define BUFFER_OVERHEAD (int) 0     // overhead amount associated with each buffer




#pragma mark -
#pragma mark Function Prototypes

/* Function Prorotypes */
inline void appendBufferToStream(Buffer* buff, stream_t &stream);
inline monom_t pop (stream_t& stream1, stream_t& stream2);



// computations of sizes and capacities
inline void getKAndS (long values [], int i);
inline heap_size_t computeLinkRequiredStorage(int i);



#pragma mark -
#pragma mark Header Functions

/*
 * initializeQueue: Allocates memory for the queue, and returns a pointer to 
 *                  its head element
 */
FunnelHeap::FunnelHeap (int numberOfLinks) {  // "numOfLinks" should be derived depending on input size (DO LATER)
	numOfLinks = numberOfLinks;
	
	_size = 0;
	// create link pointers
	links = new Link* [numberOfLinks];
	
    const int s1 = 8;
    
    
	// compute the required storage amount for the entire heap
    int I_size = s1 + BUFFER_OVERHEAD;       // size of I is s1 elements and the buffer overhead
    
    heap_size_t storage = I_size + 2;		// size of buffer I + 2 (for size of queue and I_NULL element)
    
    for (int i = 1; i <= numOfLinks; i++) {
        storage += computeLinkRequiredStorage(i) + LINK_OVERHEAD;
    }
    
    
	if (verboseLevel(VERBOSE_LOW))
		printf("Queue size (in MB) = %llu\n", (monom_t) ((storage / ( 1024.0 * 1024.0) * sizeof(monom_t) )));
														
	// reserve space for the heap
	heap.reserve(storage);
		
	// setup the heap
    setupHeap();  
}


bool FunnelHeap::isEmpty() {
//	return _size == 0;
	return I->isEmpty() && getLink(1)->A->isExhausted();
}


monom_t FunnelHeap::peek() {
    
    Buffer* A_1 = I->getLeftChild();
    
    if (A_1->isEmpty() && !A_1->isExhausted()) {
        A_1->fill();
    }
    
	bool I_empty = I->isEmpty();
	bool A_1_empty = A_1->isEmpty();
	
	// if both buffers are empty, return 0 (ie nothing)
    if (I_empty && A_1_empty) {
        return 0;
    }
    
	// if neither is empty, then take the max element from both
    if (!I_empty && !A_1_empty) {
		if (GET_DEGREE(I->front()) > GET_DEGREE(A_1->front())) {		// if I contains the max element
        	return I->front();
    	}
		else 		// if A_1 contains the max element
		    return A_1->front();
	}
	// if I is empty but A_1 isn't, take max from A_1
	else if (I_empty && !A_1_empty) {
		return A_1->front();
	}
	// if A_1 is empty but I isn't, take max from I
	else
		return I->front();
  
}


void FunnelHeap::insert (deg_t degree, coef_t coef, ID_t id) {
	insert (degree, coef);
}

void FunnelHeap::insert (deg_t degree, coef_t coef) {
    
    /*
     * General Code Sketch:
     *
     *  1)  Start by inserting the new element in the buffer I.
     *  2)  If I becomes full, then sweep to free up space for
     *      new elements to be inserted in the heap
     */
	
	
    // insert new element in I
	I->insertSorted(createMonomial(degree, coef, 0));
	
	// increment the size of the heap
	_size++;
    
    // if I is full, then return. Otherwise, move on to the sweeping
    if (!I->isFull()) {
        return;
    }

	// cout << "I-Size = " << I->getSize() << ", I-Cap = " << I->getCapacity() << endl;
    
    // locate empty leaf to sweep to
    
    // get first empty leaf
    long c_i = 1;   
    
    
    Link* link_i = getLink(1);
    
    
    // We sweep in order to free up space in the top of the heap so
    // that we could insert new elements. This is done by shifting
    // elements further down into the heap while maintaining their
    // sorted order, thus freeing up space at the top.
    
    // The buffer to which we sweep is selected based on one of two conditions:
    //      1)  if c_i <= the number of leafs in a link i
    //      2)  if the number of elements in the lower part of a link is
    //          less than or equal to (k_i * s_i)/2
    int i = 1;
    while (true) {
        c_i = link_i->c;
        
//		long k_i = constantK[i-1];
		long k_i = getK(i);
        
        if (c_i <= k_i) {      // then we found the first empty buffer
            break;
        }
    
        // go to next link (if it exists)
		i++;
        if (i > numOfLinks) {           // if there are no more link, don't insert
            printWarning("Queue is full. Element cannot be inserted!");
            return;
        }
		else
        	link_i = getLink(i);
    }

    sweep(i);
}


monom_t FunnelHeap::poll () {
	assert(_size > 0);			// sanity check
	
	monom_t elem = pollInternal();  
	
	_size--;		// decrement size
	
	return elem;
}

heap_size_t FunnelHeap::size () {
	return _size;
}


#pragma mark -
#pragma mark Initilization of Queue (sub methods)



/*
 * setQueueDefaultValues: sets all the default values of links and buffers in their correct
 *                          location in the heap structure
 */
void FunnelHeap::setupHeap () {
    const int s1 = 8;
    
	mem_t it = heap.begin();
	// let 'mem' point to the start of the queue (i.e. to the queue size element)
	//heap_t::iterator mem = heap.begin();
	mem_t mem = it;
	(*mem) = 0;			// set the initial queue size to 0
	mem++;			// let mem point to the start of buffer I
		
	// create buffer I
	I = new Buffer (mem, s1, I_BUF_T);

	

	// setup the links
	int linkIndex;
    for (int i = 0; i < numOfLinks; i++) {
		linkIndex = i+1;
		
		// if (linkIndex == 6)
		// 	setGlobalBoolean(true);
			
		links[i] = new Link (mem, linkIndex);
		
		// for the first link, set the parent of A_1 to I, and the children of I to A_1
		if (i == 0) {
			links[0]->A->setParent(I);
			I->setLeftChild(links[0]->A);
			I->setRightChild(links[0]->A);
		}
		// for the remaining links (i >= 1), set the parent of A_i to A_i-1, and the right child of A_i-1 to A_i
		else {
			links[i]->A->setParent(links[i-1]->A);
			links[i-1]->A->setRightChild(links[i]->A);
		}
    }
}



#pragma mark -
#pragma mark Traversing Queue

/**
 *	Returns a pointer to the specified link. Note: the first link is indexed by 1
 *
 *	@param linkIndex the index of the link (1 <= linkIndex <= numberOfLinks)
 */
Link* FunnelHeap::getLink(int linkIndex) {
	assert(linkIndex >= 1 && linkIndex <= numOfLinks);			// sanity check
	
	return links[linkIndex-1];
}




#pragma mark -
#pragma mark Computing Capacities

/* 
 * getSAndK: computes s_i and k_i of link i, and puts the value of k_i at 
 *           index 0, and s_i at index 1 in the values array
 */
inline void getKAndS (long values [], int i) {
    
    /************************ Definition of k_i and s_i *************************
     *                                                                          *
     * (k_1, s_1) = (2,8)                                                       *
     * s_{i+1} = s_i (k_i + 1)                                                  *
     * k_{i+1} = 2^{log((s_{i+1})^{1/3})}   (where log(x^{1/3}) is rounded up)  *
     *                                                                          *
     ****************************************************************************/
    
    long k_i = 2, s_i = 8;
    
    for (int count = 2; count <= i; count++) {
        s_i = s_i * (k_i + 1);
        
        long temp = ceil(log2(pow(s_i, 1/3.0)));
        k_i = (long) pow(2.0, (double) temp);
    }
    
    values[0] = k_i;
    values[1] = s_i;
}

/*
 * computeLinkSize: computes the size of the i'th link (the number of monom_t elements)
 *                  NOTE: The link capacity is the number of elements that can fit 
 *                  from A_i to s_iki
 */
unsigned long long computeLinkRequiredStorage(int i) {
    
    /********************* Definition of Link i *************************
     *                                                                  *
     * Link capacity and index                                          *
     * c_i                                                              *
     * k_i buffers, each of size s_i                                    *
     * K_i merger (funnel) that mergers the buffers                     *
     * an output buffer B_i of size (k_i)^3                             *
     * v_i merger, that mergers b_i with A_{i+1}                        *
     * A_i buffer, the output of the v_i merger, of size (k_i)^3        *
     *                                                                  *
     * Structure of link i in memory:                                   *
     *      _________________________________________________________   *
     *     | Cap | Index | c_i | A_i | B_i |  K_i  | S_i1 |...| S_ik |  *
     *      ---------------------------------------------------------   *
     *                                                                  *
     ********************************************************************/
    
   
    long kAndS [2];
    getKAndS(kAndS, i);
    
    unsigned long long k_i = kAndS[0];
    unsigned long long s_i = kAndS[1];
    
 
    unsigned long long totalSize = k_i * s_i + (k_i * BUFFER_OVERHEAD);     // num of elem in each buffer plus buffer overhead for each buffer
    totalSize += 2 * ((k_i*k_i*k_i) + BUFFER_OVERHEAD) ;   // num of elements in A_i and B_i plus the packing and index
    
    totalSize += kSize(k_i) + (k_i - 2) * BUFFER_OVERHEAD;      // compute size of k-merger
    
    
    return totalSize;      
}



#pragma mark -
#pragma mark Auxiliary Funnel Heap Functions



/*
 *  sweep: used to empty buffer I and make room for new elements to be inserted
 *
 *  @param q: pointer to the q 
 *  @param linkIndex: the index of the link where the values in I will be "sweeped" into
 */
void FunnelHeap::sweep (int linkIndex) {
    /* Sweep Algorithm *\
     *
     * - Effect: + Move the content of links 1,...,i-1 to the 
     *             buffer S_ici, where i = linkIndex.
     *           + This means that all links less than i will 
     *             be empty.
     *
     * - Code Sketch:
     *      1) Traverse path 'p' from A_1 to S_ici, and:
     *          a) Record how many elements each buffer on this
     *             path currently contains.
     *          b) From sorted stream_1 of the elements in the 
     *             part of path 'p' from A_i to S_ici
     *      2) From sorted stream_2 of all elements in links less
     *         than i and in buffer I, by marking A_i as exhausted
     *         and callining poll() repeatedly.
     *      3) Merge the two streams, and traverse 'p' again while
     *         inserting the front elements of the new stream in
     *         the buffers on 'p' in such a way that these buffers
     *         contain the same numbers of elements as before the
     *         insertion, and then insert the remaining elements
     *         in S_ici.
     *      4) Reset the 'c' counter for all the links less than i
     *         back to 1, and increment c_i by 1.
     *  
     *
     */

   	// store all sizes of the A-buffers of all links up to (BUT EXCLUDING) A_sweep
	buf_size_t ASizes [linkIndex-1];
	Link* currLink;
	for (int i = 0; i < linkIndex-1; i++) {
		currLink = getLink(i+1);
		ASizes[i] = currLink->A->getSize();
	}
	
	// if (verboseLevel(VERBOSE_HIGH)) cout << "Forming stream 1" << endl;
    // store all sizes of path from A_sweep to S_ci (in which case S_ci is the S buffer being sweeped to)
	stream_t stream1;
 	Link* link_sweep = getLink(linkIndex);


	vector<buf_size_t> stream1BufferSizes;		// includes A_sweep
	
	// start with A_sweep
	Buffer* currBuff = link_sweep->A;		// point to A_sweep
	stream1BufferSizes.push_back(currBuff->getSize());

	// append elements to stream1
	appendBufferToStream(currBuff, stream1);

	// move on to B_sweep
	currBuff = currBuff->getLeftChild();		// point to B_sweep
	
	// traverse the path from B_sweep to S_ci, appending the buffers to stream1 and recording their sizes.
	int currC = link_sweep->c;
//	long currK = constantK[linkIndex-1];
	long currK = getK(linkIndex);
	
	// loop until we reach S_ci
    while (currBuff->getType() != LEAF_BUF_T) {
		stream1BufferSizes.push_back(currBuff->getSize());
		appendBufferToStream(currBuff, stream1);
		step(currBuff, currC, currK);
	}
	
	// set A_sweep as exhausted (possibly break invariant temporarily, intentionally);
	link_sweep->A->setExhausted(true);

    // STATE: stream1 contains all elements in path  A_sweep to S_ci
    
    // Forming stream2: Perform an internal extraction (i.e. does not affect the size parameter of the heap).
	// Having already marked A_sweep as exhausted, this extraction will move all elements in the links
	// before link_sweep to stream2.
	stream_t stream2;
	
	// while I is not empty, or A_1 is not exhausted
	while (!(I->isEmpty() && getLink(1)->A->isExhausted())) {
		stream2.push(pollInternal());
	}
	
	
	// Re-inserting elements: Conceptually, the two streams should be merged, and then reinserted into the heap.
	// But for performance, the streams will not be merged together into a thrid stream, rather they will act 
	// together as one virtual stream, meanining that when we need an element from this virtual (third) stream, 
	// we take the maximum element of the two streams.
	// The insertion involves inserting elements into the A-buffers of the links before link_sweep, preserving
	// the sizes they had before this sweep operation. Then, the remaining elements are stored in the path
	// from A_sweep to S_ci (also preserving previous sizes), and the elements that remain are stored in S_ci
	
	// start by filling the A-buffers of the links before link_sweep
	for (int i = 1; i < linkIndex; i++) {
		currLink = getLink(i);
		refillBuffer(currLink->A, stream1, stream2, ASizes[i-1]);
		currLink->A->setExhausted(false);			// set not exhausted
	}
	
	// fill A_sweep
	refillBuffer(link_sweep->A, stream1, stream2, stream1BufferSizes[0]);
	
	link_sweep->A->setExhausted(false);		// reset A_sweep to not be exhausted
	
	// fill the path from B_sweep to S_ci
	currBuff = link_sweep->A->getLeftChild();		// points to B_sweep
	currC = link_sweep->c;
//	currK = constantK[linkIndex-1];
	currK = getK(linkIndex);
	
	for (size_t i = 1; i < stream1BufferSizes.size(); i++) {
		refillBuffer(currBuff, stream1, stream2, stream1BufferSizes[i]);
		currBuff->setExhausted(false);			// set not exhausted
		step(currBuff, currC, currK);
	}
	

	// fill the leaf with the remaining elements
	refillBuffer(currBuff, stream1, stream2, stream1.size() + stream2.size());
	currBuff->setExhausted(false);			// set not exhausted
        
    // reset c_i and r_i (for i < linkindex) to 1
    for (int i = 1; i <= linkIndex-1; i++) {
		currLink = getLink(i);
		currLink->c = 1;
    }
    
    
    // increment c_i (for i = linkIndex) by 1
	link_sweep->c = link_sweep->c + 1;
}

monom_t FunnelHeap::pollInternal() {
	Buffer* A_1 = I->getLeftChild();

    if (A_1->isEmpty() && !A_1->isExhausted()){
        A_1->fill();
	}

	bool I_empty = I->isEmpty();
	bool A_1_empty = A_1->isEmpty();;

	// if both buffers are empty, return 0 (ie nothing)
    if (I_empty && A_1_empty) {
        return 0;
    }

	monom_t popped = 0;
	// if neither is empty, then take the max element from both
    if (!I_empty && !A_1_empty) {
		if (GET_DEGREE(I->front()) >= GET_DEGREE(A_1->front())) {		// if I contains the max element
			popped = I->pop_front();
    	}
		else {		// if A_1 contains the max element
			popped = A_1->pop_front();
		}
	}
	// if I is empty but A_1 isn't, take max from A_1
	else if (I_empty && !A_1_empty) {
		popped = A_1->pop_front();
	}
	// if A_1 is empty but I isn't, take max from I
	else {
		popped = I->pop_front();
	}

	return popped;
}


#pragma mark -
#pragma mark Stream Operations

/* 	peek: returns the maximum element from the two streams, without
 *			actually extracting it or performing merging of any kind
 */
inline monom_t peek (stream_t& stream1, stream_t& stream2) {
	
	bool empty1 = stream1.empty();
	bool empty2 = stream2.empty();
	
	if (empty1 && empty2)
		return 0;
		
	// if both streams are empty, compare their maxes
	if (!empty1 && !empty2) {
		return (GET_DEGREE(stream1.front()) > GET_DEGREE(stream2.front()) ? stream1.front() : stream2.front());
	}
	// if stream1 empty but stream2 isn't
	else if (empty1 && !empty2) {
		return stream2.front();
	}
	// if stream2 empty but stream1 isn't
	else
		return stream1.front();
}


/*
 *  pop: extracts the maximum element from the two streams.
 *
 *  @param stream1: the first stream of elements
 *  @param stream2: the second stream of elements
 */
inline monom_t pop (stream_t& stream1, stream_t& stream2) {
    /*
     * General Code Sketch:
     *
     *  1) Extract the element that has the highest degree in both streams
     *  2) Return it
     */
    
    if (stream1.empty() && stream2.empty())
        return 0;
    
	
    // Find the max element in both streams
    monom_t maxElement = 0;
    if (!stream1.empty())
        maxElement = stream1.front();
    
    if (!stream2.empty()) {
        if (GET_DEGREE(stream2.front()) >= GET_DEGREE(maxElement)) {
            // the max is from stream2
            maxElement = stream2.front();
            stream2.pop();
        }
        else {
            // the max is from stream1
            stream1.pop();
        }
    }
    else {      // the max was from stream1 (since stream2 is empty)
        stream1.pop();
    }
	
    return maxElement;      
}

inline bool empty(stream_t& stream1, stream_t& stream2) {
    return stream1.empty() && stream2.empty();
}


void appendBufferToStream(Buffer* buff, stream_t &stream) {
	// append the buffer to the end of the stream
	while (!buff->isEmpty()) {
		stream.push(buff->pop_front());
	}
}

void FunnelHeap::refillBuffer(Buffer* buff, stream_t &stream1, stream_t &stream2, buf_size_t numberOfElements) {
	// make sure the buffer can fit <numberOfElements> and the stream has enough elements
	// cout << "NOE = " << numberOfElements << ", B_SIZE = " << buff->getSize() << ", B_CAP = " << buff->getCapacity() << endl;
	assert(numberOfElements <= buff->getCapacity());
	assert(numberOfElements <= stream1.size() + stream2.size());
	
	for (buf_size_t i = 0; i < numberOfElements; i++) {
		buff->push_back(pop(stream1, stream2));		
	}
}






// PRINTING THE HEAP
void FunnelHeap::print() {
	// print the I-buffer
	cout << "I: "; I->print();
	
	// print all links
	for (long i = 1; i <= numOfLinks; i++) {
		Link* currLink = getLink(i);
		
		// if the link is exhausted, stop here
		if (currLink->A->isExhausted())
			break;
		
		// otherwise, print the link
		cout << "Link " << i << ": " << endl;
		currLink->print(); cout << endl;
		
	}
	
	cout << endl;
	
}