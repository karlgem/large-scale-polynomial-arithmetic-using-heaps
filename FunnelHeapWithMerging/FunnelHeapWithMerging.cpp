//
//  CO_PriorityQueue.cpp
//  CO Without IDs
//
//  Created by Karl Gemayel on 05/03/12.
//  Modifed by Karl Gemayel on 05/03/12.
//  Copyright 2012 American University of Beirut. All rights reserved.

/*
 *	Specifications:
 *		1) Funnel Heap-based implementation of a (max) priority queue.
 *		2) Merges elements with equal degrees together
 *
 */



/********************************** MEMORY LAYOUT **********************************\
 *          ________________________________________________________
 * QUEUE:  | Queue Size | I | Link 1 | Link 2 | ... | Link i | ... |     
 *          -------------------------------------------------------
 *          ___________________________________________________________________
 * Link i: | Capacity | L_i | c_i | r_i | A_i | B_i | K_i | S_i1 | ... | S_iki |
 *          -------------------------------------------------------------------
 *         |_____________ ______________|
 *                       |
 *                 Link Overhead
 *
 *                                  ____________________________________________________
 * I, A_i, B_i, S_i1...S_iki:      | P_i | H_i | Left | Right | Parent | HT |  elements |
 *                                  ----------------------------------------------------
 *                                 |___________________ ____________________|
 *                                                     |
 *                                              Buffer Overhead
 *
 * NOTES:
 *  1) The link Capacity is the number of elements that can fit from c_i to S_iki
 *  2) P_i contains the size, capacity, isLeaf, and isExhausted bool of a buffer (details below)
 *  3) H_i contains the height of the tree, and the index of the buffer (details below)
 *  4) L_i contains the current link index, and the total number of links (details below)
 *  5) HT contains the the indices of the head and tail of the elements
 *  6) r_i represents the number of elements in the lower part of the link (ie B_i till leafs)
 *  7) Buffer Index is the index of the buffer:
 *          i) Index of I is the max number that can fit
 *         ii) Index of A_i is 0
 *        iii) Indices of B_i to S_ik (including K_i) are defined as follows:
 *                  - B_i is 1
 *                  - S_ik is the number of buffers from B_i to S_ik (including K_i)
 *                  - The rest are in sequential order, from left to right, by level.
 *
 * General Info:
 *  1) Links are indexed starting from 1
 *  2) Leafs are indexed starting from 1
 *  3) The root of every K-merger (the B buffer) is at height 0
 *
 ***********************************************************************************/


#include <iostream>
#include <cstdlib>
#include <math.h>
#include <stack>
#include <queue>
#include "FunnelHeapWithMerging.h"

#include <assert.h>

using namespace std;


// Enable assertions to ensure that the functions are operating correctly on a technical level
#define ASSERTIONS true

// Enable the PROFILE_ADAPTIVE_PERF to optimize the insertion (when a sweep is needed)
#define PROFILE_ADAPTIVE_PERF false

/* Function Prorotypes */
inline void appendBufferToStream(Buffer* buff, stream_t &stream);
inline monom_t pop (stream_t& stream1, stream_t& stream2, buf_size_t &numOfElementsMerged);
inline monom_t top (stream_t& stream1, stream_t& stream2);
inline bool empty(stream_t& stream1, stream_t& stream2);
inline void checkStream(stream_t &stream);
void checkBufferInvariant(Buffer* b);
void printStream(stream_t &stream);
#define LINK_OVERHEAD (int) 0       // number of elements in beginning of a link or
#define BUFFER_OVERHEAD (int) 0     // overhead amount associated with each buffer




#define I_INDEX 4294967295   //(long) pow(2, 31*8) - 1    // computes the index of the buffer I


#pragma mark -
#pragma mark Function Prototypes

/* Function Prototypes */



// computations of sizes and capacities
inline void getKAndS (long values [], int i);
inline unsigned long long computeLinkRequiredStorage(int i);



#pragma mark -
#pragma mark Header Functions

/**
 * 	A constructor that initializes the heap with the specified number of links
 *
 *	@param numOfLinks the number of links
 */
FunnelHeapWithMerging::FunnelHeapWithMerging (int numberOfLinks) {  // "numOfLinks" should be derived depending on input size (DO LATER)
	numOfLinks = numberOfLinks;
	
	_size = 0;
	// create link pointers
	links = new Link* [numberOfLinks];

	
	
    const int s1 = 8;
    
    
	// compute the required storage amount for the entire heap
    int I_size = s1 + BUFFER_OVERHEAD;       // size of I is s1 elements and the buffer overhead
    
    unsigned long long storage = I_size + 2;		// size of buffer I + 2 (for size of queue and I_NULL element)
    
    for (int i = 1; i <= numOfLinks; i++) {
        storage += computeLinkRequiredStorage(i) + LINK_OVERHEAD;
    }
    
    
	if (verboseLevel(VERBOSE_LOW))
		printf("Queue size (in MB) = %llu\n", (monom_t) ((storage / ( 1024.0 * 1024.0) * sizeof(monom_t) )));
														
	// reserve space for the heap
	heap.reserve(storage);
		
	// setup the heap
    setupHeap();  


	// INITIALIZE HEAP STATISTICS
	STAT_linkSweeps = new unsigned long long [numOfLinks];
	for (int i = 0; i < numOfLinks; i++) {
		STAT_linkSweeps[i] = 0;
	}
}


bool FunnelHeapWithMerging::isEmpty() {
//	return _size == 0;
	return I->isEmpty() && getLink(1)->A->isExhausted();
}


monom_t FunnelHeapWithMerging::peek() {
    
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

void FunnelHeapWithMerging::insert (deg_t degree, coef_t coef, ID_t id) {
	insert (degree, coef);
}


void FunnelHeapWithMerging::insert (deg_t degree, coef_t coef) {
    
    /*
     * General Code Sketch:
     *
     *  1)  Start by inserting the new element in the buffer I.
     *  2)  If I becomes full, then sweep to free up space for
     *      new elements to be inserted in the heap
     */
	
	bool merged = false;		// indicates if element was merged in I (ie dont increase size)

    // insert new element in I
	I->insertSortedWithMerging(createMonomial(degree, coef, 0), merged);
	
	// increment the size of the heap
	
	if (!merged) _size++;
    
	
    // if I is full, then return. Otherwise, move on to the sweeping
    if (!I->isFull()) {
        return;
    }

    
    // locate empty leaf to sweep to
    
    // get first empty leaf
    long c_i = 0;   
    
    
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
        
		long k_i = getK(i);
        
        if (c_i <= k_i) {      // then we found the first empty buffer
            break;
        }
    
        // go to next link (if it exists)

        if (i >= numOfLinks) {           // if there are no more link, don't insert
            printWarning("Queue is full. Element cannot be inserted!");
            return;
        }
		
        link_i = getLink(i+1);
		i++;
    }

	// cout << "Sweeping to link " << i << "..." << endl;
	// print();
	// pausePrompt();
    sweep(i);
	// if (verboseLevel(VERBOSE_HIGH)) cout << "Sweeped successfully" << endl << endl;
	
	// print();
	
	updateStatistics(i);
	
	// cout << "1" << endl;
	// if (lamp == 34155) checkInvariant();
	// cout << "2" << endl << endl;
}


size_t lamp = 0;

monom_t FunnelHeapWithMerging::poll () {
	assert(_size > 0);			// sanity check
	
	++lamp;
	
	buf_size_t numberOfElementsRemoved = 0;
	monom_t elem = pollInternal(numberOfElementsRemoved);  
	
	assert(numberOfElementsRemoved <= _size);		// sanity check
	_size -= numberOfElementsRemoved;
	
	return elem;
}


heap_size_t FunnelHeapWithMerging::size () {
	return _size;
}

#pragma mark -
#pragma mark Initilization of Queue (sub methods)



/*
 * setQueueDefaultValues: sets all the default values of links and buffers in their correct
 *                          location in the heap structure
 */
void FunnelHeapWithMerging::setupHeap () {
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
Link* FunnelHeapWithMerging::getLink(int linkIndex) {
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
#pragma mark Sweep, Fill...


/**
 *	Forms the first stream of the sweep method; this involves extracting elements from the
 *	link we're "sweeping" to, starting from the A-buffer on a path down to the parent of the
 *	leaf we're sweeping to.
 *
 *	@param stream the stream to be filled
 */
void FunnelHeapWithMerging::formStream1 (stream_t &stream1, int linkIndex, vector<buf_size_t> &stream1BufferSizes) {
	// if (verboseLevel(VERBOSE_HIGH)) cout << "Forming stream 1" << endl;
	Link* link_sweep = getLink(linkIndex);
	
	// start with A_sweep
	Buffer* currBuff = link_sweep->A;		// point to A_sweep
	stream1BufferSizes.push_back(currBuff->getSize());

	// append elements to stream1
	appendBufferToStream(currBuff, stream1);

	// move on to B_sweep
	currBuff = currBuff->getLeftChild();		// point to B_sweep
	
	// traverse the path from B_sweep to S_ci, appending the buffers to stream1 and recording their sizes.
	int currC = link_sweep->c;
	long currK = getK(linkIndex);
	
	// loop until we reach S_ci
    while (currBuff->getType() != LEAF_BUF_T) {
		stream1BufferSizes.push_back(currBuff->getSize());
		appendBufferToStream(currBuff, stream1);
		step(currBuff, currC, currK);
	}
	
}


/**
 *	Forms the second stream of the sweep method; this performs an internal extraction (i.e. does 
 *	not affect the size parameter of the heap). Having already marked A_sweep as exhausted, this 
 *	extraction will move all elements in the links before link_sweep to stream2.
 */
void FunnelHeapWithMerging::formStream2 (stream_t &stream2) {
	// if (verboseLevel(VERBOSE_HIGH)) cout << "Forming stream 2" << endl;
    // Forming stream2: 	
	
	buf_size_t numOfElementsRemoved = 0;		// indicates the number of elements removed during a pollInternal call
	
	// while I is not empty, or A_1 is not exhausted
	while (!(I->isEmpty() && getLink(1)->A->isExhausted())) {
		stream2.push(pollInternal(numOfElementsRemoved));
	}
	
}


void FunnelHeapWithMerging::reinsertStreams(stream_t &stream1, stream_t &stream2, int sweepIndex, vector<buf_size_t> &ASizes, vector<buf_size_t> &stream1BufferSizes) {
	// if (verboseLevel(VERBOSE_HIGH)) cout << "Re-inserting elements" << endl;
	// Re-inserting elements: Conceptually, the two streams should be merged, and then reinserted into the heap.
	// But for performance, the streams will not be merged together into a thrid stream, rather they will act 
	// together as one virtual stream, meanining that when we need an element from this virtual (third) stream, 
	// we take the maximum element of the two streams.
	// The insertion involves inserting elements into the A-buffers of the links before link_sweep, preserving
	// the sizes they had before this sweep operation. Then, the remaining elements are stored in the path
	// from A_sweep to S_ci (also preserving previous sizes), and the elements that remain are stored in S_ci
	//
	// NOTE: The whole duplicate-merging process (i.e. merging together elements of equal degrees while in the heap)
	// may cause a break in the heap invariant if not dealt with cautiously. Specifically, merging duplicates might
	// place "lower-ranked" elements in buffers in which they do not belong. This is because we have previously 
	// recorded the sizes of the buffers before the sweep method, and are refilling these buffers with that same 
	// amount of elements. In particular, this can happen when dealing with buffers that are in the link to which
	// we are sweeping to. To solve this, the method responsible for refilling the buffers with elements after
	// they have been removed (i.e. the refillBuffer method) will perform a check to make sure its not adding
	// an element to the buffer that is smaller than the largest element in one of its children. If it does reach
	// such a case, it will stop filling that buffer, and will move on to the next one
	
	Link* currLink;
	
	// start by filling the A-buffers of the links before link_sweep
	for (int i = 1; i < sweepIndex; i++) {
		if (empty(stream1, stream2)) 
			return;
		
		currLink = getLink(i);
		refillBuffer(currLink->A, stream1, stream2, ASizes[i-1]);
		currLink->A->setExhausted(false);			// set not exhausted
	}

	Link* link_sweep = getLink(sweepIndex);

	if (empty(stream1, stream2)) 
		return;
		
	// fill A_sweep
	refillBuffer(link_sweep->A, stream1, stream2, stream1BufferSizes[0]);
	link_sweep->A->setExhausted(false);

	Buffer* currBuff;
	int currC;
	long currK;
	
	// fill the path from B_sweep to S_ci
	currBuff = link_sweep->A->getLeftChild();		// points to B_sweep
	currC = link_sweep->c;
	currK = getK(sweepIndex);

	for (size_t i = 1; i < stream1BufferSizes.size(); i++) {
		if (empty(stream1, stream2)) 
			return;
			
		refillBuffer(currBuff, stream1, stream2, stream1BufferSizes[i]);
		currBuff->setExhausted(false);
		
		step(currBuff, currC, currK);
	}


	if (empty(stream1, stream2)) 
		return;
		
	// // fill the leaf with the remaining elements
	refillBuffer(currBuff, stream1, stream2, stream1.size() + stream2.size());
	currBuff->setExhausted(false);

}


/*
 *  sweep: used to empty buffer I and make room for new elements to be inserted
 *
 *  @param q: pointer to the q 
 *  @param linkIndex: the index of the link where the values in I will be "sweeped" into
 */
void FunnelHeapWithMerging::sweep (int linkIndex) {
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
	assert (linkIndex >= 1 && linkIndex <= numOfLinks);		// sanity check
	
	// form stream 1
	stream_t stream1;
	vector<buf_size_t> stream1BufferSizes;
	
	formStream1(stream1, linkIndex, stream1BufferSizes);

	// store all sizes of the A-buffers of all links up to (BUT EXCLUDING) A_sweep
	vector<buf_size_t> ASizes;
	Link* currLink;
	for (int i = 1; i <= linkIndex-1; i++) {
		currLink = getLink(i);
		ASizes.push_back(currLink->A->getSize());
	}

	Link* link_sweep = getLink(linkIndex);
    // set A_sweep as exhausted (possibly break invariant temporarily, intentionally);
	bool A_exhausted = link_sweep->A->isExhausted();
	link_sweep->A->setExhausted(true);
	
	// form stream 2
	stream_t stream2;
	formStream2(stream2);
	
	link_sweep->A->setExhausted(A_exhausted);		// set the exhausted value of A_sweep back to what it was
	
	
	// reinsert the two streams back into the heap
	reinsertStreams(stream1, stream2, linkIndex, ASizes, stream1BufferSizes);
	
    // reset c_i and r_i (for i < linkindex) to 1
    for (int i = 1; i <= linkIndex-1; i++) {
		currLink = getLink(i);
		currLink->c = 1;
    }
    
    // increment c_i (for i = linkIndex) by 1
	link_sweep->c = link_sweep->c + 1;
	
}


/* 	peek: returns the maximum element from the two streams, without
 *			actually extracting it or performing merging of any kind
 */
inline monom_t top (stream_t& stream1, stream_t& stream2) {
	
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
 *	@param qSize: points to the size of the queue
 */
inline monom_t pop (stream_t& stream1, stream_t& stream2, buf_size_t &numOfElementsMerged) {
    /*
     * General Code Sketch:
     *
     *  1) Extract the element that has the highest degree in both streams
     *  2) Return it
     */
    
    if (stream1.empty() && stream2.empty())
        return 0;
    
	bool chosenFromStream1 = false;		// true if element is chosen from stream1
	
    // Find the max element in both streams
    monom_t maxElement = 0;
    if (!stream1.empty())
        maxElement = stream1.front();
    
    if (!stream2.empty()) {
        if (GET_DEGREE(stream2.front()) >= GET_DEGREE(maxElement)) {
            // the max is from stream2
            maxElement = stream2.front();
            stream2.pop();
			chosenFromStream1 = false;
        }
        else {
            // the max is from stream1
            stream1.pop();
			chosenFromStream1 = true;
        }
    }
    else {      // the max was from stream1 (since stream2 is empty)
        stream1.pop();
		chosenFromStream1 = true;
    }

	// now that we removed the max element, check if there are any more
	// elements of that same degree in any of the two streams. If found,
	// remove them from the streams and merge them with the current max
	// (i.e. add up their coefficients). As elements are being merged 
	// together, remember to decrease the overall size of the queue
	
	deg_t maxDegree = GET_DEGREE(maxElement);
	coef_t maxCoef = GET_COEF(maxElement);
	
	numOfElementsMerged = 0;
	
	// check if there are any 'max' elements in stream1
	while (!stream1.empty()) {
		// if it exists, then merge it 
		if (GET_DEGREE(stream1.front()) == maxDegree) {
			maxCoef += GET_COEF(stream1.front());	// add up the coefs
			stream1.pop();		// remove it from the stream
			numOfElementsMerged++;
		}
		else {
			break;
		}
	}
	
	
	// check if there are any 'max' elements in stream2
	while (!stream2.empty()) {
		// if it exists, then merge it 
		if (GET_DEGREE(stream2.front()) == maxDegree) {
			maxCoef += GET_COEF(stream2.front());	// add up the coefs
			stream2.pop();		// remove it from the stream
			numOfElementsMerged++;
		}
		else {
			break;
		}
	}
	
	//maxElement = packDC(maxDegree, maxCoef);		// re-pack the degree and coef	
	maxElement = createMonomial(maxDegree, maxCoef, 0);
	
    return maxElement;      
}

inline bool empty(stream_t& stream1, stream_t& stream2) {
    return stream1.empty() && stream2.empty();
}



#pragma mark -
#pragma mark profile_adaptive_performance


inline void appendBufferToStream(Buffer* buff, stream_t &stream) {
	// append the buffer to the end of the stream
	while (!buff->isEmpty()) {
		stream.push(buff->pop_front());
	}
}

unsigned long long counter = 0;

void FunnelHeapWithMerging::refillBuffer (Buffer* buff, stream_t &stream1, stream_t &stream2, buf_size_t numberOfElements) {
	// make sure the buffer can fit <numberOfElements> and the stream has enough elements
	assert(numberOfElements <= buff->getCapacity());
	assert(buff->isEmpty());
	
	// NOTE: due to the fact that we are merging elements of equal degrees while in the heap, it is possible
	// that the <numberOfElements> to be inserted into the current buffer will be greater than the actual 
	// remaining number of elements contained in the two streams. Therefore, we add a check that will stop
	// inserting elements into the buffer if the streams are empty. 
	// Also, the merging process adds the possiblity of breaking a heap invariant, where we insert an element
	// in the current buffer whose degree is smaller than the degree of the largest element in one of its children.
	// We therefore stop inserting elements at that point as well
	
	buff->checkInvariant();
	
	buf_size_t numOfElementsMerged = 0;
	buf_size_t tempMerged = 0;
	for (buf_size_t i = 0; i < numberOfElements; i++) {
		if (empty(stream1, stream2))
			break;
		
		monom_t elemDegree = GET_DEGREE(top(stream1, stream2));
		Buffer* left = buff->getLeftChild();
		Buffer* right = buff->getRightChild();
		
		if (left != NULL && !left->isExhausted()) {
			if (elemDegree < GET_DEGREE(getLargestElementInChild(left))) {
				break;
			}
		}
		if (right != NULL && !right->isExhausted()) {
			if (elemDegree < GET_DEGREE(getLargestElementInChild(right))) {
				break;
			}
		}
			
		monom_t popped = pop(stream1, stream2, tempMerged); 
		assert (GET_DEGREE(popped) == elemDegree);
		buff->push_back(popped);		
		numOfElementsMerged += tempMerged;
		
	}
	

	_size -= numOfElementsMerged;

	buff->checkInvariant();
}



monom_t FunnelHeapWithMerging::pollInternal(buf_size_t &numOfElementsRemoved) {
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

	// indicate the number of elements that are removed from the heap
	numOfElementsRemoved = 0;

	monom_t popped = 0;
	// if neither is empty, then take the max element from both
    if (!I_empty && !A_1_empty) {
		// if I and A_1 contain the max, merge the elements together
		if (GET_DEGREE(I->front()) == GET_DEGREE(A_1->front())) {
			// extract the top elements that have the same key in A_1 as I

			// first start by extracting the first element from A_1 and I
			deg_t maxDegree = GET_DEGREE(I->front());
			coef_t maxCoef = GET_COEF(I->front()) + GET_COEF(A_1->front());
			I->pop_front();
			A_1->pop_front();

			// indicate that elements have been removed
			numOfElementsRemoved = 2;

			// fill A_1 if its empty
			if (A_1->isEmpty()) 
				A_1->fill();


			popped = createMonomial(maxDegree, maxCoef, 0);
		}
		else if (GET_DEGREE(I->front()) >= GET_DEGREE(A_1->front())) {		// if I contains the max element
			popped = I->pop_front();
			numOfElementsRemoved = 1;
    	}
		else {		// if A_1 contains the max element
			popped = A_1->pop_front();
			numOfElementsRemoved = 1;
		}
	}
	// if I is empty but A_1 isn't, take max from A_1
	else if (I_empty && !A_1_empty) {
		popped = A_1->pop_front();
		numOfElementsRemoved = 1;
	}
	// if A_1 is empty but I isn't, take max from I
	else {
		popped = I->pop_front();
		numOfElementsRemoved = 1;
	}


	deg_t maxDegree = GET_DEGREE(popped);
	coef_t maxCoef = GET_COEF(popped);

	// NOTE: There might be (consecutive) elements in A_1 having
	// that max key. These are not merged because they might be sitting
	// in a leaf that is not affected by the sweep method at some point 
	// in time

	if (A_1->isEmpty()) 
		A_1->fill();

	// look for remaining elements in A_1 of that degree
	while (!A_1->isExhausted() && GET_DEGREE(A_1->front()) == maxDegree) {
		maxCoef += GET_COEF(A_1->front());
		A_1->pop_front();

		// indicate that an element has been removed
		numOfElementsRemoved++;

		// make sure A_1 is not empty
		if (A_1->isEmpty()) 
			A_1->fill();
	}

	popped = createMonomial(maxDegree, maxCoef, 0);

	return popped;
}



void checkBufferInvariant(Buffer* b) {
	if (b->isEmpty())
		return;
		
	Buffer *left = b->getLeftChild();
	Buffer *right = b->getRightChild();
	Buffer *parent = b->getParent();
	
	// check left child
	if (left != NULL) {
		if (!left->isEmpty()) {
			assert (GET_DEGREE(b->back()) >= GET_DEGREE(left->front()));
		}
	}
	// check right child
	if (right != NULL) {
		if (!right->isEmpty()) {
			assert (GET_DEGREE(b->back()) >= GET_DEGREE(right->front()));
		}
	}
	
	// check parent
	if (parent != NULL) {
		if (!parent->isEmpty() && parent->getType() != I_BUF_T) {
			assert (GET_DEGREE(b->front()) <= GET_DEGREE(parent->back()));
		}
	}
	
}



// PRINTING THE HEAP
void FunnelHeapWithMerging::print() {
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



/**
 *	Updates the statistics of the heap
 *
 *	@param sweepIndex the index of the link at to which the sweep operation occurred
 */
inline void FunnelHeapWithMerging::updateStatistics(int sweepIndex) {
	assert (sweepIndex >= 1 && sweepIndex <= numOfLinks);		// sanity check
	
	if (STAT_linkSweeps != NULL) {
		STAT_linkSweeps[sweepIndex-1]++;		// increment the sweep 
	}
}






void FunnelHeapWithMerging::checkInvariant() {
	for (int i = 1; i <= numOfLinks; i++) {
		getLink(i)->checkInvariant();
	}
}



inline void checkStream(stream_t &stream) {
	if (stream.empty())
		return;
		
	stream_t s (stream);

	monom_t elem1 = s.front();
	s.pop();
	while (!s.empty()) {
		assert (GET_DEGREE(elem1) >= GET_DEGREE(s.front()));
		elem1 = s.front();
		s.pop();
	}
}



inline void printStream(stream_t &stream) {
	stream_t s (stream);
	cout << "Stream (size = " << s.size() << "): [";
	
	size_t limit = s.size();
	for (size_t i = 0; i < limit; i++) {
		if (i > 0) {
			cout << " ";
		}
		
		cout << GET_DEGREE(s.front()) << " ";
		s.pop();
	}
	cout << "]" << endl;
}

