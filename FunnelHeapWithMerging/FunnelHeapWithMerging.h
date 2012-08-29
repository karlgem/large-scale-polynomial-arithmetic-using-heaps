//
//  Funnel Heap
//  Funnel Heap 
//
//  Created by Karl Gemayel on 16/05/2012.
//  Copyright 2012 American University Of Beirut. All rights reserved.
//

/** 
 * \class FunnelHeap
 * \brief A Max Heap implementation using the Funnel Heap design.
 *
 * The Funnel Heap is list of consecutive links, where each link contains its own K-Funnel,
 * whose properties depend on the index of the link. The K-Funnel is a binary tree-like 
 * structure whose nodes represent mergers, and edges and leafs represent buffers. Each buffer
 * stores a sorted sequence of elements, and therefore invoking a 'fill' operation on the 
 * root would fill that root with a sorted sequence of the elements present in the tree. 
 * Consecutive links are also connected together with mergers, which would result in a sorted
 * extraction from the whole queue when we wish to extract elements.
 */

#ifndef FunnelHeapWithMerging_h
#define FunnelHeapWithMerging_h

#include "../FunnelHeapComponents/Link.h"
#include "../Heap/Heap.h"
#include "../Polynomial/Monomial.h"
#include "../Options/Options.h"


typedef queue<monom_t> stream_t;



class FunnelHeapWithMerging : public Heap {
private:
	int numOfLinks;		/**< number of links */
	heap_t heap;		/**< the actual heap itself */
	
	heap_size_t _size;
	
	// links
	Link** links;
	
	// I-buffer
	Buffer* I;
	
	void setupHeap();
	
	Link* getLink(int linkIndex);
	
	void sweep (int linkIndex);
	
	monom_t pollInternal(buf_size_t &numOfElementsRemoved);
	
	void refillBuffer(Buffer* buff, stream_t &stream1, stream_t &stream2, buf_size_t numberOfElements);
	
	// HEAP STATISTICS
	
	// statistics variables
	unsigned long long *STAT_linkSweeps;		// array that indicates the number of sweeps that occured at each link
	
	
	/**
	 *	Updates the statistics of the heap
	 *
	 *	@param sweepIndex the index of the link at to which the sweep operation occurred
	 */
	inline void updateStatistics(int sweepIndex);
	
public:
	/**
	 * A constructor that initializes the heap with the specified number of links
	 *
	 * \param numOfLinks the number of links
	 */
	FunnelHeapWithMerging(int numOfLinks);
	
	/**
	 *	Destructor that deallocates memory for variables belonging to this class
	 */
	~FunnelHeapWithMerging() {
		delete I;		// delete the I-buffer
		
		// delete all links
		for (int i = 0; i < numOfLinks; i++) {
			delete links[i];
			links[i] = NULL;
		}
		
		delete links;
		
		
		// print and clear statistics
		if (STAT_linkSweeps != NULL) {
			unsigned long long totalNumberOfSweeps = 0;
			cout << "Heap Statistics:" << endl;
			for(int i = 1; i <= numOfLinks; i++) {
				cout << "\tLink " << i << ": " << STAT_linkSweeps[i-1] << " sweeps" << endl;
				totalNumberOfSweeps += STAT_linkSweeps[i-1];
			}
			
			cout << "\tTotal Number of Sweeps = " << totalNumberOfSweeps << endl;
			cout << endl;
			
			delete STAT_linkSweeps;
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
	virtual void insert (deg_t degree, coef_t coef);
	virtual void insert (deg_t degree, coef_t coef, ID_t id);
	
	/** 
	 * Extracts the max element of the heap and returns it.
	 *
	 * \return The value of the max element
	 */
	virtual monom_t poll ();
	
	/** 
	 * Returns the value of the max element of the heap without actually extracting it.
	 *
	 * @return The value of the max element 
	 */
	virtual monom_t peek();
	
	/** 
	 * Returns true if the heap is empty, and false otherwise
	 *
	 * \return if heap is empty - true. Otherwise - false
	 */
	virtual bool isEmpty();
	
	/**
	 * Returns the size of the heap.
	 *
	 * \return The size of the heap.
	 */
	virtual heap_size_t size ();
	
	/**
	 *	Prints the heap
	 */
	virtual void print();
};

#endif
