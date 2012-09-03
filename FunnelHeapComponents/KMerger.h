/*
 *	\class KMerger
 *	\brief A KMerger is a structure that merges K sorted input streams into one, preserving their sorted order
 *
 *	Implementation Note: Given the context that this KMerger will be used in and the requirement
 *	that all its components must be stored contiguously in memory, the KMerger is passed
 *	(at initialization) a memory location which represents the start of where this structure is located
 *	in memory. All subsequent components of this KMerger must be stored contiguously starting at the
 *	specified location.
 */

#ifndef KMERGER_H
#define KMERGER_H

#include <math.h>
#include <assert.h>

#include "../Options/Options.h"
#include "Buffer.h"



class KMerger {
private:
	long k_value;		/**< the value of k */
	
	Buffer* output;		/**< the output buffer, to which the K-merger utilmately merges all input streams to */
	
	/*
	 *	Sets up the K-merger with all its buffers.
	 *
	 * 	@param k the number of input streams that the K-merger can merge
	 *	@param mem a pointer to the memory location where the subtree will be placed
	 *  @param parent a pointer to the output buffer of the K-merger      
	 */
	void setupKMerger (long k, mem_t& mem, Buffer* parent);
	
public:
	
	/**
	 *	Constructor: Initializes a K-merger of the corresponding height.
	 *	
	 *	@param mem a pointer to the memory location where the k-merger will be placed
	 *	@param k the number of input streams that the K-merger can merge
	 *	@param outputBuffer the buffer to which the K-merger will send its merged output
	 */
	KMerger(mem_t& mem, long k, Buffer* outputBuffer);

	~KMerger ();

	/**
	 *	Adds the input sources to the K-merger
	 *
	 *	@param S an array of pointers to input buffers
	 *	@param numberOfInputSources the number of input sources
	 */
	void addInputSources(Buffer **S, long numberOfInputSources);
	
	/**
	 *	Prints the KMerger
	 */
	void print();
	
	void checkInvariant();
};


// Static Methods
// Computes the size needed for buffers in a k-merger. 
int kSize(int k);

// steps through a K-merger, one step at a time

void step(Buffer*& pos, int& currC, long& currK);


#endif