/*
 *	\class Link
 *	\brief A Link that holds buffers and a K-merger
 *
 *
 *	Implementation Note: Given the context that this link will be used in and the requirement
 *	that all its componenets must be stored contiguously in memory, the link is passed
 *	(at initialization) a memory location which represents the start of where this link is located
 *	in memory. All subsequent components of this link must be stored contiguously starting at the
 *	specified location.
 */

#ifndef LINK_H
#define LINK_H

#include "Buffer.h"
#include "KMerger.h"
#include "../Options/Options.h"


// static int constantS[7] = {8, 24, 120, 1080, 18630, 605880, 78158520};
// static int constantK[7] = {2,  4,   8,   16,    32,    128,      512};

class Link {
private:
	int index;		// the index of the link
	
	
	Buffer* B;
	
	KMerger* merger;
	
	Buffer** S;
	
public:
	
	Buffer* A;
	
	// the empty leaf counter
	unsigned long c;
	
	
	Link (mem_t& mem, int linkIndex);
	
	~Link () {
		delete A;
		delete B;
		delete merger;
		delete S;
	}
	
	/**
	 * Returns the index of the link
	 */
	int getIndex();
	
	/**
	 *	Prints the link
	 */
	void print();
};


long getK(int linkIndex);

buf_size_t getS(int linkIndex);

#endif