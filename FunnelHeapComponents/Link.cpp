#include "Link.h"

// Tables of S and K
int constantS[7] = {8, 24, 120, 1080, 18630, 605880, 78158520};
int constantK[7] = {2,  4,   8,   16,    32,    128,      512};

Link::Link (mem_t& mem, int linkIndex) {
	
	index = linkIndex;
	
	long k = constantK[index-1];
	
	
	c = 1;
	
	// compute the capacities of the Buffer A and B
	buf_size_t bufferSize = (buf_size_t) (pow(static_cast<long double>(k),3));
	
	// allocate space for A
	A = new Buffer (mem, bufferSize, A_BUF_T);
	
	// allocate space for B
	B = new Buffer (mem, bufferSize, B_BUF_T);
	
	// set parent/children of A and B
	A->setLeftChild(B);
	B->setParent(A);
	
	// Create the K-merger
	merger = new KMerger(mem, k, B);

	// Create the K-merger input buffers	
	S = new Buffer* [k];
	
	buf_size_t S_size = constantS[index-1];
	for (long i = 0; i < k; i++) {
		S[i] = new Buffer (mem, S_size, LEAF_BUF_T);
	}
	
	// add the input buffers to the K-merger
	merger->addInputSources(S, k);
}


#include <iostream>
using namespace std;

void Link::print() {
	// print A 
	cout << "A: "; A->print(); cout << endl;
	
	// print B
	cout << "B: "; B->print(); cout << endl;
	
	// print the KMerger
	merger->print();
	
	// print the leafs
	long numOfLeafs = constantK[index-1];
	
	for (long i = 1; i <= numOfLeafs; i++) {
		cout << "Leaf " << i << ": "; S[i-1]->print(); cout << endl;
	}
}





long getK(int linkIndex) {
	assert(linkIndex >= 1);
	return constantK[linkIndex-1];
}

buf_size_t getS(int linkIndex) {
	assert(linkIndex >= 1);
	return constantS[linkIndex-1];
}


void Link::checkInvariant() {
	if (A->isExhausted())
		return;
		
	A->checkInvariant();
	B->checkInvariant();
	
	merger->checkInvariant();
	
	
	// check invariants of leafs
	for (long i = 0; i < getK(index); i++) {
		S[i]->checkInvariant();
	}
	
}