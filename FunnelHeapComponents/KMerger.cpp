#include "KMerger.h"
#include <iostream>
#include "../Options/Options.h"

using namespace std;
/**
 *	Constructor: Initializes a K-merger of the corresponding height.
 *	
 *	@param mem a pointer to the memory location where the k-merger will be placed
 *	@param k the number of input streams that the K-merger can merge
 *	@param outputBuffer the buffer to which the K-merger will send its merged output
 */
KMerger::KMerger(mem_t& mem, long k, Buffer *outputBuffer) {
	k_value = k;		
	output = outputBuffer;
	
	// cout << "Creating a K-merger with K = " << k << endl << endl;
	setupKMerger(k, mem, output);
	
}


void postOrderFree(Buffer* n) {
	if (n == NULL)
		return;
	
	// do not free leafs
	if (n->isLeaf())
		return;

	postOrderFree(n->getLeftChild());		// free left child
	postOrderFree(n->getRightChild());		// free right child
	
	// free buffer if it is NOT a B-buffer
	if (n->getType() != B_BUF_T)
		delete n;			// free me
}


KMerger::~KMerger() {
	// perform postorder traversal of the tree and free every buffer (except B-buffer and leafs)
	postOrderFree(output);
}

/**
 *	Adds the input sources to the K-merger
 *
 *	@param S an array of pointers to input buffers
 *	@param numberOfInputSources the number of input sources
 */
void KMerger::addInputSources(Buffer **S, long numberOfInputSources) {
	assert (numberOfInputSources == k_value);			// sanity check
	
	// cout << "Adding Input Sources to K-merger with K = " << k_value << endl;
	// pausePrompt();
	
	// go to each buffer of the last level of the K-merger, and add its
	// children from S
	
	long numOfBuffersInLastLevel = (long) (numberOfInputSources/2);
	
	long inputCounter = 0;
	for (long i = 0; i < numOfBuffersInLastLevel; i++) {
				
		// setup variables for the step function
		Buffer* currBuff = output;		// points to B
		int destination = i+1;
		long currK = numOfBuffersInLastLevel;

		
		// if (k_value == 8) cout << "T: " << currBuff->getType() << ", d: " << destination << ", cK: " << currK<< endl;
		
		
		// go the the <destination>'th buffer at the last level of the K-merger
	 	while (currBuff->getLeftChild() != NULL) {
		// while (currK > 1){
			step(currBuff, destination, currK);
			// if (k_value == 8) cout << "T: " << currBuff->getType() << ", d: " << destination << ", cK: " << currK<< endl;
		}
		// if (k_value == 8) cout << endl;
		
		currBuff->setLeftChild(S[inputCounter]);	// setup its left child
		S[inputCounter]->setParent(currBuff);		// set parent of left child
		inputCounter++;
		
		currBuff->setRightChild(S[inputCounter]);		// setup right child
		S[inputCounter]->setParent(currBuff);		// set parent of right child
		inputCounter++;
	}
	
}





/****** Private Methods ******/

/*
 *	Sets up the K-merger with all its buffers.
 *
 * 	@param k the number of input streams that the K-merger can merge
 *	@param mem a pointer to the memory location where the subtree will be placed
 *  @param parent a pointer to the output buffer of the K-merger      
 */
void KMerger::setupKMerger (long k, mem_t& mem, Buffer *parent) {
	
	assert (k >= 2);		// sanity check
    
	if (getGlobalBoolean()) cout << "k = " << k << endl;
	
    if (k == 2) {           // no mid buffers to allocate
		return;
    }
    
    if (k == 4) {       
		// Only two "mid" buffers are needed; allocate their memory
		
		// get the capacities of the buffers
        buf_size_t bufferCapacity = (buf_size_t) ceil(pow(k, 1.5));
		

		// create first buffer, and set it as a child to its parent
		Buffer *buffer1 = new Buffer (mem, bufferCapacity, INTERIOR_BUF_T);
		parent->setLeftChild(buffer1);
		
		// STATE: mem points to the empty location right after the space for the first buffer
		
		// create second buffer, and set it as a child to its parent
		Buffer *buffer2 = new Buffer (mem, bufferCapacity, INTERIOR_BUF_T);
		parent->setRightChild(buffer2);

        // STATE: mem points to the empty location right after the space for the second buffer
		buffer1->setParent(parent);
		buffer2->setParent(parent);
		
		if (getGlobalBoolean()) cout << "Setup up two \"mid\" buffers" << endl;
		return;
    }
    

    double heightTemp = (log(static_cast<long double>(k))/log(static_cast<long double>(2)));
	long height = static_cast<long>(floor(heightTemp + 0.5));
	
		//     if (k >= 128)
		// height++;
		
    long top_i = (long) ceil(height/2.0);
    long bottom_i = height - top_i;
    
    long top_k = (long) pow(2.0, (double) top_i);      // get the k of the top subtree
    long bottom_k = (long) pow(2.0, (double) bottom_i);    // get the k of the bottom subtrees
    
    buf_size_t midBufferSize = (buf_size_t) ceil(pow(k, 1.5));
   
 if (getGlobalBoolean()) {
		cout << "height = " << height << endl;
		cout << "top_i = " << top_i << endl;
		cout << "bottom_i = " << bottom_i << endl;
		cout << "top_k = " << top_k << endl;
		cout << "bottom_k = " << bottom_k << endl;
		cout << "midBufferSize = " << midBufferSize << endl << endl;
		
		pausePrompt();
	
		cout << "----> setting up k-merger with k = " << top_k << endl;
	}
	
	// Reserve space for the top merger
	setupKMerger(top_k, mem, parent);
	
	if (getGlobalBoolean()) {
		cout << "<---- done setting k-merger with k = " << top_k << endl << endl;
		pausePrompt();
		cout << "Creating middle buffers..." << endl;
	}
	
	// Create the middle buffers
	Buffer* midBuffers[top_k];
	
	Buffer* currParent = parent;
	long topHeight = (long) (log(static_cast<long double>(top_k))/log(static_cast<long double>(2)));
	for (long i = 0; i < top_k; i++) {
		midBuffers[i] = new Buffer (mem, midBufferSize, INTERIOR_BUF_T);
		
		// get parent of current mid-buffer
		long currK = top_k;
		int currC = i+1;
		
		if (i%2 == 0) {
			currParent = parent;
			for (long j = 0; j < topHeight-1; j++) {
				step(currParent, currC, currK);
			}
		}
		
		// set left/right child of parent
        if (i % 2 == 0) {
            // set left child   
            currParent->setLeftChild(midBuffers[i]);     
        }
        else if (i % 2 == 1) {
            // set right child and go to next parent
            currParent->setRightChild(midBuffers[i]);
        }

		midBuffers[i]->setParent(currParent);
	}

	if (getGlobalBoolean()) {
		cout << "Done creating middle buffers" << endl << endl;
		pausePrompt();
		cout << "Creating bottom mergers..." << endl;
	}

	// Create bottom mergers
	for (long i = 0; i < top_k; i++) {
		parent = midBuffers[i];
		if (getGlobalBoolean())  cout << "----> setting up k-merger with k = " << bottom_k << endl;
		
		setupKMerger(bottom_k, mem, parent);
		
		if (getGlobalBoolean()) {
			cout << "<---- done setting k-merger with k = " << bottom_k << endl << endl;
			pausePrompt();
		}
	}
}


/****** Static Methods ******/

// Computes the size needed for buffers in a k-merger. 
int kSize(int k) {

  if (k<=2) {
    return 0; // No buffers exist in a 2-way merger
  } 

  double height = log(static_cast<long double>(k))/
                  log(static_cast<long double>(2));

  int heightInt = static_cast<int>(floor(height + 0.5));

  int bottomI = static_cast<int>(heightInt >> 1);
  int topI = heightInt - bottomI;

  int topK = 1 << topI;                 // = 2^(topI)
  int bottomK = 1 << bottomI;        // = 2^(bottomI) 

  int sizeOfBuffers = static_cast<int>(ceil(pow(k,1.5)));


  int totalMidSize = topK * sizeOfBuffers;

  // The size is the middlebuffers + the size of the top-merger + the
  // size of the topK bottomMergers.
  return totalMidSize + kSize(topK) + (topK * kSize(bottomK));
}


void step(Buffer*& pos, int& currC, long& currK) {
	assert(pos->getType() != A_BUF_T);
	assert(pos->getType() != I_BUF_T);		// make sure <pos> is not I or an A-buffer

	// Find the child on the path.
	if (currC <= (long)(currK / 2)) {
		// go one step to the left
		pos = pos->getLeftChild();
		currK = (long) (currK/2);		// take k of the left subtree
	} 
	
	else { 	
		// go one step to the right
		pos = pos->getRightChild();
		currK = (long) (currK/2);		// take k of the right subtree
		currC = currC - currK;			// modify c to only consider the right subtree
	}
} 




static int constantK[7] = {2,  4,   8,   16,    32,    128,      512};

void KMerger::print() {
	long numOfBuffersToPrint = (k_value/2) + ((k_value/2) - 2);
	
	if (numOfBuffersToPrint == 0)
		return;
		
	int destination = 1;		
	long level = 1;		// current level, where level 0 contains B-buffer
	long K = constantK[level-1];		// 2^level
	
	cout << "Level " << level << ": K = " << K << endl;
	
	for (long i = 0; i < numOfBuffersToPrint; i++) {
				
		// setup variables for the step function
		Buffer* currBuff = output;		// points to B
		int currDestination = destination;
		long currK = K;
		
		// move to requested buffer
		for (long level_count = 1; level_count <= level; level_count++) 
			step(currBuff, currDestination, currK);
		
		cout << "Buffer " << destination << ": "; currBuff->print(); cout << endl;
		
		// if we printed all buffers on current level, move down one level
		if (destination == K) {
			destination = 1;
			level++;
			K = 2*K; //constantK[level-1];
			
			if (i < numOfBuffersToPrint-1)
				cout << "Level " << level << ": K = " << K << endl;
		}
		// if not all buffers have been printed on current level, continue
		else {
			destination++;
		}
	}
}



void postOrderCheck(Buffer* n) {
	if (n == NULL)
		return;
	
	// do not check leafs
	if (n->isLeaf()) // || n->isExhausted())
		return;

	postOrderCheck(n->getLeftChild());		// check left child
	postOrderCheck(n->getRightChild());		// check right child
	
	// do not check B-buffer
	if (n->getType() != B_BUF_T) {
		n->checkInvariant();			// check me
	}
}


void KMerger::checkInvariant() {
	postOrderCheck(output);
}


