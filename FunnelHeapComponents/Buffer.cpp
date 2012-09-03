#include <assert.h>
#include "Buffer.h"
	
#include <iostream>
using namespace std;

/* Private Methods */
	
/*
 *	Increments the pointer in a circular manner (goes forward)
 */
void Buffer::increment (mem_t& p) {
	// if the pointer p has reached the end of the buffer, circulate
	if (p == _elemLoc + getCapacity() - 1) 
		p = _elemLoc;
		
	// otherwise, just increment
	else
		p++;
	
}

/*
 *	Decrements the pointer in a circular manner (goes backward)
 */

void Buffer::decrement (mem_t& p) {
	// if the pointer p has reached the start of the location where elements can be inserted, circulate
	if (p == _elemLoc)
		p = _elemLoc + (_capacity - 1);
	else
		p--;
}
	
	

/* Public Methods */

	
/**
 *	Constructor: initializes the buffer (at the given memory location) with the
 *	specified capacity.
 *
 *	@param mem the starting memory location at which to allocate the buffer's components
 *	@param cap the buffer's maximum capacity
 *	@param bType indicates the type of the buffer (I, Interior, Leaf)
 */

Buffer::Buffer(mem_t& mem, buf_size_t cap, int bType) {
	// initialize buffer variables
	_loc = mem;
	_elemLoc = mem;		// elemLoc points to the first memory location where elements are inserted
	
	_front = _elemLoc;
	_back = _elemLoc;
	
	// initialize buffer characteristics
	_size = 0;
	_capacity = cap;
	_bufferType = bType;
	_exhausted = true;
	
	_leftChild = NULL;
	_rightChild = NULL;
	_parent = NULL;
	
	// move memory pointer just after the memory space allocated for the buffer
	mem += 0 + _capacity;	
}
	
/**
 *	Returns the buffer's capacity
 *	@return capacity
 */

buf_size_t Buffer::getCapacity() {
	return _capacity;
}
	
/**
 *	Returns the buffer's size
 *	@return size
 */

buf_size_t Buffer::getSize() {
	return _size;
}
	
/**
 *	Checks if the buffer is empty
 *
 *	@return returns true if the buffer is empty; false otherwise
 */

bool Buffer::isEmpty() {
	return (_size == 0);
}

/**
 *	Checks if the buffer is full
 *
 *	@return returns true if the buffer is full; false otherwise
 */
bool Buffer::isFull() {
	return _size == _capacity;
}
	
/**
 *	Checks if the buffer is exhausted (i.e. its children and itself are empty)
 *
 *	@return returns true if the buffer is exhausted; false otherwise
 */
bool Buffer::isExhausted() {
	return _exhausted;
}


/**
 *	Returns the type of the buffer
 *
 *	@return type
 */
int Buffer::getBufferType() {
	return _bufferType;
}
	
	
/**
 *	Inserts an element to the end of the buffer
 *	
 *	@param element inserted element
 */
void Buffer::push_back(monom_t element) {
	assert(!isFull());		// sanity check
	if (isFull()) {
		cout << _bufferType << " full " << endl;
		pausePrompt();
		return;
	}
		
	if (_size == 0) {
		*_back = element;
	}
	else {
		
		decrement(_back);		// move the back pointer one step to the back to make space for the element
//		if (isLeaf()) cout << "o " << _back - _elemLoc << "/" << _capacity << endl;
		*_back = element;		// add the element to the end of the list of elements		
//		if (isLeaf()) cout << "a " << endl;
	}
	
	_size++;					// increment the size by 1
	
	// the buffer is not exhausted
	_exhausted = false;
	
//	checkInvariant();
}

/**
 *	Extracts the first element of the buffer
 *
 *	@return the first element
 */
monom_t Buffer::pop_front() {
	assert(!isEmpty());		// sanity check
	monom_t popped = *_front;		// get the first element
	decrement(_front);		// move the front pointer one step to the back to "delete" the front element
	_size--;					// decrement the size by 1
	
	// check if the buffer became exhausted, and make sure front and back pointers are not wrapped
//	if (isEmpty()) {
	if (_size == 0) {
		_front = _elemLoc;
		_back = _elemLoc;
		
		if (isLeaf()) 
			_exhausted = true;
		else {
			if (_leftChild->isExhausted()) {
				if (_rightChild == NULL)
					_exhausted = true;
				else if (_rightChild->isExhausted())
					_exhausted = true;
			}
		}
	}
	
//	checkInvariant();
	
	return popped;
}

/**
 *	Peeks at the front element without actually extracting it
 *
 *	@return a copy of the front element
 */
monom_t Buffer::front() {
	assert(!isEmpty());		// sanity check
	return *_front;
}

/**
 *	Peeks at the back element without actually extracting it
 *
 *	@return a copy of the back element
 */
monom_t Buffer::back() {
	assert(!isEmpty());
	return *_back;
}

/**
 *	Inserts the element into a sorted buffer, while maintaining its sorted order
 *
 *	@param element the element to be inserted
 */
void Buffer::insertSorted (monom_t element) {
	assert(_bufferType == I_BUF_T);			// sanity checks
	assert(!isFull());		
	
	if (_size == 0) {
		push_back(element);
		return;
	}
	
	mem_t pos = _front;
	
	buf_size_t counter = 0;
	
	monom_t monomial = (monom_t) element;
	// get attributes
	deg_t elemDegree = GET_DEGREE(monomial);
	
	// loop over all elements to find correct place for new element
	while (counter < _size) {
		// if the degree of the element is greater than <pos>, move up the list
 		if (elemDegree < GET_DEGREE(*pos)) {
			decrement(pos);
			counter++;
		}
		// if the degree of the element is greater or equal to <pos>, break;
		else {
			break;
		}
	}
	
	// if we reached the end of the list, append to the end
	if (counter == _size) {
		push_back(element);
		return;
	}
	
	// otherwise, insert element at position and shift elements downward
	decrement(_back);		// make room for new element
	mem_t shifter = _back;		// used to shift elements downward
	mem_t position = _back;		// used to indicate where the element should be shifted
	
	while (shifter != pos) {
		increment(shifter);
		*position = *shifter;
		
		position = shifter;
	}
	
	// add the last element to the position
	*pos = element;
	_size++;
}


/**
 *	Inserts the element while preserving the sorted order of the elements
 *
 *	NOTE: this should only be called on the I'th buffer
 */
void Buffer::insertSortedWithMerging (monom_t element, bool &merged) {
	assert(_bufferType == I_BUF_T);			// sanity checks
	assert(!isFull());		
	
	merged = false;
	
	if (_size == 0) {
		push_back(element);
		return;
	}
	
	// cout << "Degree = " << GET_DEGREE(element) << " --> "; print();
	
	mem_t pos = _front;
	
	buf_size_t counter = 0;
	
	monom_t monomial = (monom_t) element;
	// get attributes
	deg_t elemDegree = GET_DEGREE(monomial);
	coef_t elemCoef = GET_COEF(monomial);
	
	// loop over all elements to find correct place for new element
	while (counter < _size) {
		// if the degree of the element exists, merge the two
		deg_t tempDegree = GET_DEGREE(*pos);
		// cout << "elemDegree = " << elemDegree << "\ttempDegree = " << tempDegree << endl;
		
		if ((deg_t)elemDegree == (deg_t)tempDegree) {
			// cout << "YES" << endl;
			setCoef(*pos, GET_COEF(*pos) + elemCoef);	// add coefs
			merged = true;
			return;
		}
		// if the degree of the element is greater than <pos>, move up the list
		else if (elemDegree < GET_DEGREE(*pos)) {
			decrement(pos);
			counter++;
		}
		// if the degree of the element is greater than <pos>, break;
		else {
			break;
		}
	}
	
	// if we reached the end of the list, append to the end
	if (counter == _size) {
		push_back(element);
		return;
	}
	
	// otherwise, insert element at position and shift elements downward
	decrement(_back);		// make room for new element
	mem_t shifter = _back;		// used to shift elements downward
	mem_t position = _back;		// used to indicate where the element should be shifted
	
	while (shifter != pos) {
		increment(shifter);
		*position = *shifter;
		
		position = shifter;
	}
	
	// add the last element to the position
	*pos = element;
	_size++;
	
}


/******* Managing Children and Parent *******/


/**
 *	Returns a pointer to the left child
 */
Buffer* Buffer::getLeftChild() {
	return _leftChild;
}

/**
 *	Returns a pointer to the right child
 */
Buffer* Buffer::getRightChild() {
	return _rightChild;
}

/**
 *	Returns a pointer to the parent
 */
Buffer* Buffer::getParent() {
	return _parent;
}

/**
 *	Sets the left child of the buffer
 *
 *	@param left the left child
 */
void Buffer::setLeftChild(Buffer *left) {
	_leftChild = left;
}

/**
 *	Sets the right child of the buffer
 *
 *	@param left the right child
 */
void Buffer::setRightChild(Buffer *right) {
	_rightChild = right;
}

/**
 *	Sets the parent of the buffer
 *
 *	@param par the parent
 */
void Buffer::setParent(Buffer *par) {
	_parent = par;
}


/*
 *  fill: fills the current buffer by recursively filling its children until that buffer is full
 *
 *  @param buffer: a pointer to the buffer that will be filled
 */
void Buffer::fill () {
    // EFFECTS: if after filling, a child of 'buffer' is empty, then that child
    //          is set to be exhausted
    
	// if this buffer is exhausted or if its a leaf, then it can't be filled; return
    if (isExhausted() || isLeaf()) {
        return;
    }
		
	assert(_leftChild != NULL);		// sanity check. Note right child CAN be NULL (if we're the A-buffer of the last link)
	

    // if both children are exhausted AND current buffer is empty, set it as exhausted
    if (_leftChild->isExhausted()) {
		if (_rightChild == NULL) {
			if (isEmpty()) 
				_exhausted = true;
			
			return;
		}
		else if (_rightChild->isExhausted()) {
			if (isEmpty()) 
				_exhausted = true;
			
			return;
		}
    }

    
	// fill the buffer
    while (!isFull()) {
        // if the _leftChild is not exhausted and is empty, fill it (if it is not a leaf)
        if (!_leftChild->isExhausted() && _leftChild->isEmpty()) {
            if (_leftChild->isLeaf()) {
                _leftChild->setExhausted(true);  // set leaf to be exhausted since it is empty
            }
            else
                _leftChild->fill();
        }
        
        // if the _rightChild is not exhausted and is empty, fill it (if it is not a leaf)
        if (_rightChild != NULL ) {
			if (!_rightChild->isExhausted() && _rightChild->isEmpty()) {
            	if (_rightChild->isLeaf()) {
                	_rightChild->setExhausted(true);
            	}
            	else
                	_rightChild->fill();
        	}
		}
        
        if (_leftChild->isExhausted()){ 
			if (_rightChild == NULL)
				break;
			else if (_rightChild->isExhausted()) {
            	break;
        	}
		}	

        
        mergeStep();
    }

    
    // assertions
        // // 1) if buffer is empty and its children are exhausted, then buffer is also exhausted
        // assert(!(isEmpty() && _leftChild->isExhausted() && _rightChild->isExhausted() && !isExhausted()));
        // 
        // // 2) if a child is a leaf and empty, then it is exhausted
        // assert(!(_leftChild->isLeaf() && _leftChild->isBufferEmpty() && !_leftChild->isExhausted()));
        // assert(!(_rightChild-?isLeaf() && _rightChild->isBufferEmpty() && !_rightChild->isExhausted()));
 
    
}



/*
 *  mergeStep: pops the max from the children of parent and pushes it into the parent
 *
 *  @param parent: a pointer to the parent buffer
 */
void Buffer::mergeStep () {
    // EFFECTS: if a buffer becomes empty as a result of this function, then it is marked
    //          as exahusted
	checkInvariant();
	if (_leftChild == NULL && _rightChild == NULL)
		return;
		
	// if the left child is empty but it is NOT exhausted, fill it
	if(_leftChild->isEmpty() && !_leftChild->isExhausted())
		_leftChild->fill();
		
	// if the left child is empty but it is NOT exhausted, fill it
	if(_rightChild != NULL && _rightChild->isEmpty() && !_rightChild->isExhausted())
		_rightChild->fill();
		
    if (_leftChild->isExhausted()) {
		if (_rightChild == NULL)
			return;	
		else if(_rightChild->isExhausted())
        	return;
    }
	
	
	assert(_leftChild != NULL);		// sanity check. Note: the right child can be NULL if we're the A-buffer of the last link


	// initialize both values to 0
	deg_t leftMaxDegree = 0;
	deg_t rightMaxDegree = 0;
	
	// left-right valid show if the values of rightMax and leftMax are 'valid' (depending on whether or not the buffer was empty)
	bool leftValid = false;
	bool rightValid = false;
	
	// get (peek) max from non empty buffers
	if (_leftChild != NULL && !_leftChild->isExhausted()) {
		leftMaxDegree = GET_DEGREE(_leftChild->front());
		leftValid = true;
	}
	
	
	if (_rightChild != NULL) {
		if (!_rightChild->isExhausted()) {
			rightMaxDegree = GET_DEGREE(_rightChild->front());
			rightValid = true;
		}
	}
    
	// if neither value is valid, don't merge anything
	if (!leftValid && !rightValid) {
		return;
	}
	
	monom_t maxValue = 0;		// will hold the value to be popped from the child holding the max element
	Buffer* chosenBuffer = NULL;		// will point to the child buffer that has the max element (and will be popped)
	
	// if both values are valid, take the max
	if (leftValid && rightValid) {
		if (leftMaxDegree >= rightMaxDegree) {		// if the left value is greater (or equal) to the right value, pop value from the left buffer
			chosenBuffer = _leftChild;
		}
		else {			// if the right value is greater (or equal) to the left value, pop value from the right buffer
			chosenBuffer = _rightChild;
		}
	}
	// if the left value is valid, but the right is not, then pop the value from the left buffer
	else if (leftValid && !rightValid) {
		chosenBuffer = _leftChild;
	}
	// if the right value is valid, but the left is not, then pop the value from the right buffer
	else {
		chosenBuffer = _rightChild;
	}
	
		
	assert (chosenBuffer != NULL);		// sanity check
	
	// remove the max element from the chosen child and give it to the parent buffer
	maxValue = chosenBuffer->front(); 
	chosenBuffer->pop_front();
	push_back(maxValue);
	
	// if we removed an element from a B-buffer, then we need to decrement the number of elements
	// in the lower part of the current link
	if (chosenBuffer->getType() == B_BUF_T) {
    	// decrement number of elements in lower
		int amountToDecrement = 1;
		amountToDecrement++;

		// INSERT CODE HERE
	}
	
	checkInvariant();
}




void Buffer::print() {
	if (_bufferType == I_BUF_T)
		cout << "I ";
	else if (_bufferType == A_BUF_T) 
		cout << "A ";
	else if (_bufferType == B_BUF_T)
		cout << "B ";
	else if (_bufferType == INTERIOR_BUF_T) 
		cout << "INTERIOR ";
	else if (_bufferType == LEAF_BUF_T)
		cout << "LEAF ";
		
	cout << "Ex = " << _exhausted << " Size = " << _size << "/" << _capacity << ": [";
	
	mem_t pos = _front;
	
	const bool limit = true;
	buf_size_t limitValue = 10;
	
	for (buf_size_t i = 0; i < _size; i++) {
		if (i > 0) 
			cout << ", ";
			
		cout << GET_DEGREE(*pos);
		decrement(pos);
		
		if (limit && i == limitValue) {
			cout << " ... ";
			break;
		}
			
	}
	
	cout << "] (" << GET_DEGREE(*_back) << ")" << endl;
	
}


monom_t getLargestElementInChild(Buffer* b) {
	
	if (b == NULL)
		return 0;
		
	if (b->isExhausted()) 
		return 0;
		
	if (!b->isEmpty()) {
		return b->front();
	}
	
	
	// if b is empty but not exhausted, check its children
	monom_t leftElem = getLargestElementInChild(b->getLeftChild());
	monom_t rightElem = getLargestElementInChild(b->getRightChild());
	
	return (GET_DEGREE(leftElem) > GET_DEGREE(rightElem)) ? leftElem : rightElem;
	
}


monom_t getSmallestInParent(Buffer* b) {
	while (b != NULL && b->getType() != I_BUF_T) {
		if (b->isEmpty())
			b = b->getParent();
		else
			return b->back();
	}
	
	return 1000000;
}


void Buffer::checkInvariant() {
	// if the buffer is empty or it is the I-buffer, then the invariant is not broken
	if (isEmpty() || _bufferType == I_BUF_T) 
		return;
		
	// otherwise, check that the elements in the buffer are smaller than those of its parent, but larger
	// than the elements in its children
	deg_t leftLargestDegree = GET_DEGREE(getLargestElementInChild(_leftChild));		// largest degree in left path
	deg_t rightLargestDegree = GET_DEGREE(getLargestElementInChild(_rightChild));	// largest degree in right path
	deg_t parentSmallestDegree = GET_DEGREE(getSmallestInParent(_parent));			// smallest degree in parent path (i.e. upwards)

	assert (parentSmallestDegree >= GET_DEGREE(*_front));
	assert (leftLargestDegree <= GET_DEGREE(*_back));
	if (rightLargestDegree > GET_DEGREE(*_back)) {
		cout << "RLD = " << rightLargestDegree << ", ME " << GET_DEGREE(*_back) << endl;
	}
	assert (rightLargestDegree <= GET_DEGREE(*_back));
	// if (_parent != NULL && _parent->getType() != I_BUF_T && !_parent->isEmpty()) {
	// 	if (GET_DEGREE(*_front) > GET_DEGREE(_parent->back())) {
	// 		_parent->print();
	// 		print();
	// 	}
	// 	assert (GET_DEGREE(*_front) <= GET_DEGREE(_parent->back()));
	// }
	// if (_leftChild != NULL && !_leftChild->isEmpty()) {
	// 	assert (GET_DEGREE(*_back) >= GET_DEGREE(_leftChild->front()));
	// }
	// if (_rightChild != NULL && !_rightChild->isEmpty()) {
	// 	assert (GET_DEGREE(*_back) >= GET_DEGREE(_rightChild->front()));
	// } 
}




