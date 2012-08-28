/*
 *	\class Buffer.h
 *	\brief A buffer that holds elements 
 *
 *	Every buffer has a fixed capacity, a parent buffer, and a left child and right child (which can be NULL)
 *
 *	Implementation Note: The buffer is implemented as a circular buffer in order to optimize
 *	performance. Moreover, given the context that this buffer will be used in and the requirement
 *	that all elements and overhead must be stored contiguously in memory, the buffer is passed
 *	(at initialization) a memory location which represents the start of where this buffer is located
 *	in memory. All subsequent elements of this buffer must be stored contiguously starting at the
 *	specified location.
 */

#ifndef BUFFER_H
#define BUFFER_H

#include "../Polynomial/Monomial.h"
#include "../Options/Options.h"

// define the types of buffers
#define I_BUF_T 0 
#define A_BUF_T 1
#define B_BUF_T 2
#define INTERIOR_BUF_T 3
#define LEAF_BUF_T 4

// typedef representing buffer's size
typedef unsigned long long buf_size_t;



class Buffer {
private:
	
	mem_t _loc;	/**< the starting memory location for the whole buffer (including overhead)*/
	mem_t _elemLoc; /**< starting memory location where elements can be inserted */
	
	// pointers that indicate the list of elements in the buffer
	mem_t _front;		/**< points to the head of the elements */
	mem_t _back;		/**< points to the tail of the elements */
	
	// buffer characteristics
	buf_size_t _capacity;
	buf_size_t _size;
	bool _exhausted;
	int _bufferType;
	
	// parent and children
	Buffer* _parent;
	Buffer* _leftChild;
	Buffer* _rightChild;
	
	/*
	 *	Increments the pointer in a circular manner
	 */
	void increment (mem_t& p);
	
	/*
	 *	Decrements the pointer in a circular manner
	 */
	void decrement (mem_t& p);
	
	
public:	
	/**
	 *	Constructor: initializes the buffer (at the given memory location) with the
	 *	specified capacity.
	 *
	 *	@param mem the starting memory location at which to allocate the buffer's components
	 *	@param cap the buffer's maximum capacity
	 *	@param bType indicates the type of the buffer (I, Interior, Leaf)
	 */
	Buffer(mem_t& mem, buf_size_t cap, int bType);
	
	
	~Buffer() {
		_parent = NULL;
		_leftChild = NULL;
		_rightChild = NULL;
	}
	
	/**
	 *	Returns the buffer's capacity
	 *	@return capacity
	 */
	buf_size_t getCapacity();
	
	/**
	 *	Returns the buffer's size
	 *	@return size
	 */
	buf_size_t getSize();
	
	/**
	 *	Checks if the buffer is empty
	 *
	 *	@return returns true if the buffer is empty; false otherwise
	 */
	bool isEmpty();
	
	/**
	 *	Checks if the buffer is full
	 *
	 *	@return returns true if the buffer is full; false otherwise
	 */
	bool isFull();
	
	/**
	 *	Checks if the buffer is exhausted (i.e. its children and itself are empty)
	 *
	 *	@return returns true if the buffer is exhausted; false otherwise
	 */
	bool isExhausted();
	
	/**
	 *	Returns true if the buffer is a Leaf
	 */
	bool isLeaf() {
		return _bufferType == LEAF_BUF_T;
	}
	
	/**
	 *	Sets the exhausted value of the buffer
	 */
	void setExhausted(bool value) {
		_exhausted = value;
	}
	/**
	 *	Returns the type of the buffer
	 *
	 *	@return type
	 */
	int getBufferType();
	
	/**
	 *	Inserts an element to the end of the buffer
	 *	
	 *	@param element inserted element
	 */
	void push_back(monom_t element);
	
	/**
	 *	Extracts the first element of the buffer
	 *
	 *	@return the first element
	 */
	monom_t pop_front();
	
	/**
	 *	Peeks at the front element without actually extracting it
	 *
	 *	@return a copy of the front element
	 */
	monom_t front();
	
	
	/**
	 *	Peeks at the back element without actually extracting it
	 *
	 *	@return a copy of the back element
	 */
	monom_t back();
	
	/**
	 *	Inserts the element while preserving the sorted order of the elements
	 */
	void insertSorted (monom_t element);
	
	/**
	 *	Inserts the element while preserving the sorted order of the elements. If
	 *	an element of the same degree already exists in the heap, then the two elements
	 *	are merged together, and the boolean parameter 'merged' is set to true
	 *
	 *	@param element the element to be inserted
	 *	@param merged set to true if the element has been merged with another of the same degree
	 */
	void insertSortedWithMerging (monom_t element, bool &merged);
	
	/******* Managing Children and Parent *******/
	
	/**
	 *	Returns a pointer to the left child
	 */
	Buffer* getLeftChild();
	
	/**
	 *	Returns a pointer to the right child
	 */
	Buffer* getRightChild();
	
	/**
	 *	Returns a pointer to the parent
	 */
	Buffer* getParent();
	
	/**
	 *	Sets the left child of the buffer
	 *
	 *	@param left the left child
	 */
	void setLeftChild(Buffer* left);
	
	/**
	 *	Sets the right child of the buffer
	 *
	 *	@param left the right child
	 */
	void setRightChild(Buffer* right);
	
	/**
	 *	Sets the parent of the buffer
	 *
	 *	@param parent the parent
	 */
	void setParent(Buffer* parent);
		
	/*
	 *  fill: fills the current buffer by recursively filling its children until that buffer is full
	 *
	 *  @param buffer: a pointer to the buffer that will be filled
	 */
	void fill ();
	
	/*
	 *  mergeStep: pops the max from the children of parent and pushes it into the parent
	 *
	 *  @param parent: a pointer to the parent buffer
	 */
	void mergeStep ();
	
	int getType() {
		return _bufferType;
	}
	
	/**
	 *	Print the elements in the buffer
	 */
	void print();
};

#endif