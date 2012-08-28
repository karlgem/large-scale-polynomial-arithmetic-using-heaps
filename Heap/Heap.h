/*
 *  Heap.h
 *  Generic Heap Interface
 *
 *  Created by Karl Gemayel on 16/05/2012.
 *  Copyright 2012 American University of Beirut. All rights reserved.
 *
 */

/** 
 * \class Heap
 * \brief An interface defining the required functionality of a heap
 *
 * The heap will need to perform the job of a max heap; i.e. be able to return the 
 * maximum element when required.
 */

#ifndef Heap_h
#define Heap_h

#include "../Polynomial/Monomial.h"
#include "../Options/Options.h"
#include <vector>
#include <stxxl.h>

using namespace std;

//typedef std::vector<w_type> heap_t;

typedef unsigned long long heap_size_t;

class Heap {

public:
	// functions
	Heap () {};
	
	virtual ~Heap() {}
	/**	
	 * Inserts an element that has a degree, coefficient, and an id into the heap
	 * in its correct position.
	 *
	 * \param degree the degree of the element
	 * \param coef the coefficient of the element
	 * \param the id of the element
	 */
	virtual void insert (deg_t degree, coef_t coef, ID_t f_id) {}
	virtual void insert (deg_t degree, coef_t coef) {}

	/** 
	 * Returns the value of the max element of the heap without actually extracting it.
	 *
	 * @return The value of the max element 
	 */
	virtual monom_t peek() { return 0; }
	
	/** 
	 * Extracts the max element of the heap and returns it.
	 *
	 * \return The value of the max element
	 */
	virtual monom_t poll() { return 0; }
	
	/**
	 * Returns the size of the heap.
	 *
	 * \return The size of the heap
	 */
	virtual heap_size_t size() { printf("size HEAP\n"); return 2; }
	
	/** 
	 * Returns true if the heap is empty, and false otherwise
	 *
	 * \return if heap is empty - true. Otherwise - false
	 */
	virtual bool isEmpty() { printf("isEmpty HEAP\n"); return false; };
	
	
	virtual void print() {};
};

#endif

