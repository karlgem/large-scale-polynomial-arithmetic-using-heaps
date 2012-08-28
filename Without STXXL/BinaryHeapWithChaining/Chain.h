/*
 *  Chain.h
 *
 *  Created by Karl Gemayel on 16/5/2012
 *  Copyright 2012 American University of Beirut. All rights reserved.
 *
 */


#ifndef Chain_h
#define Chain_h

#include <list>

#include "../../Polynomial/Monomial.h"
#include "../../Options/Options.h"



class Chain {
private:
	bool valid;
	
	std::list<monom_t> elements;
public:
	
	
	Chain () {
		valid = false;
	}
	
	Chain (size_t cap) {
		// setCapacity(cap);		
	}
	
	void setCapacity(size_t cap);
	
	
	unsigned long long getSize();
	
	
	bool isEmpty();
	
	
	void push(monom_t monomial);
	
	monom_t pop();
	
	/*
	 *	Validation of a Chain:
	 *		A chain is said to be valid if and only if a monomial of degree equal
	 *		to that represented by the chain is already in the heap.
	 */
	bool isValid ();
	
	void setValid(bool value);
	
};


#endif