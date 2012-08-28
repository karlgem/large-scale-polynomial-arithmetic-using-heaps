#include "Chain.h"

void Chain::setCapacity(size_t cap) {
	// try {
	// 	elements.reserve(cap+1);
	// }
	// catch(bad_alloc& ba) {
	// 	cout << "Chain.setCapacity (cap = " << elements.capacity() << "): bad_alloc caught: " << ba.what() << endl;
	// }
}
	

unsigned long long Chain::getSize() {
	return elements.size();
}



bool Chain::isEmpty() {
	return elements.size() == 0;
}


void Chain::push(monom_t monomial) {		
	try {
		elements.push_back(monomial);		// add the element to the end of the chain
	}
	catch (bad_alloc& ba) {
		cout << "Chain.push: bad_alloc caught: " << ba.what() << endl;
	}
	
}

monom_t Chain::pop() {
	if (elements.size() == 0) {
		printError("Chain is empty. Can't pop element.");
		return 0;
	}
	
	monom_t popped = elements.back();		// remove the last element in the chain
	elements.pop_back();
}


/*
 *	Validation of a Chain:
 *		A chain is said to be valid if and only if a monomial of degree equal
 *		to that represented by the chain is already in the heap.
 */
bool Chain::isValid () {
	return valid;
}

void Chain::setValid(bool value) {
	valid = value;
}
