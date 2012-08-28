/*
 *  Polynomial.h
 *
 *  Created by Karl Gemayel on 16/5/2012
 *  Copyright 2012 American University of Beirut. All rights reserved.
 *
 */

/**
 * \file Polynomial.h
 * \brief A polynomial is defined to be a (degree-sorted) sequence of monomials.
 *
 *	This file describes the polynomial and how it is implemented.
 */

#ifndef Polynomial_h
#define Polynomial_h

// include header files
#include <vector>
#include "Monomial.h"

using namespace std;


// define the polynomial type
typedef vector<monom_t> poly_t; 	/**< \typedef definition of polynomial type */


#endif