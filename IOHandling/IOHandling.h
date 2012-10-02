/*
 *  IOHandling.h
 *
 *  Created by Karl Gemayel on 22/09/2012
 *  Copyright 2012 American University of Beirut. All rights reserved.
 *
 */

/* 
 *	This file defines defines methods that handle I/O operations such as
 *	reading input from files.
 */

#ifndef IOHandling_h
#define IOHandling_h

#include <vector>
#include <string>

#include "../Polynomial/Polynomial.h"

using namespace std;


/*
 *	Reads the two sets of polynomials from the input file. The file must be formatted
 *	according to the following guidelines:
 *		- The first line is a number (call it N) indicating the number of polynomial 
 *			pairs in the file.
 *		- The next N lines are polynomials constituting the first set.
 *		- The final N lines are polynomials constituting the second set.
 *
 *	NOTE: The format of the lines containing polynomials is:
 *				p = a1*x^n1 + a2*x^n2 + ... 
 *			where:
 *				- p can be anything WITHOUT the "=" sign (eg f_1, g_2, ...)
 *				- a1, a2, ... are positive integers
 *				- n1, n2, ... are positive integers (in decreasing order)
 *
 *	@param inputName the name of the input file
 *	@param set1 vector that will be filled with the first set of polynomials
 *	@param set2 vector that will be filled with the second set of polynomials
 *	@param numberOfPairs will be set to be the number of pairs
 */
void readInputFile (string inputName, vector<poly_t> &set1, vector<poly_t> &set2, int &numberOfPairs);


#endif