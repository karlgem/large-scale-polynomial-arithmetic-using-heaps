// PolynomialGenerator.h

/* 
 *	This file defines methods to generate polynomials for use in the summation of products routine
 */

#ifndef PolynomialGenerator_h
#define PolynomialGenerator_h

#include <string>

using namespace std;
	
/*
 *	Generates a file containing two sets of polynomials. The generated file will have
 *	the following format:
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
 *	@param n the number of polynomial pairs
 *	@param s the sparsity of the polynomials
 *	@param filename the name of the generated file
 */
void generatePolynomials (unsigned int n, double s, const string filename);

#endif
