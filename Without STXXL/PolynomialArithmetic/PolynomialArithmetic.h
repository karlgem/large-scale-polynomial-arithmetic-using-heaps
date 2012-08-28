//
//  Polynomial_Operations.h
//  Polynomial Operations
//
//  Created by Karl Gemayel on 3/28/12.
//  Copyright (c) 2011 American University of Beirut. All rights reserved.
//


#ifndef PolynomialArithmetic_h
#define PolynomialArithmetic_h

#include <vector>
#include <string>
#include "../../Polynomial/Polynomial.h"

using namespace std;

size_t polynomialSize(const poly_t &p);

/*
 * Initializes the prime on which the 'modulo' operation will be based on
 *
 * @param prime: the prime number
 */
void init(int prime);


// NOTE: All packed representations of defined in Monomial_Packing.h


/*
 * Multiplies multiple polynomial pairs (f_i * g_i), and sums up their 
 * results.
 *
 * @param f: the first set of polynomials
 * @param g: the second set of polynomials
 * @param result: an allocated poly_t variable where the result will be stored
 *
 */
void summationOfProducts (std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result);


/*
 * Multiplies a single pair of polynomials, f and g.
 *
 * @param f: the first polynomial
 * @param g: the second polynomial
 * @param result: an allocated poly_t variable where the result will be stored
 *
 * @return Returns the product of f and g
 */
void multiplySinglePair(const poly_t &f, const poly_t &g, poly_t &result);


/*
 * Multiplies multiple polynomial pairs (f_i * g_i), and sums up their 
 * results.
 *
 * @param f: the first set of polynomials
 * @param g: the second set of polynomials
 * @param result: an allocated poly_t variable where the result will be stored
 *
 */
void multiplyMultiplePairs(std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result);
void multiplyMultiplePairsFunnelWithMerging(std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result);


/*
 * Subtracts the polynomials f and g. The polynomials f and g both have a
 * packed representation.
 *
 * @param f: the first polynomial
 * @param g: the second polynomial
 * @param result: an allocated poly_t variable where the result will be stored
 *
 */
void subtract(const poly_t &f, const poly_t &g, poly_t &result);


/*
 * Returns the degree of the polynomial
 *
 * @param P: the polynomial.
 */
deg_t deg(const poly_t &p);


/*
 * Parses the polynomial string into a poly_t type. The polynomial string 
 * has the same form as:    a*x^n + b*x^(n-1) + ... + c*x^0  (where a != 0)
 * The resulting polynomial will have size n.
 *
 *
 * @param p: the polynomial to be parsed
 * @param result: a pointer to the output generated from this method (in
 *                  packed representation)
 */
void parsePolynomialString (string p, poly_t &result);


/*
 *  Converts the polynomial representation to a string of the form:
 *      a*x^n + b*x^(n-1) + ... + c*x^0  (where a != 0)
 *
 *  @param p: the polynomial 
 *  @param output: string representation of the polynomial p
 */
void toString(const poly_t &p, string &output);

/*
 * Returns the number of monomials in a given string polynomial
 *
 * @param f: the polynomial
 */
int getNumberOfMonomials (const char *f);

void printPolynomial(const poly_t &p);



#endif