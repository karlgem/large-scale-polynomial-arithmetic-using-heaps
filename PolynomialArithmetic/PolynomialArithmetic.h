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
#include "../Polynomial/Polynomial.h"
#include "../Heap/Heap.h"

using namespace std;


/*
 * Initializes the prime on which the 'modulo' operation will be based on
 *
 * @param prime: the prime number
 */
void init(int prime);


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




#endif