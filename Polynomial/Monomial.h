/*
 *  Monomial.h
 *
 *  Created by Karl Gemayel on 16/5/2012
 *  Copyright 2012 American University of Beirut. All rights reserved.
 *
 */

/**
 * \file Monomial.h
 * \brief A monomial is defined by its degree and coefficient.
 *
 * This file describes the monomial and how it is implemented. It also provides
 * some methods that present read/write/set functionalities on monomials.
 *
 *		MONOMIAL REPRESENTATION AND MANIPULATION:
 *
 * 		Below are the masks, bits-per-section, and the word type.
 * 		How to use: \newline
 *
 *		1) The first line defines the word type (monom_t). Note that if 
 *			the type is changed, some modifications should be made to 
 *			the masks and the number of bits per section.
 *
 *		2) The next three lines define the bit-masks:
 *			a) DEGREE_MASK: represents the degree section in a word.
 *			b) COEF_MASK:	represents the coefficient section in a word.
 *			c) ID_MASK:     represents the id section in a word.
 *
 *		3) The last two lines define the amount to shift for each variable:
 *			a) COEF_SHIFT:  displacement of coefficient from the right.
 *			b) ID_SHIFT:	displacement of the id from the right.
 *
 *		CAUTION: The masks, bits-per-section, and word type are all dependent. Therefore, 
 *			 a change to one of them requires a change to all.
 *			 Here are some guidelines:
 *				1) First define monom_t to get maximum number of bits per word
 *				2) Concatenating DEGREE_MASK, COEF_MASK, and ID MASK should result 
 *				   in a word of length = (sizeof(monom_t) * 8), and its binary digits 
 *				   are all 1's 
 *
 *		NOTE: a deviation from any of these rules might result in undefined behavior.
 *
 */

#ifndef Monomial_h
#define Monomial_h

// FOR NOW
typedef unsigned long long w_type;

// define the monomial, degree, coefficient, and ID types

typedef unsigned long long monom_t;		/**< \typedef definition of monomial type */
typedef unsigned long long deg_t;		/**< \typedef definition of degree type */
typedef unsigned long long coef_t;		/**< \typedef definition of coefficient type */
typedef unsigned long long ID_t;		/**< \typedef definition of ID type */


#define DEGREE_MASK 0x00000000FFFFFFFFULL	/**< mask representing degree */
#define COEF_MASK   0x00000FFF00000000ULL	/**< mask representing coef	*/
#define ID_MASK     0xFFFFF00000000000ULL   /**< mask representing id */    

#define COEF_SHIFT  32          /**< displacement of coef from right */
#define ID_SHIFT    44          /**< displacement of id from right */

// the below three macros automatically unpack id, coef, or degree
// based on the 'masks' and 'shifts' above
#define GET_ID(W)      ((ID_t) ((W & ID_MASK) >> ID_SHIFT))			/**< gets the ID of the monomial */
#define GET_COEF(W)    ((coef_t) ((W & COEF_MASK) >> COEF_SHIFT))		/**< gets the coefficient of the monomial */
#define GET_DEGREE(W)  ((deg_t) ((W & DEGREE_MASK)))					/**< gets the degreeof the monomial */



/**
 * Creates a monomial with a degree, coefficient, and ID.
 *
 * @param degree the degree 
 * @param coef the coefficient
 * @param ID the ID
 * @return Returns a monomial containing a degree, coefficient and an ID
 */
inline monom_t createMonomial(deg_t degree, coef_t coef, ID_t ID) {
    monom_t word = degree;
    
    word = word | (coef << COEF_SHIFT);
    word = word | (ID << ID_SHIFT);
    
    return word;
}


/**
 * Sets the degree of a monomial.
 *
 * \param monomial the monomial which we want to set the degree of
 * \param degree the value of the new degree to be set
 */
inline void setDegree(monom_t &monomial, monom_t degree) {
    monomial = monomial & (ID_MASK | COEF_MASK);
    monomial = monomial | (degree);
}

/**
 * Sets the coefficient of a monomial.
 *
 * \param monomial the monomial which we want to set the coefficient of
 * \param coef the value of the new coefficient to be set
 */
inline void setCoef(monom_t &monomial, monom_t coef) {
    monomial = monomial & (ID_MASK | DEGREE_MASK);
    monomial = monomial | (coef << COEF_SHIFT);
}

#endif