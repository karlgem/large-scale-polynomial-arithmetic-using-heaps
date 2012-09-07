/*
 *  Chain.h
 *
 *  Created by Karl Gemayel on 16/5/2012
 *  Copyright 2012 American University of Beirut. All rights reserved.
 *
 */

/**
 * \file Chain.h
 * \brief A chain contains a set of elements (usually having a common charactersitic)
 *
 * This file describes the chain and how it is implemented. It also provides
 * some methods that present read/write/set functionalities on chains.
 *
 *		CHAIN REPRESENTATION AND MANIPULATION:
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

#ifndef Chain_h
#define Chain_h


typedef heap_t::iterator chain_t;		/**< \typedef definition of chain type */
typedef unsigned long long chain_size_t;	/**< \typedef definition of chain size type */



/**
 * THE BELOW METHODS DEAL WITH CHAINING
 *
 * Implementation Details:
 * 	The container that holds the binary heap is expanded to hold the chains. The global 
 * 	pointer 'chains' points to the start of the first chain. The set of chains is layed
 *	out consecutively in memory. The first element of each chain contains the size and
 * 	capacity of that chain.
 */

/** Reading/Writing the Size and Capacity of a Chain */
#define CHAIN_CAP_MASK  0xFFFFFFFF00000000ULL	/**< mask representing chain size */
#define CHAIN_SIZE_MASK 0x00000000FFFFFFFFULL	/**< mask representing chain capacity */

#define CHAIN_CAP_SHIFT  32          /**< displacement of capacity from right */

// the below two macros automatically unpack the size and capacity of a chain
// based on the 'masks' and 'shifts' above
#define GET_CHAIN_CAP(W) ((chain_size_t) ((W & CHAIN_CAP_MASK) >> CHAIN_CAP_SHIFT))	/**< gets the capacity of the chain */
#define GET_CHAIN_SIZE(W)  ((chain_size_t) ((W & CHAIN_SIZE_MASK)))					/**< gets the size of the chain */


#define CHAIN_OVERHEAD 1
#define CHAIN_PARAMETERS_INDEX 0

/**
 * Creates a word representing chain parameters: size and capacity
 *
 * @param chainSize the chain's size
 * @param chainCapacity the chain's capacity
 * @return Returns a word containing the size and capacity of the chain
 */
inline chain_size_t createSizeCapPair(chain_size_t chainSize, chain_size_t chainCapacity) {
    chain_size_t word = chainSize;
    word = word | (chainCapacity << CHAIN_CAP_SHIFT);

    return word;
}

inline void setChainSize (chain_t &chain, chain_size_t chainSize) {
	chain_size_t chainParameters = (*(chain + CHAIN_PARAMETERS_INDEX));
	chainParameters = chainParameters & CHAIN_CAP_MASK;		// reset the chain size to 0
	chainParameters = chainParameters | (chainSize & CHAIN_SIZE_MASK);		// set the new chain size
	*(chain + CHAIN_PARAMETERS_INDEX) = chainParameters;
}

inline chain_size_t getChainSize(chain_t chain) {
	return GET_CHAIN_SIZE(chain[CHAIN_PARAMETERS_INDEX]);
}

inline void incrementChainSize (chain_t &chain) {
	chain_size_t newChainSize = getChainSize(chain) + 1;		// increment chain size by 1
	setChainSize(chain, newChainSize);
}

inline void decrementChainSize (chain_t &chain) {
	chain_size_t newChainSize = getChainSize(chain) -1;		// increment chain size by 1
	setChainSize(chain, newChainSize);
}




inline chain_size_t getChainCapacity(chain_t chain) {
	return GET_CHAIN_CAP(chain[CHAIN_PARAMETERS_INDEX]);
}

inline bool chainIsEmpty(chain_t chain) {
	return (getChainSize(chain) == 0 ? true : false);
}


inline void setChainParameters (chain_t &chain, chain_size_t chainSize, chain_size_t chainCapacity) {
	chain[CHAIN_PARAMETERS_INDEX] = createSizeCapPair(chainSize, chainCapacity);
}

#endif