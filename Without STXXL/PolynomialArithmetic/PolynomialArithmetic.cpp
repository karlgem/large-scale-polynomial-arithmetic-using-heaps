#include <iostream>
#include <sstream>
#include <string>
#include <assert.h>
#include <vector>
#include "PolynomialArithmetic.h"
#include "../../Options/Options.h"

// heap
#include "../BinaryHeapWithChaining/BinaryHeapWithChaining.h"

#define ASSERTIONS false

#define ENABLE_STATE_PRINTING false
#define TURN_ON_PAUSES false
#define ASSERT_REL_BOUNDS false

using namespace std;

int p = 113;                      // prime 
void init(int prime) {
    p = prime;
}

// function prototypes
inline void skipSpaces (char &current, size_t &index, string line);

// OSI stands for optimized sequence of inserts
void OSImultiply(std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result);



// returns the degree of the polynomial
deg_t deg(const poly_t &p) {
	return GET_DEGREE(p[0]);
}












void printPolynomial(const poly_t &p) {
    for (int i = 0; i < p.size(); i++) {
        
        if (i != 0) {
            printf(" + ");
        }
        

        if (GET_COEF(p[i]) != 0) {
            printf("%llu", GET_COEF(p[i]));
        }
        
        if (GET_DEGREE(p[i]) > 0) {
            printf("*x");
        }
        
        if (GET_DEGREE(p[i]) > 1) {
            printf("^%llu", GET_DEGREE(p[i]));
        }
    }
    
    printf("\n");
}


void toString(const poly_t &p, string &output) {
		
    stringstream ss (stringstream::in | stringstream::out);
        
    output = "";
    
    for (int i = 0; i < p.size(); i++) {
        
        if (i != 0) {
            ss << " + ";
        }
        
        
        if (GET_COEF(p[i]) != 0) {
            ss << GET_COEF(p[i]);
        }
        
        if (GET_DEGREE(p[i]) > 0) {
            ss << "*x";
        }
        
        if (GET_DEGREE(p[i]) > 1) {
            ss << "^"  << GET_DEGREE(p[i]);
        }
    }
    
    output = ss.str();
}


// multiply multiple pairs
void multiplyMultiplePairsFunnelWithMerging(std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result) {
	/*
	 * General Code Sketch:
	 *
	 *  1)  Let S = n_1 + n_2 + ... + n_k, where n_i is the number of 
     *      monomials in f_i (for 1 <= i <= k).
     *  2)  Insert the S maximum products of monomials (from all (f_i, g_i) 
     *      pairs) into the heap.
     *  3)  Repeat until all pairs have been fully traversed:
     *          i)  Extract the maximum elemement from the heap
     *              and place it in the result
     *         ii)  Insert the next maximum element into the heap
     *  4)  When all pairs have been traversed, extract all elements
     *      from the max heap and place them in the result
	 *
	 * Note on managing ID's:
	 *	1)	An array of stacks will be used to keep record of what elements
	 *	   	inserted into the heap correspond to which monomial multiplication:	
	 *			a) Each row indicates a certain degree (i.e. 0'th row = degree 0)
	 *			b) Each column in a row contains the id of an element that was 
	 *			   inserted.
	 *	2)	Another array will hold what the next multiplication corresponding
	 *		to a certain ID should be.
	 *	
	 *
	 */
    
	if (verboseLevel(VERBOSE_LOW)) 
		cout << "\\* Starting Multiplication Procedure \"Non-Optimized\"! (" << g_polynomials.size() << " pairs) */" << endl;


    size_t numOfPairs = g_polynomials.size();
    
	// get maximum number of multiplications that will be done 
	size_t maxMultiplications = 0;
	for (int i = 0; i < numOfPairs; i++) {
		maxMultiplications += f_polynomials[i].size() * g_polynomials[i].size();
	}
    
    size_t currentMultiplications = 0;		// used to track the number of multiplications already done
	// set n_f as the number of monomials in f
	

    // get the max number of monomials in the f's, and the total heap capacity
	size_t heap_cap = 0;	// capacity of heap
	size_t f_max = 0;
	for (int i = 0; i < numOfPairs; i++) {
		heap_cap += f_polynomials[i].size();
		if (f_polynomials[i].size() > f_max) {
			f_max =(size_t) f_polynomials[i].size();
		}
	}
	
	if (verboseLevel(VERBOSE_MED)) cout << "Retrieving the largest degree of a monomial multiplication..." << endl;
	
	// get the largest degree that can be produced by a mutliplication of f and g
	deg_t largestDegree = 0;
	deg_t possibleLargestDegree = 0;
	for (int i = 0; i < numOfPairs; i++) {
		possibleLargestDegree = GET_DEGREE((f_polynomials[i])[0]) + GET_DEGREE((g_polynomials[i])[0]);
		if (largestDegree < possibleLargestDegree) {
			largestDegree = possibleLargestDegree;
		}
	}

	if (verboseLevel(VERBOSE_MED)) cout << "Largest degree retrieved: " << largestDegree <<  endl << endl;
	
	size_t numberOfChains = largestDegree+1;
		
	if (verboseLevel(VERBOSE_LOW)) cout << "Using Binary Heap (with chaining)" << endl;
	// allocate heap
	BinaryHeapWithChaining *A = new BinaryHeapWithChaining (heap_cap, numberOfChains);
	
	
    
    monom_t f_element, g_element;    // represent an element of f and g
	

	

	if (verboseLevel(VERBOSE_MED)) cout << "Allocating id-monomial product mapping array.." << endl;

	// create an array of vectors that would hold a mapping of a monomial inserted
	// into the heap to the actual f and g monomials whose product is that monomial
//stack<unsigned long> ID_Monom_pair [largestDegree];
	vector<unsigned long> ID_Monom_pair [largestDegree+1];
	

	if (verboseLevel(VERBOSE_MED)) cout << "ID-monomial product mapping array allocated" << endl;



	if (verboseLevel(VERBOSE_LOW)) cout << "\nPerforming the initial insertions..." << endl;
	
	unsigned long long maxNumberOfElementsInHeap = 0;
	// Insert (f1_1*g1_1, f1_2*g1_2,.... f1_n*g1_n).....(fn_1*gn_1....fn_n*gn_n) into the heap
	int id = 0;
	for (int i = 0; i < numOfPairs; i++) {			// loop on all polynomials
		
		poly_t f = f_polynomials[i];
		poly_t g = g_polynomials[i];
		
		for (int j = 0; j < f_polynomials[i].size(); j++) {		// loop on monomials of i'th polynomial
            f_element = f[j];
			g_element = g[0];
            
            coef_t coef = GET_COEF(f_element) * GET_COEF(g_element);
            deg_t degree = GET_DEGREE(f_element) + GET_DEGREE(g_element);
			
			// if (verboseLevel(VERBOSE_HIGH)) 
			// 	cout << "SIZE=" << A->size() << "\tInserting: coef=" << 1 << ", deg=" << degree << endl;
				
			A->insert(degree, 1);
			
			
			if (ENABLE_STATE_PRINTING) {
				cout << "Inserted " << degree << endl << endl;
				// A->printQueueState();
				// 			pausePrompt();
			}
			
			// record the id-monomial relationship
//			stack<unsigned long> rel = ID_Monom_pair[degree];
			if (ASSERT_REL_BOUNDS) {cout << degree << "/" << largestDegree << endl; assert(degree < largestDegree);}
			vector<unsigned long> *rel = &ID_Monom_pair[degree];
			
			//rel.push(id);
			rel->push_back(id);
						
			maxNumberOfElementsInHeap++;
			id++;
			
			currentMultiplications++;
			
		}
		
		if (f_polynomials[i].size() < f_max) {			// go to next level in record table
			id += f_max-f_polynomials[i].size();
		}
	}


	if (verboseLevel(VERBOSE_LOW)) cout << "Done with initial insertions" << endl << endl;


	// print the max number of elements in the heap at any one time
	if (verboseLevel(VERBOSE_LOW)) 
		cout << "The max number of elements in the heap is: " << maxNumberOfElementsInHeap << endl;
	
	/*
	 * The following creates an array of linked lists to reference all monomials of f_i. This is
	 * changed from the previously 2D matrix to save on memory usage. Note that this is important
	 * when dealing with sparse polynomials.
	 */
	
	struct ID {
		int g_counter;
		ID *next;
	};
	
	int n_f = (int) numOfPairs;
	ID *recArray [n_f];
	for (int i = 0; i < n_f; i++) {
		if (f_polynomials[i].size() > 0) {		// check that it has at least one monomial
			recArray[i] = new ID;			
		}
		else {
			continue;
		}
		
		
		ID *current = recArray[i];
		current->g_counter = 0;
		current->next = NULL;
		
		for (int j = 1; j < f_polynomials[i].size(); j++) {		// create link for every monomials of f_i
			current->next = new ID;
			current = current->next;
			current->g_counter = 0;			// added new
			current->next = NULL;
		}
	}
		
	if (verboseLevel(VERBOSE_LOW)) cout << "Performing second set of insertions..." << endl;
	
	
	while (currentMultiplications < maxMultiplications) {	
		monom_t max = A->poll();
		
		
        
		if (max == 0) {
			break;
		}
		
		coef_t maxCoef = GET_COEF(max);
		deg_t maxDegree = GET_DEGREE(max);
		
		if (verboseLevel(VERBOSE_HIGH)) 
			cout << "SIZE=" << A->size() << "\tExtracted: coef=" << maxCoef << ", deg=" << maxDegree << endl;
		
	
		// using the table defining the relationship between ID's and monomials
		// inserted into the heap, check the degree of this max monomial to obtain
		// from what monomial products it resulted.		
		
		// first get the number of ids that correspond to this max
//		stack<unsigned long> relation = ID_Monom_pair[maxDegree];
		if (ASSERT_REL_BOUNDS) {
			assert(maxDegree < largestDegree);
		}
			
		vector<unsigned long> *relation = &ID_Monom_pair[maxDegree];
		size_t numOfIds = relation->size();
				
		// second, for each of the id's, insert the next monomial multiplication
		// that should occur by consulting the table defining the relation between
		// ID's and monomial multiplications
		for (size_t id_counter = 0; id_counter < numOfIds; id_counter++) {
			// get one of the ids corresponding to this relation
			// unsigned long currentMaxID = relation.top();
			// relation.pop();
			
			if (TURN_ON_PAUSES) {
				cout << "1" << endl;
				pausePrompt();	
			}
			
			unsigned long currentMaxID = relation->back();
			relation->pop_back();

			// increment recTable of f(id)
			int maxRow = (int) (currentMaxID / f_max);
			int maxCol = (int) (currentMaxID % f_max);

			if (TURN_ON_PAUSES) {
				cout << "2" << endl;
				pausePrompt();
			}
			
			// find current reference and update g_counter
			ID *maxRef = recArray[maxRow];
			for (int i = 0; i < maxCol; i++) {
				maxRef = maxRef->next;
			}

			if (maxRef->g_counter == g_polynomials[maxRow].size() - 1) 
				maxRef->g_counter = -1;
			else
				maxRef->g_counter++;

			if (TURN_ON_PAUSES) { 
				cout << "3" << endl;
				pausePrompt();
			}

	        if (maxRef->g_counter != -1) {
				// add the successor
				poly_t f = f_polynomials[maxRow];
				poly_t g = g_polynomials[maxRow];

				//			int g_successor = recTable[maxRow][maxCol];
				int g_successor = maxRef->g_counter;

	            f_element = f[maxCol];
	            g_element = g[g_successor];

	            coef_t coef = GET_COEF(f_element) * GET_COEF(g_element);
	            deg_t degree = GET_DEGREE(f_element) + GET_DEGREE(g_element);

//				if (coef % p != 0) {
					if (verboseLevel(VERBOSE_HIGH)) 
						cout << "SIZE=" << A->size() << "\tInserting: coef=" << 1 << ", deg=" << degree << endl;
						
					A->insert (degree, 1);
					
					
					
					// record the id-monomial relationship
					// stack<unsigned long> rel = ID_Monom_pair[degree];
					// rel.push(currentMaxID);
					
					if (TURN_ON_PAUSES) {
						cout << "4" << endl;
						pausePrompt();
					}
					
					if (ASSERT_REL_BOUNDS) assert(degree < largestDegree);
					vector<unsigned long> *rel = &ID_Monom_pair[degree];
					rel->push_back(currentMaxID);
					
					if (TURN_ON_PAUSES) {
						cout << "5" << endl;
						pausePrompt();
					}
					
//				}			
				currentMultiplications++;
			}
			
		}

		// if (maxCoef % p == 0) {
		//             continue;
		//         }
        
        	
//        setCoef(max, maxCoef % p);
		setCoef(max, maxCoef);
        result.push_back(max);
	}
		
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Done with second set of insertions" << endl << endl;
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Extracting all remaining elements from heap" << endl;
	// extract all remaining elements from heap
	
	
	while (!A->isEmpty()) {
		monom_t max = A->poll();
        		
		coef_t maxCoef = GET_COEF(max);
		deg_t maxDegree = GET_DEGREE(max);
		
        
//        setCoef(max, maxCoef % p);
		setCoef(max, maxCoef);
        result.push_back(max);
	}
    
    
    delete A;
	
}






// multiply using an optimized sequence of inserts
void OSImultiply(std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result) {
	if (verboseLevel(VERBOSE_LOW)) 
		cout << "Procedure: Multiple Polynomial Pairs Multiplication using Optimized Sequence of Inserts" << endl;
		
	// make sure that f and g have the same number of polynomials
	if (f_polynomials.size() != g_polynomials.size()) {
		printError("The number of polynomials in the sets f and g should be equal!");
		exit(0);
	}
	
	// get the number of polynomial pairs
	size_t numOfPairs = f_polynomials.size();
	
	
	// FOR NOW: for every pair in f and g, let f hold the polynomial with the larger number of monomials
	for (size_t i = 0; i < numOfPairs; i++) {
		if (f_polynomials[i].size() < g_polynomials[i].size()) {
			// swap the two polynomials
			poly_t temp = g_polynomials[i];
			g_polynomials[i] = f_polynomials[i];
			f_polynomials[i] = temp;
		}
	}
	
	size_t largestSizeReached = 0;
	
	
	// get the (starting) heap capacity
	size_t heap_cap = numOfPairs;

	
	// calculate the largest degree
	deg_t largestDegree = 0;
	deg_t tempDegree = 0;
	for (size_t i = 0; i < numOfPairs; i++) {
		tempDegree = GET_DEGREE((f_polynomials[i])[0]) + GET_DEGREE((g_polynomials[i])[0]);
		if (largestDegree < tempDegree)
			largestDegree = tempDegree;
	}
	
	size_t numberOfChains = largestDegree+1;
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Using Binary Heap (with chaining)" << endl;		

	// allocate heap
	BinaryHeapWithChaining *A = new BinaryHeapWithChaining (heap_cap, numberOfChains);				

	
	
	/* IDs and Counters Management
	 *	- Id's are used to link every monomial product to its orignal polynomial pair. This is 
	 *		needed to identify the source of the successor that needs to be inserted
	 *	- Counters are used to preserve the 'state' of every polynomial pair in terms of where
	 *		it had last reached in its multiplication process. Two counters will be assigned
	 *		to each pair (one for the f polynomial, and the other for the g). 
	 *		
	 * How it all works:
	 *	Once a monomial is extracted, identify the polynomial pair from which it was calculated
	 *	(using ID's). Then, associate the polynomial pair with its respective 'counters'. Proceed
	 *	to perform the calculation (according to the optimized sequence of inserts algorithm) based
	 *	on the state defined by these two counters.
	 */
	
	vector<size_t> degreeToSource [largestDegree+1];		// associates a degree of a monomial to a set of polynomial pairs from which it was generated
	
	// setup the counters
	const unsigned int FCOUNTER_INDEX = 0;
	const unsigned int GCOUNTER_INDEX = 1;
	
	size_t counters [numOfPairs] [2];
	
	// Setup a parallel array of bools that indicates whether a pair of counters is still valid for its
	// corresponding pair. In other words, this states whether or not all monomial multiplications of
	// this pair have been already done.
	bool validCounters [numOfPairs];
	size_t numOfValidPairs = numOfPairs;
	
	// set initial states of counters to 0's, and those of the validCounters array to true
	for (size_t i = 0; i < numOfPairs; i++) {
		counters[i][FCOUNTER_INDEX] = 0;
		counters[i][GCOUNTER_INDEX] = 0;
		
		validCounters[i] = true;
	}
	
	/* Insertions:
	 *	The insertions (for each polynomial pair) can be divided up into two consecutive sets, the 
	 *	'initial' and the 'final' insertions. The way each works out is as follows:
	 *
	 *	Initial insertions:
	 *		Let f = x^3 + x^2 + x + 1 and g = x^2 + x + 1. The initial insertions algorithm sets up
	 *		the counters for the two polynomials to initially be equal 0's. The first step would then
	 *		be to calculate f[0] * g[0] = x^3 * x^2 = x^5. After that, the g counter is increment by 1.
	 *		The second step would then comprise of a loop that multiplies f[0]*g[1] and f[1]*g[0]. Again,
	 *		the g counter is increment (making it equal to 2). And the third would involve multiplying
	 *		f[0]*g[2], f[1]*g[1], and f[2]*g[0]. This would go on until the g counter becomes greater 
	 *		than the length of the f polynomial; in other words, assuming that 'n' equals the length
	 *		of f, the last multiplication set to be performed in the initial set of insertions is:
	 *		f[0]*g[n-1], f[1]*g[n-2], ..., f[n-2]*g[1], and f[n-1]*g[0] (NOTE: since the length of 
	 *		f is greater or equal to the length of g, then g[n-1], and others, might not be valid. This
	 *		would of course result in a monomial equal to 0). After this last step, we move to the second
	 *		set of insertions; the final insertions.
	 *
	 *	Final Insertions:
	 *		In this set, the g counter always starts of as equal to n-1, where n is the length of f (NOTE:
	 *		it is n-1 and not n since array indexing starts from 0, i.e. the last monomial in f is at n-1).
	 *		Also, the first set of insertions in the 'final' insertions sets the f counter to 1. Given the
	 *		above two polynomial examples, the 'final' insertion starts by multiplying the following:
	 *		f[1]*g[n-1], f[2]*g[n-2], ..., f[n-1]*g[1]. After that, the f counter is increment by 1. Then,
	 *		the multiplications of the second set would be: f[2]*g[n-1], f[3]*g[n-2], ..., f[n-1]*g[2].
	 *		This will keep on going until the f counter starts of as larger than the length of f
	 *
	 *	After all these insertions, some monomial might still be in the heap. So loop until the heap is
	 *	empty, and extract them all.
	 */
	
	/* First set of insertions
	 *	This set involves inserting only the largest monomial product of every polynomial pair.
	 *	What this does is it makes every polynomial pair a 'member' of the heap. When an element
	 *	is later extracted from the heap, we can identify from which polynomial pair it was
	 *	calculated, and insert its successor from that pair. In the end, this will lead to all
	 *	monomial products in all polynomial pairs to be inserted in the heap, and thus part of the 
	 *	final result.
	 */
	if (verboseLevel(VERBOSE_LOW)) cout << "Performing the first set of insertions" << endl;
		
	deg_t currDegree = 0;		// to be used anywhere we need to hold the current degree	
	coef_t currCoef = 0;		// to be used anywhere we need to hold the current coefficient
	for (size_t i = 0; i < numOfPairs; i++) {
		currDegree = GET_DEGREE((f_polynomials[i])[0]) + GET_DEGREE((g_polynomials[i])[0]);
		currCoef = GET_COEF((f_polynomials[i])[0]) * GET_COEF((g_polynomials[i])[0]);
		
		if (verboseLevel(VERBOSE_HIGH)) cout << "Size = " << A->size() << ", Inserting degree " << currDegree << ", coef " << currCoef << endl;
		
//		A->insert(currDegree, currCoef, 0);
		A->insert(currDegree, 1, 0);
		
		if (largestSizeReached < A->size())
			largestSizeReached = A->size();
		// pauseProfiling();
		
		// set the IDs and update the counters
		degreeToSource[currDegree].push_back(i);		// add an index to the polynomial pair that generated this degree
		
		counters[i][GCOUNTER_INDEX]++;		// increment the g counter for the next set of insertions
		
		// if the g counter is greater or equal than the length of the f polynomial, then the 'initial' insertions are
		// done for that pair, and we move on to the 'final' insertions set.
		if (counters[i][GCOUNTER_INDEX] >= f_polynomials[i].size()) {
			counters[i][FCOUNTER_INDEX]++;
			counters[i][GCOUNTER_INDEX] = f_polynomials[i].size() - 1;
			
			// if we performed all multiplications for that pair, then set it to be invalid
			if (counters[i][FCOUNTER_INDEX] >= f_polynomials[i].size()) {
				validCounters[i] = false;
				numOfValidPairs--;
			}
		}
		

		
		// continueProfiling();
	}
	
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Done with first set of insertions" << endl << endl;
	
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Performing remaining insertions" << endl;
	// perform all insertions
	bool doneWithInsertions = (numOfValidPairs == 0) ? true : false;		// true when all insertions into the heap have been done

	// REMOVED FOR OPTIMIZATION
	// poly_t currF;		// to be used anywhere we need to hold the current f polynomial
	// poly_t currG;		// to be used anywhere we need to hold the current g polynomial
	
	while (!doneWithInsertions) {
		// start by extracting the maximum
		monom_t maxMonomial = A->poll();
		deg_t maxDegree = GET_DEGREE(maxMonomial);
		coef_t maxCoef = GET_COEF(maxMonomial);
		
		// extract all elements having the same degree as the one to be extracted
		while (!A->isEmpty() && GET_DEGREE(A->peek()) == maxDegree)
			maxCoef += GET_COEF(A->poll());
			
		result.push_back(createMonomial(maxDegree, maxCoef, 0));			// extract the max from the heap and place it in the result
		
		
		if (verboseLevel(VERBOSE_HIGH)) cout << "Size = " << A->size() << ", Extracted degree = " << maxDegree << ", coef " << maxCoef << endl;
		
		// perform the insertion(s) of the successor(s) of the extracted element(s)
		size_t source = 0;
		while (degreeToSource[maxDegree].size() > 0) {		// while there are more pairs associated to this degree
			source = degreeToSource[maxDegree].back();
			degreeToSource[maxDegree].pop_back();
			
			// if pair is all used up, move on to the next
			if (!validCounters[source]) {
				continue;
			}
			
			// get the f and g polynomial of that pair
			poly_t &currF = f_polynomials[source];
			poly_t &currG = g_polynomials[source];
			
			/* Perform 'initial' or 'final' insertions based on the state of that pair */
			
			// if it is in the 'initial' state
			if (counters[source][FCOUNTER_INDEX] == 0) {
				
				size_t f_counter = 0;
				size_t g_counter = counters[source][GCOUNTER_INDEX];
				
				size_t g_length = currG.size();
				
				// skip 'invalid' g's
				if (g_counter >= g_length) {
					g_counter = g_counter - (g_counter - g_length);
					f_counter = f_counter + (g_counter - g_length);
				}

				while (true) {
					if (g_counter < g_length) {
						// insert the new monomial
						currDegree = GET_DEGREE(currF[f_counter]) + GET_DEGREE(currG[g_counter]);
						currCoef = GET_COEF(currF[f_counter]) * GET_COEF(currG[g_counter]);

						if (verboseLevel(VERBOSE_HIGH)) cout << "Initial: Size = " << A->size() << ", Inserting degree " << currDegree << ", coef " << currCoef << endl;
						
//						A->insert(currDegree, currCoef, 0);
						A->insert(currDegree, 1, 0);
						
						if (largestSizeReached < A->size())
							largestSizeReached = A->size();
						
						degreeToSource[currDegree].push_back(source);		// add an index to the polynomial pair that generated this degree
					}

					if (g_counter == 0)
						break;
					else {
						// update values of the counters
						f_counter++;
						g_counter--;
					}

				}
				
				// set the IDs and update the counters
				degreeToSource[currDegree].push_back(source);		// add an index to the polynomial pair that generated this degree

				counters[source][GCOUNTER_INDEX]++;		// increment the g counter for the next set of insertions

				// if the g counter is greater or equal than the length of the f polynomial, then the 'initial' insertions are
				// done for that pair, and we move on to the 'final' insertions set.
				if (counters[source][GCOUNTER_INDEX] >= currF.size()) {
					counters[source][FCOUNTER_INDEX]++;
					counters[source][GCOUNTER_INDEX] = currF.size() - 1;

					// if we performed all multiplications for that pair, then set it to be invalid
					if (counters[source][FCOUNTER_INDEX] >= currF.size()) {
						validCounters[source] = false;
						numOfValidPairs--;
						
						if (verboseLevel(VERBOSE_MED)) cout << "Valid pairs = " << numOfValidPairs << endl;
						
						if (numOfValidPairs == 0) {
							doneWithInsertions = true;
						}
					}
				}
				
				
			}

			// else, if it is in the 'final' state
			else {

				size_t f_counter = counters[source][FCOUNTER_INDEX];
				size_t g_counter = currF.size() - 1;
				
				size_t g_length = currG.size();
				
				// skip 'invalid' g's
				if (g_counter >= g_length) {
					g_counter = g_counter - (g_counter - g_length);
					f_counter = f_counter + (g_counter - g_length);
				}

				size_t g_limit = f_counter;
				while (g_counter >= g_limit) {
					if (g_counter < g_length) {
						// insert the new monomial
						currDegree = GET_DEGREE(currF[f_counter]) + GET_DEGREE(currG[g_counter]);
						currCoef = GET_COEF(currF[f_counter]) * GET_COEF(currG[g_counter]);

						if (verboseLevel(VERBOSE_HIGH)) cout << "Final: Size = " << A->size() << ", Inserting degree " << currDegree << ", coef " << currCoef << endl;
//						A->insert(currDegree, currCoef, 0);
						A->insert(currDegree, 1, 0);
						
						if (largestSizeReached < A->size())
							largestSizeReached = A->size();
							
						degreeToSource[currDegree].push_back(source);		// add an index to the polynomial pair that generated this degree
					}

					f_counter++;
					g_counter--;
				}
				
				// increment the f counter
				counters[source][FCOUNTER_INDEX]++;
				
				if (verboseLevel(VERBOSE_HIGH))cout << "F-Counter = " << counters[source][FCOUNTER_INDEX] << ", f-size = " << currF.size() << endl;
				
				// check if pair is done with final state (i.e. becomes invalid)
				if (counters[source][FCOUNTER_INDEX] >= currF.size()) {
					validCounters[source] = false;
					numOfValidPairs--;
					
					if (verboseLevel(VERBOSE_MED)) cout << "Valid pairs = " << numOfValidPairs << endl;
					
					if (numOfValidPairs == 0) {
						doneWithInsertions = true;
					}
				}
			}
			
		}
	}


	if (verboseLevel(VERBOSE_LOW)) cout << "Done with all insertions" << endl << endl;
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Extracting remaining monomials from the heap" << endl;
	// up till now, all insertions have been done. What remains is to extract the rest of the monomials
	// from the heap and add them to the result
	while (!A->isEmpty()) {
		// start by extracting the maximum
		monom_t maxMonomial = A->poll();
		deg_t maxDegree = GET_DEGREE(maxMonomial);
		coef_t maxCoef = GET_COEF(maxMonomial);
		
		// extract all elements having the same degree as the one to be extracted
		while (!A->isEmpty() && GET_DEGREE(A->peek()) == maxDegree)
			maxCoef += GET_COEF(A->poll());
		
		result.push_back(createMonomial(maxDegree, maxCoef, 0));			// extract the max from the heap and place it in the result
	}
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Done with multiplication process" << endl << endl;
	
	cout << "Largest size reached by the heap is: " << largestSizeReached << endl << endl;
	
	
	delete A;

}




// call procedure depending on the heap chosen
void summationOfProducts (std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result) {
	
	
	if (usingOptimizedSequenceOfInserts()) {
		OSImultiply(f_polynomials, g_polynomials, result);
	}
	else {
		multiplyMultiplePairsFunnelWithMerging (f_polynomials, g_polynomials, result);
	}
}







// parse the polynomial string into a poly_t
void parsePolynomialString (string p, poly_t &result) {
    
    long res_size = getNumberOfMonomials(p.c_str());
    p.append("?");      // ? is a sentinel value
    
	// reserve the size for the result
   	int resultPosition = 0;    
	
	string num = "";
    
    coef_t coef = 0;
    deg_t degree = 0;
	
	size_t i = 0;               // to traverse the poly string]
    char current = p[i++];
    skipSpaces(current, i, p);
    
	while (current != '?') {
        if (current == 'x') {
            coef = 1;
            
            current = p[i++];
            
            if (current == '^') {
                current = p[i++];       // current is (start of) a number
                
                while (current >= '0' && current <= '9') {
                    //num.append(current + "");
                    num += current;
                    current = p[i++];
                }
                
                degree = atoll(num.c_str());
                num.clear();
            }
            else {      // degree is 1
                degree = 1;
            }
            
            skipSpaces(current, i, p);   // go to nearest '+' or end
        }
        else {      // current is a number
            num.clear();
            while (current >= '0' && current <= '9') {
                //num.append(current + "");
                num += current;
                current = p[i++];
                
            }
            
            coef = atoll(num.c_str());
            num.clear();
            
            // coef read; check for degree
            
            if (current == '*') {           // degree exists
                current = p[i++];           // current = x
                
                current = p[i++];           // current is ^ or end of monomial
                
                if (current == '^') {       
                    current = p[i++];       // current is (start of) a number
                    
                    while (current >= '0' && current <= '9') {
                        //num.append(current + "");
                        num += current;
                        current = p[i++];
                    }
                    
                    degree = atoll(num.c_str());
                    num.clear();
                }
                else {      // degree is 1
                    degree = 1;
                }
                
            }
            else {              // degree is 0
                degree = 0;
            }
            
            skipSpaces(current, i, p);
        }
        
        //(*result)[resultPosition++] = coef;
        //(*result)[resultPosition++] = degree;
        //result[resultPosition++] = createMonomial(degree, coef, 0);
		result.push_back(createMonomial(degree, coef, 0));
        
        
        if (current == '?') {       // done reading poly
            break;
        }
        else {      // current is '+'
            current = p[i++];
            skipSpaces(current, i, p);
        }
	}
	
//	(*result)[resultPosition++] = 0;           // add two 0's that act
//    (*result)[resultPosition++] = 0;           // as sentinel values
}


inline void skipSpaces (char &current, size_t &index, string line) {
    while (index <= line.length() && current == ' ') {
        current = line[index++];
    }
}


/*
 getNumberOfMonomials: gets the number of monomials in a (SORTED) encoded polynomial
 */
int getNumberOfMonomials (const char *f) {
	const char *current = f;
	
	int numberOfMonomails = 1;
	while ((*current) != '\0') {
		if ((*current) == '+') {		// found an exponent
			numberOfMonomails++;
		}
		
		current++;
	}
	
	return numberOfMonomails;
}

size_t polynomialSize(const poly_t &p) {
    return p.size();
   // return GET_DEGREE(p[0]);
}
