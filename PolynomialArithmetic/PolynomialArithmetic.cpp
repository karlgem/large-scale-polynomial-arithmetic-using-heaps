#include <iostream>
#include <sstream>
#include <string>
#include <assert.h>
#include <vector>
#include "PolynomialArithmetic.h"
#include "../Options/Options.h"

// include heaps
#include "../Heap/Heap.h"
#include "../FunnelHeap/FunnelHeap.h"
#include "../BinaryHeap/BinaryHeap.h"
#include "../FunnelHeapWithMerging/FunnelHeapWithMerging.h"
#include "../BinaryHeapWithChaining/BinaryHeapWithChaining.h"

#define ASSERTIONS false

#define ENABLE_STATE_PRINTING false
#define TURN_ON_PAUSES false
#define ASSERT_REL_BOUNDS false

using namespace std;

/***** Function Prototypes *****/
void sopPackedIDs (std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result);
void sopUnpackedIDs (std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result);
void singleOSImultiply (poly_t &f, poly_t &g, poly_t &result);
void sopOSI (std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result);

void bhWithChainingPreprocessing (std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials);


int p = 113;                      // prime 
void init(int prime) {
    p = prime;
}

/*
 *	Computes the summation of the products of polynomial pairs. This version packs ID's into
 *	each monomial, making the ID management easier, but limits the flexibility of the data
 *	structures that are being used.
 *
 *	@param f_polynomials the set of f polynomials
 *	@param g_polynomials the set of g polynomials
 *	@param result the variable that will hold the result
 */
void sopPackedIDs(std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result) {
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
	 */
    

    size_t numOfPairs = g_polynomials.size();
    
#ifdef VERBOSE
	cout << "Number of pairs = " << numOfPairs << endl;
#endif
	// get maximum number of multiplications that will be done 
	size_t maxMultiplications = 0;
	for (size_t i = 0; i < numOfPairs; i++) {
		maxMultiplications += f_polynomials[i].size() * g_polynomials[i].size();
	}
    
#ifdef VERBOSE
	cout << "Max number of multiplications = " << maxMultiplications << endl;
#endif
    size_t currentMultiplications = 0;		// used to track the number of multiplications already done
	// set n_f as the number of monomials in f
	

    // get the max number of monomials in the f's, and the total heap capacity
	size_t heap_cap = 0;	// capacity of heap
	size_t f_max = 0;
	for (size_t i = 0; i < numOfPairs; i++) {
		heap_cap += f_polynomials[i].size();
		if (f_polynomials[i].size() > f_max) {
			f_max =(size_t) f_polynomials[i].size();
		}
	}
	
	
	if (chosenHeap() == BINARY_HEAP) {
		#define CHOSE_BINARY
	}
	else if (chosenHeap() == FUNNEL_HEAP) {
		#define CHOSE_FUNNEL
	}
	else {
		printError("No heap was initialized! Exiting...");
		exit(1);
	}
	
	#ifdef CHOSE_BINARY
	BinaryHeap bh (heap_cap);
	Heap &A = bh;
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Using Binary Heap" << endl;
	
	#elif CHOSE_FUNNEL
 	FunnelHeap fh (6);
	Heap &A = fh;
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Using Funnel Heap" << endl;
	
	#endif
	

    monom_t f_element, g_element;    // represent an element of f and g
	
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Starting first set of mults" << endl;
	
	unsigned long long maxNumberOfElementsInHeap = 0;
	
	// Insert (f1_1*g1_1, f1_2*g1_2,.... f1_n*g1_n).....(fn_1*gn_1....fn_n*gn_n) into the heap
	int id = 0;
	for (size_t i = 0; i < numOfPairs; i++) {			// loop on all polynomials
		
		poly_t f = f_polynomials[i];
		poly_t g = g_polynomials[i];
		
		for (size_t j = 0; j < f_polynomials[i].size(); j++) {		// loop on monomials of i'th polynomial
            f_element = f[j];
			g_element = g[0];
            
            //coef_t coef = GET_COEF(f_element) * GET_COEF(g_element);
            deg_t degree = GET_DEGREE(f_element) + GET_DEGREE(g_element);
			
			A.insert(degree, 1, id); 
			
			
			maxNumberOfElementsInHeap++;
			id++;
			
			currentMultiplications++;
		}
		
		if (f_polynomials[i].size() < f_max) {			// go to next level in record table
			id += f_max-f_polynomials[i].size();
		}
	}
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Done with first set of mults" << endl;
	
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
		
		for (size_t j = 1; j < f_polynomials[i].size(); j++) {		// create link for every monomials of f_i
			current->next = new ID;
			current = current->next;
			current->g_counter = 0;			// added new
			current->next = NULL;
		}
	}


	if(verboseLevel(VERBOSE_LOW)) cout << "Starting second set of Mults" << endl;
	
	while (currentMultiplications < maxMultiplications) {
		//cout << "i" << endl;
		monom_t max = A.poll();
		//cout << "o" << endl;
		
        
		
		if (max == 0) {
			break;
		}
		
		ID_t maxID = GET_ID(max);
		
		
		
		// increment recTable of f(id)
		int maxRow = (int) (maxID / f_max);
		int maxCol = (int) (maxID % f_max);
		
		// find current reference and update g_counter
		ID *maxRef = recArray[maxRow];
		for (int i = 0; i < maxCol; i++) {
			maxRef = maxRef->next;
		}
		
		if (maxRef->g_counter == (int) g_polynomials[maxRow].size() - 1) 
			maxRef->g_counter = -1;
		else
			maxRef->g_counter++;
		
		coef_t maxCoef = GET_COEF(max);
		deg_t maxDegree = GET_DEGREE(max);
		
        if (maxRef->g_counter != -1) {
			// add the successor
			poly_t f = f_polynomials[maxRow];
			poly_t g = g_polynomials[maxRow];
			
			//			int g_successor = recTable[maxRow][maxCol];
			int g_successor = maxRef->g_counter;
            
            f_element = f[maxCol];
            g_element = g[g_successor];
            
            //coef_t coef = GET_COEF(f_element) * GET_COEF(g_element);
            deg_t degree = GET_DEGREE(f_element) + GET_DEGREE(g_element);
			
			
			A.insert(degree, 1, maxID);
			
			currentMultiplications++;
		}
        
		// check if max of heap has degree equal to maxDegree
		
		//cout << "i" << endl;
		monom_t current = A.peek();
		//cout << "o" << endl;
		
		//cout << "i" << endl;
		while (!A.isEmpty()) {                              // REVIEW: was current != NULL
			//cout << "o" << endl;
			if (GET_DEGREE(current) == maxDegree) {

				current = A.poll();
                
               // result.push_back(current);		// wrong
				
				long currentID = GET_ID(current);
				
				int currentRow = (int) (currentID / f_max);
				int currentCol = (int) (currentID % f_max);
				
				// find current reference and update g_counter
				ID *currentRef = recArray[currentRow];
				for (int i = 0; i < currentCol; i++) {
					currentRef = currentRef->next;
				}
				
				if (currentRef->g_counter == (int) g_polynomials[currentRow].size() - 1) 
					currentRef->g_counter = -1;
				else
					currentRef->g_counter++;
				
				
				if (currentRef->g_counter != -1) {
					// add the successor
					poly_t f = f_polynomials[currentRow];
					poly_t g = g_polynomials[currentRow];
					
					//int g_successor = recTable[currentRow][currentCol];
					int g_successor = currentRef->g_counter;
                    
                    f_element = f[currentCol];
                    g_element = g[g_successor];
                    
                    //coef_t coef = GET_COEF(f_element) * GET_COEF(g_element);
                    deg_t degree = GET_DEGREE(f_element) + GET_DEGREE(g_element);
                    
					A.insert(degree, 1, (int) currentID);
					
					currentMultiplications++;
				}
                
				
				maxCoef += GET_COEF(current);
				
				current = A.peek();
			}
			else {
				break;
			}
			
		}
        
        // if (maxCoef % p == 0) {
        //     continue;
        // }
        
        	
//        setCoef(max, maxCoef % p);
		setCoef(max, maxCoef);
        result.push_back(max);
	}
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Done with second set of Mults" << endl << endl;	

	if (verboseLevel(VERBOSE_LOW)) cout << "Starting third set of mults" << endl;
	
	// extract all remaining elements from heap
	while (!A.isEmpty()) {
		
		if (A.size() == 0) {
			cout << "STUPID" << endl;
		}
		
		monom_t max = A.poll();
        		
		coef_t maxCoef = GET_COEF(max);
		deg_t maxDegree = GET_DEGREE(max);
		
		monom_t current = A.peek();

		while (!A.isEmpty()) {
			if (A.size() == 0) cout << "Size = " << A.size() << endl;
			
			if (GET_DEGREE(current) == maxDegree) {
				
				current = A.poll();
				maxCoef += GET_COEF(current);
                
				current = A.peek();				
			}
			else {
				break;
			}
			
		}
        
        // if (maxCoef % p == 0) {
        //     continue;
        // }
        
//        setCoef(max, maxCoef % p);
		setCoef(max, maxCoef);
        result.push_back(max);
	}
	
	//cout << "o" << endl;
	if (verboseLevel(VERBOSE_LOW)) cout << "Done with third set of mults" << endl;
    
	
}



#ifdef COMPILE_MERGING
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
		cout << "\\* Starting Multiplication Procedure \"FH-Merging\"! (" << g_polynomials.size() << " pairs) */" << endl;


    size_t numOfPairs = g_polynomials.size();
    
	// get maximum number of multiplications that will be done 
	size_t maxMultiplications = 0;
	for (size_t i = 0; i < numOfPairs; i++) {
		maxMultiplications += f_polynomials[i].size() * g_polynomials[i].size();
	}
    
    size_t currentMultiplications = 0;		// used to track the number of multiplications already done
	// set n_f as the number of monomials in f
	

    // get the max number of monomials in the f's, and the total heap capacity
	size_t heap_cap = 0;	// capacity of heap
	size_t f_max = 0;
	for (size_t i = 0; i < numOfPairs; i++) {
		heap_cap += f_polynomials[i].size();
		if (f_polynomials[i].size() > f_max) {
			f_max =(size_t) f_polynomials[i].size();
		}
	}
	

	bool choseFunnelWithMerging = false;
	bool choseBinary = false;
	bool choseFunnel = false;
	bool choseBinaryWithChaining = false;

	if (chosenHeap() == FUNNEL_HEAP_WITH_MERGING) {
		choseFunnelWithMerging = true;
	}
	else if (chosenHeap() == BINARY_HEAP) {
		choseBinary = true;
	}
	else if (chosenHeap() == FUNNEL_HEAP) {
		choseFunnel = true;
	}
	else if (chosenHeap() == BINARY_HEAP_WITH_CHAINING) {
		choseBinaryWithChaining = true;
	}
	else {
		printError("No heap was initialized! Exiting...");
		exit(1);
	}
	
	
	Heap *A = NULL;
	
	if (choseBinaryWithChaining) {	
		// for the binaryheap with chaining, we first need to perform some preprocessing
		// to determine some characteristics of the chains. Thie preprocessing should not
		// be included in the profiling sections.

		pauseTimer();
		bhWithChainingPreprocessing(f_polynomials, g_polynomials);
		continueTimer();
		
		// allocate heap
//		BinaryHeapWithChaining bhc (heap_cap);
		//	Heap &A = bhc;
		A = new BinaryHeapWithChaining (heap_cap);
	
	
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Binary Heap (with chaining)" << endl;
	}
	else if (choseFunnelWithMerging) {
	
//		FunnelHeapWithMerging fhm (6);
//		Heap &A = fhm;
		A = new FunnelHeapWithMerging (6);
	
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Funnel Heap (with merging)" << endl;
	
	}
	else if (choseBinary) {
		// BinaryHeap bh (heap_cap);
//		Heap &A = bh;
		A = new BinaryHeap (heap_cap);
	
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Binary Heap (capacity = " << heap_cap << ")" << endl;
	}
	else if (choseFunnel) {
	 // FunnelHeap fh (6);
	//	Heap &A = fh;
		A = new FunnelHeap (6);
	
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Funnel Heap" << endl;
	}
	
    
    monom_t f_element, g_element;    // represent an element of f and g
	
	
	if (verboseLevel(VERBOSE_MED)) cout << "Retrieving the largest degree of a monomial multiplication..." << endl;
	
	// get the largest degree that can be produced by a mutliplication of f and g
	deg_t largestDegree = 0;
	deg_t possibleLargestDegree = 0;
	for (size_t i = 0; i < numOfPairs; i++) {
		possibleLargestDegree = GET_DEGREE((f_polynomials[i])[0]) + GET_DEGREE((g_polynomials[i])[0]);
		if (largestDegree < possibleLargestDegree) {
			largestDegree = possibleLargestDegree;
		}
	}

	if (verboseLevel(VERBOSE_MED)) cout << "Largest degree retrieved: " << largestDegree <<  endl;

	

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
	for (size_t i = 0; i < numOfPairs; i++) {			// loop on all polynomials
		
		poly_t f = f_polynomials[i];
		poly_t g = g_polynomials[i];
		
		for (size_t j = 0; j < f_polynomials[i].size(); j++) {		// loop on monomials of i'th polynomial
            f_element = f[j];
			g_element = g[0];
            
            //coef_t coef = GET_COEF(f_element) * GET_COEF(g_element);
            deg_t degree = GET_DEGREE(f_element) + GET_DEGREE(g_element);
			
			if (verboseLevel(VERBOSE_HIGH)) 
				cout << "SIZE=" << A->size() << "\tInserting: coef=" << 1 << ", deg=" << degree << endl;
				
			A->insert(degree, 1);
			
			
			if (ENABLE_STATE_PRINTING) {
				cout << "Inserted " << degree << endl << endl;
				// A->print();
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
		
		for (size_t j = 1; j < f_polynomials[i].size(); j++) {		// create link for every monomials of f_i
			current->next = new ID;
			current = current->next;
			current->g_counter = 0;			// added new
			current->next = NULL;
		}
	}
		
	if (verboseLevel(VERBOSE_LOW)) cout << "Performing second set of insertions..." << endl;
	
	
	while (currentMultiplications < maxMultiplications) {	
		monom_t max = A->poll();
		
		
		if (ENABLE_STATE_PRINTING) {
			cout << "Removed " << GET_DEGREE(max) << endl << endl;
			A->print();
			pausePrompt();
		}
        
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

			if (maxRef->g_counter == (int) g_polynomials[maxRow].size() - 1) 
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

	            //coef_t coef = GET_COEF(f_element) * GET_COEF(g_element);
	            deg_t degree = GET_DEGREE(f_element) + GET_DEGREE(g_element);

//				if (coef % p != 0) {
					if (verboseLevel(VERBOSE_HIGH)) 
						cout << "SIZE=" << A->size() << "\tInserting: coef=" << 1 << ", deg=" << degree << endl;
						
					A->insert (degree, 1);
					
					
					if (ENABLE_STATE_PRINTING) {
						cout << "Inserted " << degree << endl << endl;
						A->print();
					 	pausePrompt();
					}
					
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
		

		if (ENABLE_STATE_PRINTING) {
			A->print();
			pausePrompt();
			cout << "Removed " << GET_DEGREE(max) << endl << endl;
		}
        		
        result.push_back(max);
	}
    
	delete A;
	
}

#endif


/*
 *	Computes the summation of the products of polynomial pairs. This version does not pack ID's
 *	into each monomial, but rather manages ID's in an external data structure. This data structure
 *	is slightly more complex than when ID's are packed, but it does allow for optimized versions
 *	of the data structures.
 *
 *	@param f_polynomials the set of f polynomials
 *	@param g_polynomials the set of g polynomials
 *	@param result the variable that will hold the result
 */
void sopUnpackedIDs(std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result) {
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
		cout << "\\* Starting Multiplication Procedure \"NO ID\"! (" << g_polynomials.size() << " pairs) */" << endl;


    size_t numOfPairs = g_polynomials.size();
    
	// get maximum number of multiplications that will be done 
	size_t maxMultiplications = 0;
	for (size_t i = 0; i < numOfPairs; i++) {
		maxMultiplications += f_polynomials[i].size() * g_polynomials[i].size();
	}
    
    size_t currentMultiplications = 0;		// used to track the number of multiplications already done
	// set n_f as the number of monomials in f
	

    // get the max number of monomials in the f's, and the total heap capacity
	size_t heap_cap = 0;	// capacity of heap
	size_t f_max = 0;
	for (size_t i = 0; i < numOfPairs; i++) {
		heap_cap += f_polynomials[i].size();
		if (f_polynomials[i].size() > f_max) {
			f_max =(size_t) f_polynomials[i].size();
		}
	}
	
	bool choseFunnelWithMerging = false;
	bool choseBinary = false;
	bool choseFunnel = false;
	bool choseBinaryWithChaining = false;

	if (chosenHeap() == FUNNEL_HEAP_WITH_MERGING) {
		choseFunnelWithMerging = true;
	}
	else if (chosenHeap() == BINARY_HEAP) {
		choseBinary = true;
	}
	else if (chosenHeap() == FUNNEL_HEAP) {
		choseFunnel = true;
	}
	else if (chosenHeap() == BINARY_HEAP_WITH_CHAINING) {
		choseBinaryWithChaining = true;
	}
	else {
		printError("No heap was initialized! Exiting...");
		exit(1);
	}
	
	
	Heap *A = NULL;
	
	if (choseBinaryWithChaining) {	
		// for the binaryheap with chaining, we first need to perform some preprocessing
		// to determine some characteristics of the chains. Thie preprocessing should not
		// be included in the profiling sections.

		pauseTimer();
		bhWithChainingPreprocessing(f_polynomials, g_polynomials);
		continueTimer();
		
		// allocate heap
//		BinaryHeapWithChaining bhc (heap_cap);
		//	Heap &A = bhc;
		A = new BinaryHeapWithChaining (heap_cap);
	
	
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Binary Heap (with chaining)" << endl;
	}
	else if (choseFunnelWithMerging) {
	
//		FunnelHeapWithMerging fhm (6);
//		Heap &A = fhm;
		A = new FunnelHeapWithMerging (6);
	
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Funnel Heap (with merging)" << endl;
	
	}
	else if (choseBinary) {
		// BinaryHeap bh (heap_cap);
//		Heap &A = bh;
		A = new BinaryHeap (heap_cap);
	
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Binary Heap (capacity = " << heap_cap << ")" << endl;
	}
	else if (choseFunnel) {
	 // FunnelHeap fh (6);
	//	Heap &A = fh;
		A = new FunnelHeap (6);
	
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Funnel Heap" << endl;
	}
	
    
    monom_t f_element, g_element;    // represent an element of f and g
	
	
	if (verboseLevel(VERBOSE_MED)) cout << "Retrieving the largest degree of a monomial multiplication..." << endl;
	
	// get the largest degree that can be produced by a mutliplication of f and g
	deg_t largestDegree = 0;
	deg_t possibleLargestDegree = 0;
	for (size_t i = 0; i < numOfPairs; i++) {
		possibleLargestDegree = GET_DEGREE((f_polynomials[i])[0]) + GET_DEGREE((g_polynomials[i])[0]);
		if (largestDegree < possibleLargestDegree) {
			largestDegree = possibleLargestDegree;
		}
	}

	if (verboseLevel(VERBOSE_MED)) cout << "Largest degree retrieved" << endl;

	

	if (verboseLevel(VERBOSE_MED)) cout << "Allocating id-monomial product mapping array.." << endl;

	// create an array of vectors that would hold a mapping of a monomial inserted
	// into the heap to the actual f and g monomials whose product is that monomial
	vector<unsigned long> ID_Monom_pair [largestDegree+1];
	

	if (verboseLevel(VERBOSE_MED)) cout << "ID-monomial product mapping array allocated" << endl;

	

	if (verboseLevel(VERBOSE_LOW)) cout << "\nPerforming the initial insertions..." << endl;
	
	unsigned long long maxNumberOfElementsInHeap = 0;
	// Insert (f1_1*g1_1, f1_2*g1_2,.... f1_n*g1_n).....(fn_1*gn_1....fn_n*gn_n) into the heap
	int id = 0;
	for (size_t i = 0; i < numOfPairs; i++) {			// loop on all polynomials
		
		poly_t f = f_polynomials[i];
		poly_t g = g_polynomials[i];
		
		for (size_t j = 0; j < f_polynomials[i].size(); j++) {		// loop on monomials of i'th polynomial
            f_element = f[j];
			g_element = g[0];
            
            //coef_t coef = GET_COEF(f_element) * GET_COEF(g_element);
            deg_t degree = GET_DEGREE(f_element) + GET_DEGREE(g_element);
			
			if (verboseLevel(VERBOSE_HIGH)) 
				cout << "SIZE=" << A->size() << "\tInserting: coef=" << 1 << ", deg=" << degree << endl;
				
			A->insert(degree, 1);
			
			// PAUSE Timer: to access ID-monomial relationship table
			pauseTimer();
			
			// record the id-monomial relationship
//			stack<unsigned long> rel = ID_Monom_pair[degree];
			vector<unsigned long> *rel = &ID_Monom_pair[degree];
			
			assert (rel != NULL);		// sanity check
			//rel.push(id);
			rel->push_back(id);
			
			// CONTINUE Timer: finished accessing ID-monomial relationship table
			continueTimer();
						
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
		
		for (size_t j = 1; j < f_polynomials[i].size(); j++) {		// create link for every monomials of f_i
			current->next = new ID;
			current = current->next;
			current->g_counter = 0;			// added new
			current->next = NULL;
		}
	}
		
	if (verboseLevel(VERBOSE_LOW)) cout << "Performing second set of insertions..." << endl;
	
	while (currentMultiplications < maxMultiplications) {	
		if (verboseLevel(VERBOSE_HIGH)) cout << "Mult # " << currentMultiplications << "/" << maxMultiplications << endl;	
		monom_t max = A->peek();
        
		if (max == 0) {
			break;
		}
		
		coef_t maxCoef = 0;
		deg_t maxDegree = GET_DEGREE(max);
		
		// PAUSE Timer: to access ID-monomial relationship table
		pauseTimer();
					
		// using the table defining the relationship between ID's and monomials
		// inserted into the heap, check the degree of this max monomial to obtain
		// from what monomial products it resulted.		
		
		// first get the number of ids that correspond to this max
		vector<unsigned long> *relation = &ID_Monom_pair[maxDegree];
		size_t numOfIds = relation->size();
		
		// CONTINUE Timer: finished accessing ID-monomial relationship table
		continueTimer();
				
		// second, extract all the monomials of that degree and add them to the result
		while (!A->isEmpty() && maxDegree == GET_DEGREE(A->peek())) {
			// extract the max element from the heap, and insert its successor
			max = A->poll();
			
			// add the coefficient of the extracted monomial with the previous coefficients
			// of the monomials of that same degree
			maxCoef += GET_COEF(max);
			
			if (verboseLevel(VERBOSE_HIGH)) 
				cout << "SIZE=" << A->size() << "\tExtracted: coef=" << maxCoef << ", deg=" << maxDegree << endl;
		}
		
		// now insert the successor for each of the pairs that generated monomials of having maxDegree
		for (size_t id_counter = 0; id_counter < numOfIds; id_counter++) {
			
			// PAUSE Timer: to access ID-monomial relationship table
			pauseTimer();

			// get one of the ids corresponding to this relation
			unsigned long currentMaxID = relation->back();
			relation->pop_back();
			
			// CONTINUE Timer: finished accessing ID-monomial relationship table
			continueTimer();

			// increment recTable of f(id)
			int maxRow = (int) (currentMaxID / f_max);
			int maxCol = (int) (currentMaxID % f_max);

			// find current reference and update g_counter
			ID *maxRef = recArray[maxRow];
			for (int i = 0; i < maxCol; i++) {
				maxRef = maxRef->next;
			}

			if (maxRef->g_counter == (int) g_polynomials[maxRow].size() - 1) 
				maxRef->g_counter = -1;
			else
				maxRef->g_counter++;



	        if (maxRef->g_counter != -1) {
				// add the successor
				poly_t f = f_polynomials[maxRow];
				poly_t g = g_polynomials[maxRow];

	            f_element = f[maxCol];
	            g_element = g[maxRef->g_counter];

	            //coef_t coef = GET_COEF(f_element) * GET_COEF(g_element);
	            deg_t degree = GET_DEGREE(f_element) + GET_DEGREE(g_element);


				if (verboseLevel(VERBOSE_HIGH)) 
					cout << "SIZE=" << A->size() << "\tInserting: coef=" << 1 << ", deg=" << degree << endl;
					
				A->insert (degree, 1);

				
				// record the id-monomial relationship				
				pauseTimer(); 				// PAUSE Timer: to access ID-monomial relationship table

				vector<unsigned long> *rel = &ID_Monom_pair[degree];
				rel->push_back(currentMaxID);
				
				continueTimer();	// CONTINUE Timer: finished accessing ID-monomial relationship table
				
		
				currentMultiplications++;
			}
			
		}
		

		setCoef(max, maxCoef);
		if (verboseLevel(VERBOSE_HIGH)) cout << "Add to result: " << GET_DEGREE(max) << endl;
        result.push_back(max);
	}
		
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Done with second set of insertions" << endl << endl;
	
	if (verboseLevel(VERBOSE_LOW)) cout << "Extracting all remaining elements from heap" << endl;

	// extract all remaining elements from heap
	while (!A->isEmpty()) {
		monom_t max = A->poll();
        		
		coef_t maxCoef = GET_COEF(max);
		deg_t maxDegree = GET_DEGREE(max);
		
		while (!A->isEmpty() && (GET_DEGREE(A->peek()) == maxDegree)) {
			maxCoef += GET_COEF(A->poll());
		}
		
		setCoef(max, maxCoef);
        result.push_back(max);
	}
	
	// free the heap
	delete A;
    
	
}



/*
 *	The "controller" method that decides which "summation of products" procedure should be used,
 *	based on some set options.
 *
 *	@param f_polynomials the set of f polynomials
 *	@param g_polynomials the set of g polynomials
 *	@param result the variable that will hold the result
 */
void summationOfProducts (std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result) {
	
	
	if (usingOptimizedSequenceOfInserts()) {
		sopOSI(f_polynomials, g_polynomials, result);
	}
	else {		// if using non-optimized sequence of inserts
		if (chosenHeap() == FUNNEL_HEAP_WITH_MERGING) {
			//multiplyMultiplePairsFunnelWithMerging (f_polynomials, g_polynomials, result);
			sopUnpackedIDs(f_polynomials, g_polynomials, result);
		}
		else if (chosenHeap() == BINARY_HEAP) {
			sopUnpackedIDs(f_polynomials, g_polynomials, result);
		}
		else if (chosenHeap() == FUNNEL_HEAP) {
			sopUnpackedIDs(f_polynomials, g_polynomials, result);
		}
		else if (chosenHeap() == BINARY_HEAP_WITH_CHAINING) {
			//multiplyMultiplePairsFunnelWithMerging (f_polynomials, g_polynomials, result);
			sopUnpackedIDs(f_polynomials, g_polynomials, result);
		}
	}
}





void bhWithChainingPreprocessing (std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials) {
	/* CODE SKETCH
	 *	1) Setup two (2D) arrays of numbers, each of size equal to the numbers of pairs of polynomials
	 *		a) The first array holds counters for the first set of polynomials, and the second array
	 *			holds counters for the second set of polynomials.
	 *		b) The i'th counter in the first array indexes the i'th monomial in the first set of 
	 *			polynomials. Similarly for the second array.
	 *		c) The first column of the arrays would hold the actual counters, while the second columns
	 *			would hold the maximum value for each counter (i.e. the size of the polynomial that will
	 *			be indexed by this counter)
	 *	2) Setup another array having a size equal to the largest degree (plus 1). This array (initialized
	 *		to 0's) will be modified accordingly to represent the number of monomials of a certain
	 *		degree. For instance, the element at index 5 will hold a number showing the number of monomial
	 *		products that will have a degree equal to 5.
	 *	3) Write the information to a file: 
	 *		The info file contains lines of numbers. The first line contains the number of chains	
	 *		that will need in the whole multiplication process. The remaining lines each contain 
	 *		two numbers seperated by the space character. The first of the two numbers represents 
	 *		the degree for which the chain will be created, and the second number represents the 
	 *		capacity of the chain of that degree.
	 *		NOTE: the lines containing degree-capacity pairs are MUST be sorted in increasing order of
	 *			degrees, and must contain EVERY degree, starting from 0 up to the largest degree used
	 */
	
	const unsigned int COUNTER_INDEX = 0;
	const unsigned int MAX_INDEX = 1;
	
	size_t numberOfPairs = f_polynomials.size();
	// setup the two arrays of counters
	size_t fCounters [numberOfPairs][2];
	size_t gCounters [numberOfPairs][2];
	
	// initialize all counters to 0's, and the max values to their corresponding values
	for (size_t i = 0; i < numberOfPairs; i++) {
		// initialize counters
		fCounters[i][COUNTER_INDEX] = 0;
		gCounters[i][COUNTER_INDEX] = 0;
		
		
		// initialize the max values
		fCounters[i][MAX_INDEX] = (f_polynomials[i]).size();
		gCounters[i][MAX_INDEX] = (g_polynomials[i]).size();
	}
	
	
	// get the largest degree that can be generated from all monomial products
	deg_t maxDegree = 0;
	deg_t tempDegree = 0;
	for (size_t i = 0; i < numberOfPairs; i++) {
		tempDegree = GET_DEGREE((f_polynomials[i])[0]) + GET_DEGREE((g_polynomials[i])[0]);
		
		if (maxDegree < tempDegree) 
			maxDegree = tempDegree;
	}
	
	
	// setup the array that will contain the number of monomial products for each degree
	size_t numberOfMonomials [maxDegree+1];		// the extra '1' is for the monomial of degree 0
	
	// initialize all 'numberOfMonomials' elements to 0's
	for (size_t i = 0; i <= maxDegree; i++) {
		numberOfMonomials[i] = 0;
	}
	
	
	/* Find the number of monomial products that will be generated for each degree */
	
	// for each polynomial pair, do the following:
	// 		for each monomial in the f-polynomial, multiply its degree
	//		by every degree of the monomials in the g_polynomial, and for
	//		each degree product, increment the number of monomials for that
	// 		degree in the array holding the number of monomials
	
	poly_t f_i, g_i;
	for (size_t i = 0; i < numberOfPairs; i++) {
		
		// get the i'th f and g polynomials
		f_i = f_polynomials[i];
		g_i = g_polynomials[i];
		
		// loop on the monomials of the f-polynomial, and multiply their degrees
		// with those of the current g-polynomial
		deg_t currDegree = 0;
		for (size_t j = 0; j < f_i.size(); j++) {
			currDegree = GET_DEGREE(f_i[j]);
			
			// loop on the monomials of the g-polynomial
			for (size_t k = 0; k < g_i.size(); k++) {
				
				// increment the number of monomials for the degree product
				numberOfMonomials[currDegree + GET_DEGREE(g_i[k])]++;
			}
		}
	}
	
	/* Write the preprocessing information to the info file */
	
	// open the file to write the preprocessing info to
	string cInfoFilename;
	getBHWithChainingFilename(cInfoFilename);		// get the filename of the info file
	
	ofstream cInfo (cInfoFilename.c_str());		// open the info file
	
	// the first line of the file is the total number of chains
	cInfo << maxDegree+1 << endl;
	
	// the remaining lines contain the degree-numberOfMonomials pair, sorted in increasing 
	// order, containing ALL degrees starting from 0 up to 'maxDegree'
	for (size_t i = 0; i <= maxDegree; i++) {
		cInfo << i << " " << numberOfMonomials[i];
		
		if (i < maxDegree)
			cInfo << endl;
	}
	
	
	cInfo.close();					// close the info file
}




// multiply using an optimized sequence of inserts
void sopOSI(std::vector<poly_t> &f_polynomials, std::vector<poly_t> &g_polynomials, poly_t &result) {
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
	

	
	
	// get the max number of monomials in the f's, and the total heap capacity
	size_t heap_cap = 0;
	for (size_t i = 0; i < numOfPairs; i++) {
		for (size_t j = 0; j < f_polynomials[i].size(); j++) {
			heap_cap += 2*j;
		}
	}
	// for (size_t i = 0; i < numOfPairs; i++) {
	// 	heap_cap += GET_DEGREE((f_polynomials[i])[0]) * GET_DEGREE((g_polynomials[i])[0]) + 1;		// REVIEW
	// }
	
	
	// get the heap to be used in computation
	Heap *A = NULL;

	if (chosenHeap() == FUNNEL_HEAP_WITH_MERGING) {
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Funnel Heap (with merging)" << endl;
		
		// allocate heap
		A = new FunnelHeapWithMerging (6);		
	}
	else if (chosenHeap() == BINARY_HEAP) {
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Binary Heap (capacity = " << heap_cap << ")" << endl;
		
		A = new BinaryHeap (heap_cap);		
	}
	else if (chosenHeap() == FUNNEL_HEAP) {
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Funnel Heap" << endl;
		
		// allocate heap
		A = new FunnelHeap (6);
	}
	else if (chosenHeap() == BINARY_HEAP_WITH_CHAINING) {
		// cout << "The method singleOSImultiply() is not yet compatible with Binary Heap with Chaining" << endl;
		// exit(0);
		
		// for the binaryheap with chaining, we first need to perform some preprocessing
		// to determine some characteristics of the chains. Thie preprocessing should not
		// be included in the profiling sections.

		pauseTimer();
		bhWithChainingPreprocessing(f_polynomials, g_polynomials);
		continueTimer();

		if (verboseLevel(VERBOSE_LOW)) cout << "Using Binary Heap (with chaining)" << endl;
		
		// allocate heap
		A = new BinaryHeapWithChaining (heap_cap);				
	}
	else {
		printError("No heap was initialized! Exiting...");
		exit(1);
	}
	
	
	/* IDs and Counters Management
	 *	- Id's are used to link every monomial product to its orignal polynomial pair. This is 
	 *		needed to identify the source of the successor that needs to be inserted
	 *	- Counters are used to preserve the 'state' of every polynomial pair in terms of where
	 *		it had last reached in its multiplication process. Two counters will be assigned
	 *		to each pair (one for the f polynomial, and the other for the g). In simpler terms,
	 *		the counters will be used to identify the successor of the extracted monomial product. 
	 *		
	 * How it all works:
	 *	Once a monomial is extracted, identify the polynomial pair from which it was calculated
	 *	(using ID's). Then, associate the polynomial pair with its respective 'counters'. Proceed
	 *	to perform the calculation (according to the optimized sequence of inserts algorithm) based
	 *	on the state defined by these two counters.
	 */
	
	// calculate the largest degree
	deg_t largestDegree = 0;
	deg_t tempDegree = 0;
	for (size_t i = 0; i < numOfPairs; i++) {
		tempDegree = GET_DEGREE((f_polynomials[i])[0]) + GET_DEGREE((g_polynomials[i])[0]);
		if (largestDegree < tempDegree)
			largestDegree = tempDegree;
	}
	
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

	delete A;
}



// multiplies two polynomials using an optimized sequence of inserts
void singleOSImultiply (poly_t &f, poly_t &g, poly_t &result) {
	if (verboseLevel(VERBOSE_LOW)) 
		cout << "Procedure: Polynomial Pair Multiplication using Optimized Sequence of Inserts" << endl;

   
	

    // get the max number of monomials in the f's, and the total heap capacity
	size_t heap_cap = GET_DEGREE(f[0]) * GET_DEGREE(g[0]);	// REVIEW
	
	
	Heap *A = NULL;

	if (chosenHeap() == FUNNEL_HEAP_WITH_MERGING) {
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Funnel Heap (with merging)" << endl;
		
		// allocate heap
		A = new FunnelHeapWithMerging (6);		
	}
	else if (chosenHeap() == BINARY_HEAP) {
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Binary Heap (capacity = " << heap_cap << ")" << endl;
		
		A = new BinaryHeap (heap_cap);		
	}
	else if (chosenHeap() == FUNNEL_HEAP) {
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Funnel Heap" << endl;
		
		// allocate heap
		A = new FunnelHeap (6);
	}
	else if (chosenHeap() == BINARY_HEAP_WITH_CHAINING) {
		cout << "The method singleOSImultiply() is not yet compatible with Binary Heap with Chaining" << endl;
		exit(0);
	}
	else {
		printError("No heap was initialized! Exiting...");
		exit(1);
	}
	
	
	
	// get the lengths of the two polynomials
	size_t f_length = f.size();
	size_t g_length = g.size();
	
	// switch f and g so that f is the polynomial with more monomials
	if (f_length < g_length) {
		// swap f and g
		poly_t temp = g;
		g = f;
		f = temp;
		
		// update the new lengths
		f_length = f.size();
		g_length = g.size();
	}
	
	
	// set counters and point into the polynomials
	size_t f_counter = 0;
	size_t g_counter = 0;
	
	for (size_t i = 0; i < f_length; i++) {
		if (i != 0) {
			deg_t degreeExtract = GET_DEGREE(A->peek());
			deg_t coefExtract = GET_COEF(A->peek());
			
			A->poll();
			
			// extract all elements having the same degree as the one to be extracted
			while (!A->isEmpty() && GET_DEGREE(A->peek()) == degreeExtract)
				coefExtract += GET_COEF(A->poll());
				
			result.push_back(createMonomial(degreeExtract, coefExtract, 0));			// extract the max from the heap and place it in the result
		}
		
		f_counter = 0;
		g_counter = i;
		
		while (true) {
			if (g_counter < g_length) {
				// insert the new monomial
				deg_t monomDegree = GET_DEGREE(f[f_counter]) + GET_DEGREE(g[g_counter]);
				coef_t monomCoef = GET_COEF(f[f_counter]) * GET_COEF(g[g_counter]);
				
				A->insert(monomDegree, monomCoef, 0);
			}
			
			if (g_counter == 0)
				break;
			else {
				// update values of the counters
				f_counter++;
				g_counter--;
			}
				
		}
		
	}
	
	
	for (size_t i = 1; i < f_length; i++) {
		
		deg_t degreeExtract = GET_DEGREE(A->peek());
		deg_t coefExtract = GET_COEF(A->peek());
		
		A->poll();
		
		// extract all elements having the same degree as the one to be extracted
		while (!A->isEmpty() && GET_DEGREE(A->peek()) == degreeExtract)
			coefExtract += GET_COEF(A->poll());
			
		result.push_back(createMonomial(degreeExtract, coefExtract, 0));			// extract the max from the heap and place it in the result
		
		f_counter = i;
		g_counter = f_length - 1;
		
		while (g_counter >= i) {
			if (g_counter < g_length) {
				// insert the new monomial
				deg_t monomDegree = GET_DEGREE(f[f_counter]) + GET_DEGREE(g[g_counter]);
				coef_t monomCoef = GET_COEF(f[f_counter]) * GET_COEF(g[g_counter]);
				
				A->insert(monomDegree, monomCoef, 0);
			}
			
			f_counter++;
			g_counter--;
		}
	}
	
	while (!A->isEmpty()) {
		deg_t degreeExtract = GET_DEGREE(A->peek());
		deg_t coefExtract = GET_COEF(A->peek());

		A->poll();

		// extract all elements having the same degree as the one to be extracted
		while (!A->isEmpty() && GET_DEGREE(A->peek()) == degreeExtract)
			coefExtract += GET_COEF(A->poll());
			
		result.push_back(createMonomial(degreeExtract, coefExtract, 0));			// extract the max from the heap and place it in the result
	}
}






