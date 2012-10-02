//
//  BenchmarkProcedure.cpp
//
//  Created by Karl Gemayel on 09/09/11.
//  Copyright (c) 2012 American University of Beirut. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <time.h>
#include <assert.h>

#include "../BinaryHeap/BinaryHeap.h"
#include "../IOHandling/IOHandling.h"
#include "../Polynomial/Polynomial.h"
#include "../Options/Options.h"

using namespace std;

void usage();
void readArgs(int argc, char *argv[]);
void sum (vector<poly_t> polynomials, poly_t &result);
void addPolynomials (poly_t p1, poly_t p2, poly_t &result);

void iteratedSOP (vector<poly_t> &f_polynomials, vector<poly_t> &g_polynomials, poly_t &result);
void simultaneousSOP (vector<poly_t> &f_polynomials, vector<poly_t> &g_polynomials, poly_t &result);

/****************************************/

// global arguments
string inputName = "mult.in";
bool iterated_sop = false;
bool simult_sop = false;

int main (int argc, char * argv[]) {

	// read/update arguments
	readArgs(argc, argv);
    
    
	vector<poly_t> f_polynomials, g_polynomials;
	int numOfPairs;
	
	readInputFile(inputName, f_polynomials, g_polynomials, numOfPairs);
    
	// choose the binary heap for the operations
	chooseHeap("binary");
	
	
    cout << "Benchmark Procedure (" << numOfPairs << " pairs)" << endl;

	poly_t result;
	
	startTimer();
	
	if (iterated_sop) 
		iteratedSOP (f_polynomials, g_polynomials, result);
	else 
		simultaneousSOP (f_polynomials, g_polynomials, result);
	
	stopTimer();
    
	
	cout << "result size = " << result.size() << endl;

	// compute elapsed time
	cout << "Time = " << getElapsedTime() << endl;   

	// open a file to write the execution time to
	ofstream timesFile;
	string timeFilename;
	getTimeFilename(timeFilename);
	timesFile.open (timeFilename.c_str(), ios::out | ios::app);
	
	// write the input size with the time
	timesFile << numOfPairs << ": " << getElapsedTime() << endl;
	
	// close the times file
	timesFile.close();
	
	
	// write the multiplication output to the output file
	if (printMultOutput()) {
		string result_string;
	    polyToString(result, result_string);
	
		string outputFilename;
		getOutputFilename(outputFilename);
		
		ofstream mult_output;
		mult_output.open (outputFilename.c_str(), ios::out | ios::app);
		mult_output << result_string.c_str() << endl;
		mult_output.close();
	}
}


// read/update the arguments for running the test case
void readArgs(int argc, char *argv[]) {
	// if there are no arguments (except for the program name), then print usage
	if (argc == 1) 
		usage();	
		
	bool choseInputFile = false;
	
	iterated_sop = false;
	simult_sop = false;
	
	// if there are arguments, read them
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-i") == 0) {	// choosing input file
			inputName = argv[++i];
			choseInputFile = !choseInputFile;
		}
		// else if (strcmp(argv[i], "-v") == 0) {
		// 	int vLevel = atoi(argv[++i]);
		// 	setVerbose(vLevel);
		// }
		else if (strcmp(argv[i], "-m") == 0) {
			i++;
			if (strcmp(argv[i], "iterated") == 0) {
				iterated_sop = true;
			}
			else if (strcmp(argv[i], "simultaneous") == 0) {
				simult_sop = true;
			}
		}
		else if (strcmp(argv[i], "-o") == 0) {
			printMultOutput(true);
			setOutputFilename(argv[++i]);
		}
		else if (strcmp(argv[i], "-help") == 0) {
			usage();
		}
		else {
			string invArg = argv[i];
			printError("Invalid argument \"" + invArg + "\"");
			exit(1);
		}
	}
	
	// if any "required" parameter was not given, then print message
	if (!(choseInputFile)) 
		usage();
		
}


/*
 *	Adds two polynomials.
 *
 *	@param p1 the first polynomial
 *	@param p2 the second polynomial
 *	@param result the variable that will hold the result
 */
void addPolynomials (poly_t p1, poly_t p2, poly_t &result) {
	result.clear();		// make sure the result is empty
	
	// if both polynomials are empty, return
	if (p1.size() == 0 && p2.size() == 0) {
		return;
	}
	
	// if one polynomial is empty, return the other
	if (p1.size() == 0) {
		result = p2;
		return;
	}
	else if (p2.size() == 0) {
		result = p1;
		return;
	}
		
	size_t p1_counter = 0;
	size_t p2_counter = 0;
	
	
	// loop until we pass over at least one of the polynomials
	while ((p1_counter < p1.size()) && (p2_counter < p2.size())) {
		monom_t p1_monomial = p1[p1_counter];
		monom_t p2_monomial = p2[p2_counter];
		
		deg_t p1_degree = GET_DEGREE(p1_monomial);
		deg_t p2_degree = GET_DEGREE(p2_monomial);
		
		// if the monomial degrees are equal, add their coefficients
		if (p1_degree == p2_degree) {
			coef_t coef = GET_COEF(p1_monomial) + GET_COEF(p2_monomial);
			result.push_back(createMonomial(p1_degree, coef, 0));
			
			// increment both counters
			p1_counter++;
			p2_counter++;			
		}
		// else if the monomial from p1 has the larger degree, then add it to the result
		else if (p1_degree > p2_degree) {
			result.push_back(p1_monomial);
			
			// increment counter for p1
			p1_counter++;
		}
		// else if the monomial from p2 has the larger degree, then add it to the result
		else {
			result.push_back(p2_monomial);
			
			// increment counter for p2
			p2_counter++;
		}		
	}
	
	// At this point, at least one of the polynomials has been fully traversed. If we have fully 
	// observed both polynomials, we are done. If only one has been fully observed, then we append 
	// the remaining monomials of the ther polynomial to the result.
	bool p1_done = p1_counter == p1.size();
	bool p2_done = p2_counter == p2.size();
	
	assert (p1_done || p2_done);		// assert that at least one is done
	
	poly_t rem_poly;		// the polynomial with the remaining (unobserved) monomials
	size_t rem_counter;		// the counter for rem_poly
	
	// if both are done, quit
	if (p1_done && p2_done) {
		return;
	}
	// if p1 is done but p2 isn't
	else if (p1_done && !p2_done) {
		rem_poly = p2;
		rem_counter = p2_counter;
	}
	// if p1 is not done but p2 is done
	else {
		rem_poly = p1;
		rem_counter = p1_counter;
	}
	
	// append the "remaining" monomials to the result
	while (rem_counter < rem_poly.size()) {
		result.push_back(rem_poly[rem_counter++]);
	}
	
}


/*
 *	Computes the sum of a set of polynomials
 *
 *	@param polynomials the set of polynomials
 *	@param result the variable that will hold the result
 */
void sum (vector<poly_t> polynomials, poly_t &result) {	
	result.clear();		// make sure the result start empty
	
	// setup the vector that will contain the polynomials
	vector<poly_t> currPolynomials = polynomials;
	
	// The summation of all polynmoials will be applied as a binary merger, i.e. keep merging
	// two polynomials together until we reach a single polynomial representing the sum
	
	while (currPolynomials.size() > 1) {
		vector<poly_t> newPolynomials;
			
		size_t counter = 0;
		
		while (counter < currPolynomials.size()) {
			poly_t tempResult;			
			addPolynomials(currPolynomials[counter], currPolynomials[counter+1], tempResult);

			newPolynomials.push_back(tempResult);
			counter += 2;		// increment the counter to indicate the added polynomials
		}
		
		if (newPolynomials.size() % 2 != 0 && newPolynomials.size() > 1) {
			poly_t empty; newPolynomials.push_back(empty);
		}
		currPolynomials = newPolynomials;		
	}
	
	
	assert (currPolynomials.size() == 1);
	
	result = currPolynomials[0];
	
}

void usage () {
	cerr << "usage: Polynomial Arithmetic" << endl <<endl;

	cerr << "Driver: benchmark -m (mode) -i (input_file)" << endl;
	cerr << "Where: -m (mode): the benchmarking mode (see below...)" << endl;
	cerr << "       -i (input_file): reads input from file with name 'input_file'" << endl;
	cerr << endl;
	
	cerr << "Benchmarking Modes:" << endl;
	cerr << "\titerated\tperforms individual pair multiplications, and adds" << endl;
	cerr << "\t         \tthe products as they come out" << endl;
	cerr << "\tsimultaneous\tperforms all multiplications and additions in one" << endl;
	cerr << "\t             \tsingle binary heap" << endl;
	cerr << endl;
	
	cerr << "Optional Arguments:" << endl;
	// cerr << "\t-v (num): puts the program in verbose mode (prints information)" << endl;
	// cerr << "\t\tnum: the level of verbosity (1, 2, or 3)" << endl;
	cerr << "\t-o (output_file): print the output to a file with name 'output_file'" << endl;

	cerr << endl;
	// quit program
	exit(0);
}


/*
 *	Computes the summation of the polynomial products from two sets.
 */
void sop (vector<poly_t> &f_polynomials, vector<poly_t> &g_polynomials, poly_t &result) {
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

	assert (f_polynomials.size() == g_polynomials.size());	// assert that the sets have equal sizes
		
    size_t numOfPairs = g_polynomials.size();
	
	// get maximum number of multiplications that will be done 
	size_t maxMultiplications = 0;
	for (size_t i = 0; i < numOfPairs; i++) {
		maxMultiplications += f_polynomials[i].size() * g_polynomials[i].size();
	}


    size_t currentMultiplications = 0;		// used to track the number of multiplications already done
	// set n_f as the number of monomials in f


    // get the max number of monomials in the f's, and the total heap capacity
	heap_size_t heap_cap = 0;	// capacity of heap
	size_t f_max = 0;
	for (size_t i = 0; i < numOfPairs; i++) {
		heap_cap += f_polynomials[i].size();
		if (f_polynomials[i].size() > f_max) {
			f_max = f_polynomials[i].size();
		}
	}



	if (verboseLevel(VERBOSE_LOW)) cout << "Using Binary Heap" << endl;
	
	BinaryHeap bh (heap_cap);
	Heap &A = bh;


    monom_t f_element, g_element;    // represent an element of f and g



	unsigned long long maxNumberOfElementsInHeap = 0;

	// Insert (f1_1*g1_1, f1_2*g1_2,.... f1_n*g1_n).....(fn_1*gn_1....fn_n*gn_n) into the heap
	int id = 0;
	for (size_t i = 0; i < numOfPairs; i++) {			// loop on all polynomials

		poly_t f = f_polynomials[i];
		poly_t g = g_polynomials[i];

		for (size_t j = 0; j < f_polynomials[i].size(); j++) {		// loop on monomials of i'th polynomial
            f_element = f[j];
			g_element = g[0];

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
	
	int recArray [numOfPairs] [f_max];
	for (size_t i = 0; i < numOfPairs; i++) {
		for (size_t j = 0; j < f_max; j++) {
			recArray[i][j] = 0;
		}
	}


	while (currentMultiplications < maxMultiplications) {
		assert (!A.isEmpty());
		
		monom_t max = A.poll();


		ID_t maxID = GET_ID(max);

		// increment recTable of f(id)
		int maxRow = (int) (maxID / f_max);
		int maxCol = (int) (maxID % f_max);

		// find current reference and update g_counter			
		if (recArray[maxRow][maxCol] == (int) (g_polynomials[maxRow].size() - 1))
			recArray[maxRow][maxCol] = -1;
		else
			recArray[maxRow][maxCol]++;

		coef_t maxCoef = GET_COEF(max);
		deg_t maxDegree = GET_DEGREE(max);

        if (recArray[maxRow][maxCol] != -1) {
			// add the successor
			poly_t f = f_polynomials[maxRow];
			poly_t g = g_polynomials[maxRow];

			int g_successor = recArray[maxRow][maxCol];

            f_element = f[maxCol];
            g_element = g[g_successor];

            deg_t degree = GET_DEGREE(f_element) + GET_DEGREE(g_element);


			A.insert(degree, 1, maxID);

			currentMultiplications++;
		}

		// check if max of heap has degree equal to maxDegree
		monom_t current = A.peek();

		while (!A.isEmpty()) {                              // REVIEW: was current != NULL
			if (GET_DEGREE(current) == maxDegree) {

				current = A.poll();


				long currentID = GET_ID(current);

				int currentRow = (int) (currentID / f_max);
				int currentCol = (int) (currentID % f_max);

				// find current reference and update g_counter					
				if (recArray[currentRow][currentCol] == (int) (g_polynomials[currentRow].size() - 1))
					recArray[currentRow][currentCol] = -1;
				else
					recArray[currentRow][currentCol]++;


				if (recArray[currentRow][currentCol] != -1) {
					// add the successor
					poly_t f = f_polynomials[currentRow];
					poly_t g = g_polynomials[currentRow];

					//int g_successor = recTable[currentRow][currentCol];
					int g_successor = recArray[currentRow][currentCol];

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


		setCoef(max, maxCoef);
        result.push_back(max);
	}

	// extract all remaining elements from heap
	while (!A.isEmpty()) {

		monom_t max = A.poll();

		coef_t maxCoef = GET_COEF(max);
		deg_t maxDegree = GET_DEGREE(max);

		monom_t current = A.peek();

		while (!A.isEmpty()) {
			if (GET_DEGREE(current) == maxDegree) {

				current = A.poll();
				maxCoef += GET_COEF(current);

				current = A.peek();				
			}
			else {
				break;
			}

		}


		setCoef(max, maxCoef);
        result.push_back(max);
	}

}


/*
 *	Performs an iterated version of the summation of products procedure. This is done by computing
 *	the product of the next pair, and immediately adding it to the previous result.
 *	NOTE: The i'th pair is defined to be the i'th polynomial from the first set with the i'th 
 *			polynomial of the second set (this implies that sets should be of equal length).
 *
 *	@param f_polynomials the first set of polynomials
 *	@param g_polynomials the second set of polynomials
 */
void iteratedSOP (vector<poly_t> &f_polynomials, vector<poly_t> &g_polynomials, poly_t &result) {
	assert (f_polynomials.size() == g_polynomials.size());	// assert that sets have equal lengths
	
	size_t numOfPairs = f_polynomials.size();
	
	cout << "Iterated SOP - " << numOfPairs << " pairs" << endl;
	
	for (size_t i = 0; i < numOfPairs; i++) {
		vector<poly_t> tempF; tempF.push_back(f_polynomials[i]);
		vector<poly_t> tempG; tempG.push_back(g_polynomials[i]);
		poly_t temp;
		poly_t tempResult;
		sop(tempF, tempG, temp);
		
		addPolynomials(result, temp, tempResult);
		
		result = tempResult;
	}
}


/*
 *	Performs a "simultaneous" version of the summation of products procedure. This is done using a 
 *	single heap that uses all pairs together.
 */
void simultaneousSOP (vector<poly_t> &f_polynomials, vector<poly_t> &g_polynomials, poly_t &result) {
	cout << "Simultaneous SOP - " << f_polynomials.size() << " pairs" << endl;
	sop (f_polynomials, g_polynomials, result);		// perform all pair multiplication simultaneously
}

