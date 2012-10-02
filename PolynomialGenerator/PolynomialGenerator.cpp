// #include <stdafx.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <functional>
#include <vector>
#include <ctime>
#include <cstdlib>

#include "PolynomialGenerator.h"
#include "../Options/Options.h"
using namespace std;

/***** Function Prototypes *****/

// void readArgs (int argc, char *argv[]);
// void usage();
void generateSinglePolynomial (int degree, string &polyString);
int computeNumberOfPairsNeeded1(size_t maxSizeOfQueue);


// specify default values
const int DEFAULT_NUMBER_OF_PAIRS = 10;
// const string DEFAULT_OUTPUT_FILE_NAME = "output.out";
const float DEFAULT_SPARSITY = 0.5;

// global variables
int numOfPairs = DEFAULT_NUMBER_OF_PAIRS;
// string outputFilename = DEFAULT_OUTPUT_FILE_NAME;
float sparsity = DEFAULT_SPARSITY;


// sets the highest possible value for a coefficient
const int COEF_LIMIT = 5;


/***** Function Implementations *****/

void printVector (std::vector<int> &v) {
    for (size_t i = 0; i < v.size(); i++) {
        cout << v[i] << " ";
    }
    
    cout << endl;
}



// int main (int argc, char *argv []) {
//     
// 	readArgs(argc, argv);
// 	
//     generatePolynomials ();
// }
// 

void generatePolynomials (unsigned int n, double s, const string filename) {
	numOfPairs = n;
	sparsity = s;
	
    // set the output stream
    ofstream outputFile;
    outputFile.open (filename.c_str(), ios::out);
    
    // write the number of pairs to the output file
    outputFile << numOfPairs << endl;
    
    // some processing to figure out when to print a 'progress signal'
    float rate = 0.1;
    int progress = 0;
    int control = (int) (rate * numOfPairs);                // show progress after every 'control' prints
    
    if (verboseLevel(VERBOSE_LOW)) {
        cout << "Generating f polynomials..." << endl;
    }

    string poly;
    // generate the f polynomials   
    for (int i = numOfPairs-1; i >= 0; i--) {
        // get the polynomial of degree i
        generateSinglePolynomial(i, poly);
        
        
        // print to output file
        outputFile << "f_" << i << " = " << poly.c_str() << endl;
        
        // print progress
        if (verboseLevel(VERBOSE_LOW)) {
            if (i % control == 0 && i > 0) {
                progress += rate * 100;
                cout << "\t" << progress << "\% Done" << endl;
            }
        }
    }
    
    if (verboseLevel(VERBOSE_LOW)) {
        cout << "Generating g polynomials..." << endl;
    }
    
    progress = 0;
    
    // generate the g polynomials
    for (int i = numOfPairs-1; i >= 0; i--) {
        // get the polynomial of degree i
        generateSinglePolynomial(i, poly);
        
        // print to output file
        outputFile << "g_" << i << " = " << poly.c_str() << endl;
        
        // print progress
        if (verboseLevel(VERBOSE_LOW)) {
            if (i % control == 0 && i > 0) {
                progress += rate * 100;
                cout << "\t" << progress << "\% Done" << endl;
            }
        }
    }
    
    outputFile.close();
}




void generateSinglePolynomial (int degree, string &polyString) {
    /* Code Sketch:
     *  1) Given the sequence of numbers from 0 to 'degree', shuffle
     *      them around
     *  2) Calculate the number of monomials with nonzero coefficients:
     *      numOfNonzeroCoef = degree - (s * degree)
     *  3) Take the first 'numOfNonzeroCoef' elements from the shuffled list;
     *      these will be the degrees of the monomials with nonzero coefficients
     *  4) Sort the retrieved list
     *  5) Generate random coefficients for these monomials, and place them
     *      accordingly in the polynomial string
     */
	
    srand ( unsigned ( time (NULL) ) );
    vector<int> myvector;
    vector<int>::iterator it;
    
    /* STEP (1): Create a shuffled list of numbers */
    
    // set ordered values
    for (int i=0; i<=degree; i++) myvector.push_back(i); // 0 1 2 3 4 5 ...
    
    // shuffle the list of numbers
    random_shuffle ( myvector.begin(), myvector.end() );
    
    
    /* STEP (2): Calculate number of monomials with nonzero coefs */
    unsigned int numOfNonzeroCoef = degree+1 - (int) (sparsity * degree);
    
    // if degree = 0, we want a one monomial with degree = 0, coef != 0
    if (degree == 0)
        numOfNonzeroCoef = 1;
    
    
    /* STEP (3): Get the degree of the monomials with nonzero coefs */
    vector<int> validDegrees (numOfNonzeroCoef);
    
    for (unsigned int i = 0; i < numOfNonzeroCoef; i++) {
        validDegrees[i] = myvector[i];
    }

    
    /* STEP (4): Sort the list of valid degrees */
    sort(validDegrees.begin(), validDegrees.end());
    reverse(validDegrees.begin(), validDegrees.end());
    
    /* STEP (5): Create the polynomial string */
    polyString = "";
    stringstream ss;
    
    int coef = 0;
    string temp = "";
    for (size_t i = 0; i < validDegrees.size(); i++) {
        // generate the coefficient
        coef = (int) (rand() % (COEF_LIMIT-1)) + 1;
        
        ss << coef << "*x^" << validDegrees[i];
        
        ss >> temp;
        polyString.append(temp);
        
        if (i < validDegrees.size()-1) {
            polyString.append(" + ");
        }
        
        temp = "";
        ss.clear();
    }   
    
    

}


/*
 *	Generate polynomials based on parameters S and N:
 *	Based on sparsity S and number N (where N is the number of polynomials pairs),
 *	generate a two sets of polynomials, where each set is defined as follows:
 *		- The first polynomial of the set has sparsity S and N monomials (its degree is then calculated)
 *		- The second polynomial has sparsity S and N-1 monomials
 *		- The third polynomial has sparsity S and N-2 monomials
 *		- ....
 *		- The last polynomial has sparsity S and 1 monomial
 *		
 */
// void generateSparsePolynomials (float sparsity, size_t numberOfPairs) {
// 	// set the output stream
//     ofstream outputFile;
//     outputFile.open (outputFilename.c_str(), ios::out);
// 
//     // write the number of pairs to the output file
//     outputFile << numberOfPairs << endl;
// 
//     // some processing to figure out when to print a 'progress signal'
//     float rate = 0.1;
//     int progress = 0;
//     int control = (int) (rate * numOfPairs);                // show progress after every 'control' prints
// 
//     if (verboseLevel(VERBOSE_LOW)) {
//         cout << "Generating f polynomials..." << endl;
//     }
//     string poly;
//     // generate the f polynomials   
//     for (int i = numOfPairs-1; i >= 0; i--) {
//         // get the polynomial of degree i
//         generateSinglePolynomial(i, poly);
// 
// 
//         // print to output file
//         outputFile << "f_" << i << " = " << poly.c_str() << endl;
// 
//         // print progress
//         if (verboseLevel(VERBOSE_LOW)) {
//             if (i % control == 0 && i > 0) {
//                 progress += rate * 100;
//                 cout << "\t" << progress << "\% Done" << endl;
//             }
//         }
//     }
// 
//     if (verboseLevel(VERBOSE_LOW)) {
//         cout << "Generating g polynomials..." << endl;
//     }
// 
//     progress = 0;
// 
//     // generate the g polynomials
//     for (int i = numOfPairs-1; i >= 0; i--) {
//         // get the polynomial of degree i
//         generateSinglePolynomial(i, poly);
// 
//         // print to output file
//         outputFile << "g_" << i << " = " << poly.c_str() << endl;
// 
//         // print progress
//         if (verboseLevel(VERBOSE_LOW)) {
//             if (i % control == 0 && i > 0) {
//                 progress += rate * 100;
//                 cout << "\t" << progress << "\% Done" << endl;
//             }
//         }
//     }
// 
//     outputFile.close();
// }



/*
 *      Given the maximum size of the queue we wish to fill, this method computes the number
 *      of polynomial pairs (as generated by testCase1) needed to try and acheive this size
 *      (up to a certain error)
 */
int computeNumberOfPairsNeeded1(size_t maxSizeOfQueue) {
    size_t currentSize = 0;
    
    int numberOfPairs = 0;
    unsigned long monom_i = 1;              // represents the number of monomials in f_i
    // keep looping and stop just before currentSize exceeds maxSizeOfQueue
    while (true) {
        cout << "Current Size = " << currentSize << ", |f_i| = " << monom_i << endl;
        // add |f_i| elements (ie the number of monomials in f_i)
        currentSize += monom_i;
        
        // increment the degree and number of pairs
        monom_i++;
        numberOfPairs++;
        
        // check if the next iteration will lead to currentSize > maxSizeOfQueue
        if (currentSize + monom_i > maxSizeOfQueue) {
            break;
        }
    }
    
    // STATE: the current number of pairs will yield a queue size <= maxSizeOfQueue
    
    // the following section of the code will allow the queue size to exceed the 
    // mentioned maxSizeOfQueue, by a certain percentage offset
    
    // Ex: if offset = 0.3, then the code will try (if possible) to increase currentSize to
    // maxSizeOfQueue + maxSizeOfQueue*offset (where 0 < offset < 1)
    float offset = 0.3;
    
    cout << endl << "Testing offset (" << offset << ")" << endl;
    cout  << "Current Size = " << currentSize << ", monom_i = " << monom_i << ", pairs = " << numberOfPairs <<  endl;
    
    cout << endl << "currentSize + monom_i <= maxSizeOfQueue * (1 + offset):" << endl;
    cout << currentSize + monom_i << " <= " << maxSizeOfQueue * (1 + offset) << endl;
    if (currentSize + monom_i <= maxSizeOfQueue * (1 + offset)) {
        currentSize += monom_i;
        numberOfPairs++;
    }
    
    return numberOfPairs - 1;                       // -1   FOR NOW
}


// 
// void readArgs (int argc, char *argv[]) {
// 	if (argc == 1)
// 		usage();
//     
//     // read any input parameters that the user has given
//     for (int i = 1; i < argc; i++) {
//         if (strcmp(argv[i], "-o") == 0) {               // output file name
//             outputFilename = argv[++i];
//         }
//         else if (strcmp(argv[i], "-n") == 0) {  // number of pairs
//             numOfPairs = atoi(argv[++i]);
//         }
//         else if (strcmp(argv[i], "-v") == 0) {  // enable verboseLevel(VERBOSE_LOW) feature
//             verboseLevel(VERBOSE_LOW) = true;
//         }
//         else if (strcmp(argv[i], "-s") == 0) {  // set sparsity
//             sparsity = atof(argv[++i]);
//         }
// 		else if (strcmp(argv[i], "-help") == 0) {	// print usage
// 			usage();
// 		}
//     }
// 
// 	// if maxQueueSize is found, override the number of pairs
// 	for (int i = 1; i < argc; i++) {
// 		if (strcmp(argv[i], "-qs") == 0) {      // maximum queue size we want (with error)
// 			int qs = atoi(argv[++i]);
// 
// 			// get the number of pairs from this queue size
// 			numOfPairs = computeNumberOfPairsNeeded1 (qs);
// 		}
// 	}
// }
// 
// 
// // print usage message (to standard error)
// void usage () {
//     cerr << "usage: Polynomial Generator" << endl <<endl;
//     
//     cerr << "generator  -o (output_name) -n (num) -s (num)" << endl;
//     cerr << "Where: -o (output_name) = name of output file (default = \"" << DEFAULT_OUTPUT_FILE_NAME.c_str() << "\")" << endl;
//     cerr << "       -n (num) = number of pairs (default = " << DEFAULT_NUMBER_OF_PAIRS << ")" << endl;
//     cerr << "       -s (num) = sparsity of polynomials, 0 <= s <= 1 (default = " << DEFAULT_SPARSITY << ")" << endl;
//     
//     cerr << endl;
//     
//     // quit program
//     exit(0);
// }
