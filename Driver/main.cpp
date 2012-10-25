//
//  main.cpp
//  PAL v1.0
//
//  Created by Karl Gemayel on 12/3/11.
//  Copyright (c) 2011 American University of Beirut. All rights reserved.
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
#include "../PolynomialArithmetic/PolynomialArithmetic.h"
#include "../PolynomialGenerator/PolynomialGenerator.h"
#include "../Options/Options.h"
#include "../IOHandling/IOHandling.h"

using namespace std;

void usage();
void readArgs(int argc, char *argv[]);

/****************************************/

// arguments from main
string inputName = "mult.in";
string programType = "";
bool writeMultiplicationOutputToFile = false;

bool type_gen = false;		// polynomial generator
bool type_sop = false;		// summation of products

unsigned int GEN_numOfPairs;
float GEN_sparsity;

int main (int argc, char * argv[]) {

	// read/update arguments
	readArgs(argc, argv);
	
	
	if (type_sop) {	
		int numOfPairs = 0;
	    vector<poly_t> f_polynomials;
	    vector<poly_t> g_polynomials;
    

		readInputFile(inputName, f_polynomials, g_polynomials, numOfPairs);
    
    
		poly_t result;
	
		startTimer();
		summationOfProducts(f_polynomials, g_polynomials, result);
		stopTimer();
    
	
		if (verboseLevel(VERBOSE_LOW)) cout << "result size = " << result.size() << endl;

		// compute elapsed time
		if (verboseLevel(VERBOSE_LOW)) cout << "Time = " << getElapsedTime() << endl;   

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
	
		if (verboseLevel(VERBOSE_LOW)) cout << *(stxxl::stats::get_instance());
	}
	else if (type_gen) {
		string outputFilename;
		getOutputFilename(outputFilename);
		generatePolynomials (GEN_numOfPairs, GEN_sparsity, outputFilename);
	}
}

void usage() {
	cerr << "Usage of the Polynomial Arithmetic Library" << endl <<endl;

	cerr << "Summation of Products: Computes the sum of polynomial pairwise-products" << endl;
	cerr << "usage: main -t sop -h heap_name [-osi] -i input_file [-o output_file] [-v verbose_level] " << endl;
	cerr << "Where: -h heap_name: name of heap used; names below (default binary)" << endl;
	cerr << "       -osi: enables the optimized sequence of inserts" << endl;
	cerr << "       -i input_file: name of input file" << endl;
	cerr << "       -o output_file: name of the output file" << endl;
	cerr << "       -v verbose_level: the degree of the verbose mode (1 - 3)" << endl;
	cerr << endl << endl;
	
	cerr << "Polynomial Generator: Generates polynomials of the summation of products routine" << endl;
	cerr << "usage: main -t gen -n npairs [-s sparsity] -o output_file" << endl;
	cerr << "Where: -n npairs: the number of polynomial pairs to generate" << endl;
	cerr << "       -s sparsity: the sparsity of the polynomials" << endl;
	cerr << "       -o output_file: name of the output file" << endl;
	cerr << endl << endl;
	
	cerr << "Heap Names:" << endl;
	cerr << "\tbinary --> standard binary heap" << endl;
	cerr << "\tbinary-chaining --> binary heap (with chaining)" << endl;
	cerr << "\tfunnel --> funnel heap" << endl;
	cerr << "\tfunnel-merging --> funnel heap (with merging optimization)" << endl;
	cerr << endl;

	// quit program
	exit(0);
}


// read/update the arguments for running the test case
void readArgs(int argc, char *argv[]) {
	// if there are no arguments (except for the program name), then print usage
	if (argc == 1) 
		usage();	
		
	// booleans to check validity of (number of) arguments
	bool choseHeap = false;
	bool choseInputFile = false;
	bool choseOutputFile = false;
	bool choseProgramType = false;
	bool choseNumberOfPairs = false;
	
	
	// if there are arguments, read them
	for (int i = 1; i < argc; i++) {
		if (strcasecmp(argv[i], "-t") == 0) {			// program type
			if (choseProgramType) {
				printError("Only one program type can be used.");
				exit(1);
			}
			programType = argv[++i];
			choseProgramType = true;			
		}
		else if (strcasecmp(argv[i], "-i") == 0) {		// input file
			if (choseInputFile) {
				printError("Only one input file can be used.");
				exit(1);
			}
			inputName = argv[++i];
			choseInputFile = true;
		}
		else if (strcasecmp(argv[i], "-v") == 0) {		// verbose level
			int vLevel = atoi(argv[++i]);
			setVerbose(vLevel);
		}
		else if (strcasecmp(argv[i], "-o") == 0) {		// output file
			printMultOutput(true);
			setOutputFilename(argv[++i]);
			choseOutputFile = true;
		}
		else if (strcasecmp(argv[i], "-h") == 0) {		// heap name
			if (choseHeap) {
				printError("Only one heap name can be used.");
				exit(1);
			}
			chooseHeap(argv[++i]);
			choseHeap = true;
		}
		else if (strcasecmp(argv[i], "-osi") == 0) {	// osi
			setOptimizedSequenceOfInserts(true);
		}
		
		// arguments for the generator type
        else if (strcasecmp(argv[i], "-n") == 0) {  // number of pairs
            GEN_numOfPairs = atoi(argv[++i]);
			choseNumberOfPairs = true;
        }
        else if (strcasecmp(argv[i], "-s") == 0) {  // set sparsity
            GEN_sparsity = atof(argv[++i]);
        }

		// help and invalid arguments
		else if (strcasecmp(argv[i], "-help") == 0) {	// help
			usage();
		}
		else {
			string invArg = argv[i];
			printError("Invalid argument \"" + invArg + "\"");
			exit(1);
		}
	}
	
	
	// check if any "required" parameter was not given, then print message
	if (!choseProgramType) {
		printError("Please choose a program type (ex sop, gen)");
		exit(1);
	}
	
	// check parameters for sop (i.e. summation of products)
	if (strcasecmp(programType.c_str(), "sop") == 0) {
		type_sop = true;		
		if (!choseHeap) {
			printError("Please choose a heap name.");
			exit(1);
		}
		if (!choseInputFile) {
			printError("Please choose an input filename.");
			exit(1);
		}
	}
	
	// check parameters for gen (i.e. polynomial generator)
	if (strcasecmp(programType.c_str(), "gen") == 0) {
		type_gen = true;
		if (!choseNumberOfPairs) {
			printError("Please choose the number of pairs.");
			exit(1);
		}
		if (!choseOutputFile) {
			printError("Please choose an output filename.");
			exit(1);
		}
	}
		
}