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
#include "../Options/Options.h"

using namespace std;

void usage();
void readArgs(int argc, char *argv[]);

/****************************************/

// arguments from main
string inputName = "mult.in";
bool writeMultiplicationOutputToFile = false;

// global variables for manipulating user-defined options (ex: verbose mode)
//extern bool VERBOSE = false;

int main (int argc, char * argv[]) {

	// read/update arguments
	readArgs(argc, argv);
    
    ifstream in(inputName.c_str());
    string line;
    
	
    getline(in, line);
    int numOfPairs = atoi(line.c_str());
    
    vector<poly_t> f_polynomials (numOfPairs);
    vector<poly_t> g_polynomials (numOfPairs);
    
    poly_t f, g;
	
    // loop over input file and collect all pairs

	// read the f_polynomials
    for (int i = 0; i < numOfPairs; i++) {
        
        // read the f polynomial of the i'th pair
        getline(in, line);
        line = line.substr(line.find("=") + 1);
        parsePolynomialString(line, f);
                
        
        // add f to its corresponding vectors
        f_polynomials[i] = f;

		string f_string;
		toString(f, f_string);
		
		// reset the polynomials to read new ones
		f.clear();
    }

	// read the g_polynomials
    for (int i = 0; i < numOfPairs; i++) {
	
        // read the g polynomial of the i'th pair
        getline(in, line);
        line = line.substr(line.find("=") + 1);
        parsePolynomialString(line, g);
        
        
        // add g to its corresponding vectors
        g_polynomials[i] = g;

		string g_string;
		toString(g, g_string);
		
		// reset the polynomial to read a new one
		g.clear();
    }

	poly_t result;
//	multiplyMultiplePairs(f_polynomials, g_polynomials, result);
//	multiplyMultiplePairsOptimized(f_polynomials, g_polynomials, result);
//	OSImultiply(f_polynomials, g_polynomials, result);
	
	startTimer();
	summationOfProducts(f_polynomials, g_polynomials, result);
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
	    toString(result, result_string);
	
		string outputFilename;
		getOutputFilename(outputFilename);
		
		ofstream mult_output;
		mult_output.open (outputFilename.c_str(), ios::out | ios::app);
		mult_output << result_string.c_str() << endl;
		mult_output.close();
	}
	
	cout << *(stxxl::stats::get_instance());
}

// read/update the arguments for running the test case
void readArgs(int argc, char *argv[]) {
	// if there are no arguments (except for the program name), then print usage
	if (argc == 1) 
		usage();	
		
	bool choseOneHeapType = false;
	int numOfHeapsChosen = 0;
	bool choseInputFile = false;
	
	// if there are arguments, read them
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-i") == 0) {	// choosing input file
			inputName = argv[++i];
			choseInputFile = !choseInputFile;
		}
		else if (strcmp(argv[i], "-v") == 0) {
			int vLevel = atoi(argv[++i]);
			setVerbose(vLevel);
		}
		else if (strcmp(argv[i], "-o") == 0) {
			printMultOutput(true);
			setOutputFilename(argv[++i]);
		}
		else if (strcmp(argv[i], "-h") == 0) {
			chooseHeap(argv[++i]);
			
			numOfHeapsChosen++;
		}
		else if (strcmp(argv[i], "-osi") == 0) {
			setOptimizedSequenceOfInserts(true);
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
	choseOneHeapType = numOfHeapsChosen == 1;
	if (!(choseInputFile) || !choseOneHeapType) 
		usage();
		
}

void usage() {
	cerr << "usage: Polynomial Arithmetic" << endl <<endl;

	cerr << "Driver: main -h (heap_name) -i (input_file)" << endl;
	cerr << "Where: -h (heap_name) = name of heap used; names below (default binary)" << endl;
	cerr << "       -i (input_file) = name of input file" << endl;
	cerr << endl;
	
	cerr << "Heap Names:" << endl;
	cerr << "\tbinary --> standard binary heap" << endl;
	cerr << "\tbinary-chaining --> binary heap (with chaining)" << endl;
	cerr << "\tfunnel --> funnel heap" << endl;
	cerr << "\tfunnel-merging --> funnel heap (with merging optimization)" << endl;
	cerr << endl;
	
	cerr << "Optional Arguments:" << endl;
	cerr << "\t-osi: performs the multiplication using an optimized sequence of inserts" << endl;
	cerr << "\t-v (num): puts the program in verbose mode (prints information)" << endl;
	cerr << "\t\tnum: the level of verbosity (1, 2, or 3)" << endl;
	cerr << "\t-o (output_file): print the output to a file with name 'output_file'" << endl;

	// quit program
	exit(0);
}