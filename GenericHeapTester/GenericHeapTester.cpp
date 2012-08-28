//
//  GenericHeapTester.cpp
//  PAL v1.0
//
//  Created by Karl Gemayel on 07/28/12.
//  Copyright (c) 2011 American University of Beirut. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <string>
#include <sstream>
#include <time.h>
#include "../Options/Options.h"
#include "../Heap/Heap.h"
#include "../FunnelHeapWithMerging/FunnelHeapWithMerging.h"
#include "../FunnelHeap/FunnelHeap.h"
#include "../BinaryHeap/BinaryHeap.h"
#include "../BinaryHeapWithChaining/BinaryHeapWithChaining.h"

#include "../Polynomial/Monomial.h"





using namespace std;

void usage();
void readArgs(int argc, char *argv[]);
void setupHeap(Heap *&A);

// test functions
void test1(size_t s, size_t p, Heap *A);
void test2(size_t s, size_t p, Heap *A);
void test3(size_t s, size_t p, Heap *A);


/****************************************/

// arguments from main
bool checkHeap = false;		// enable to check correctness of the heap
bool working = true;

size_t s_global;		// global variable holding the 's' parameter for tests
size_t p_global;		// global variable holding the 'p' parameter for tests
int test_id;

// main function
int main (int argc, char * argv[]) {

	// read/update arguments
	readArgs(argc, argv);
    
	Heap *A;
	
	setupHeap(A);
	
	startTimer();
	// run test
	if (test_id == 1) 
		test1(s_global, p_global, A);
	else if (test_id == 2) 
		test2(s_global, p_global, A);
	else if (test_id == 3)
		test3(s_global, p_global, A);
	
	stopTimer();
    

	free(A);		// free the heap
	
	// if 'checkHeap' is set to true, print whether or not the heap works properly
	if (checkHeap) {
		if (working) cout << "Heap is working properly" << endl;
		else printWarning("Heap did not work properly");
	}

	// compute elapsed time
	cout << "Time = " << getElapsedTime() << endl;   
	
	// print statistics
	cout << *(stxxl::stats::get_instance()) << endl;

	// open a file to write the execution time to
	ofstream timesFile;
	string timeFilename;
	getTimeFilename(timeFilename);
	timesFile.open (timeFilename.c_str(), ios::out | ios::app);
	
	// write the input size with the time

	
	// close the times file
	timesFile.close();
	

	
}

// read/update the arguments for running the test case
void readArgs(int argc, char *argv[]) {
	// if there are no arguments (except for the program name), then print usage
	if (argc == 1) 
		usage();	
		
	// variables that will be used to check if all parameters are valid
	bool choseHeap = false;
	bool choseTestID = false;
	bool choseS = false;
	bool choseP = false;
	
	// if there are arguments, read them
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-t") == 0) {
			test_id = atoi(argv[++i]);
			
			// allow only one test id to be chosen
			if (!choseTestID) choseTestID = true;
			else {
				printError("Chose more than one test id.\n");
				exit(1);
			}
		}
		else if (strcmp(argv[i], "-s") == 0) {
			s_global = atoll(argv[++i]);
			
			choseS = true;
		}
		else if (strcmp(argv[i], "-p") == 0) {
			p_global = atoll(argv[++i]);
			
			choseP = true;
		}
		else if (strcmp(argv[i], "-v") == 0) {
			int vLevel = atoi(argv[++i]);
			setVerbose(vLevel);
		}
		else if (strcmp(argv[i], "-h") == 0) {			// choose heap
			chooseHeap(argv[++i]);	
			
			// allow only one heap to be chosen
			if (!choseHeap) choseHeap = true;
			else {
				printError("Only one heap type can be selected.\n");
				exit(1);
			}
		}
		else if (strcmp(argv[i], "-check") == 0) {
			checkHeap = true;
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
	if (!choseHeap) printError("Please choose a heap type (see -h).");
	else if (!choseS) printError("Please set a test size (see -s).");
	else if (!choseP) printError("Please set test parameter (see -p).");
	else if (!choseTestID) printError("Please choose a test id (see -t).");
		
}

void usage() {
	cerr << "usage: Generic Heap Tester" << endl <<endl;

	cerr << "Driver: GenericHeaptester -h (heap_name) -t (test_id) -s (test_size) -p (test_param)" << endl;
	cerr << "Where: -h (heap_name) = name of heap used; names below (default binary)" << endl;
	cerr << "       -t (test_id) = the type of test (see below)" << endl;
	cerr << "		-s (test_size) = the size of the test" << endl;
	cerr << "		-p (test_param) = parameters for the chosen test" << endl;
	cerr << endl;
	
	cerr << "Heap Names:" << endl;
	cerr << "\tbinary --> standard binary heap" << endl;
	cerr << "\tbinary-chaining --> binary heap (with chaining)" << endl;
	cerr << "\tfunnel --> funnel heap" << endl;
	cerr << "\tfunnel-merging --> funnel heap (with merging optimization)" << endl;
	cerr << endl;
	
	cerr << "Test Types:" << endl;
	cerr << "\t1 --> Add 's' elements. Perform 'p' sequence of 's'/2 deletes and 's'/2 inserts" << endl;
	cerr << "\t		 consecutively. Finally perform 's' deletes, leaving the heap empty." << endl;
	cerr << "\t2 --> Add 's' elements. Then perform 'p' consecutive delete/insert operations," << endl;
	cerr << "\t		 and finally perform 's' deletes, leaving the heap empty." << endl;
	cerr << endl;
	
	cerr << "Optional Arguments:" << endl;
	cerr << "\t-v (num): puts the program in verbose mode (prints information)" << endl;
	cerr << "\t\tnum: the level of verbosity (1, 2, or 3)" << endl;
	cerr << "\t-check: checks that the heap is running correctly. NOTE: if this option is enabled," << endl;
	cerr << "\t		   the runtime displayed at the end will not be accurate for profiling reasons." << endl;

	// quit program
	exit(0);
}




// TESTS

/*
 * test1: Add 's' elements. Perform 'p' sequence of 's'/2 deletes and 's'/2 inserts 
 *		  consecutively. Finally perform 's' deletes, leaving the heap empty.
 */
void test1 (size_t s, size_t p, Heap *A) {
	
	if (verboseLevel(VERBOSE_LOW)) {
		cout << "---------Performing Test1----------" << endl;
		cout << "---------Heap Size = " << s << " - NOI = " << p << "------" << endl;
	}
	
	
	// FOR NOW: p is the passed as the total number of inserts
	// Calculate it again so that it represents the total number of times
	// to perform the N/2 inserts/deletes
	size_t numOfTimes;
	numOfTimes = p - s;  // remove the number of inserts done at the first
 	numOfTimes = numOfTimes / (s/2);

	p = numOfTimes;
	
	
	// setup the random generator
	const size_t UPPER_LIMIT = 1000;
	srand((size_t)time(NULL));
	size_t num = 0;
	
	// setup stxxl stats
	stxxl::stats *stats = stxxl::stats::get_instance();
	stats->reset();	
	
	
	// queue for checking heap correctness (if enabled)
	priority_queue<size_t> checker;
	
	// start by adding 's' random elements
	for (size_t i = 0; i < s; i++) {
		num = (size_t) (rand() % UPPER_LIMIT) + 1;       // generate a number to insert
		
		A->insert(num, 1, 0);
		
		// insert in checker
		if (checkHeap) checker.push(num);
	}
	
	// NOTE: In the case of the binary heap, large cases can take months (when used with STXXL).
	// Therefore, a mechanism has been implemented (only to be used for the binary heap) that stops
	// execution after a (pre-defined) number of iterations, and the total runtime will be estimated
	// from the information gathered during the executed iterations.
	unsigned int printCounter = 0;
	const unsigned int PRINT_COUNTER_MAX = 10;
	time_t counter_start = time(NULL), counter_end;
	int numOfPrintsDone = 0;	// initialize the number of prints already done to 0
	const int MAX_NUM_OF_PRINTS = 5;// the maximum number of prints to be done
	
	// the following variables will be used to calculate the average statistics of the iterations done
	// when using the binary heap, in order to better predict the total outcome
	double timeAcc = 0;
	double IOWaitAcc = 0;
	size_t readsAcc = 0, writesAcc = 0;
	size_t diskReadsAcc = 0, diskWritesAcc = 0;
	
	// now perform 'p' sequence of s/2 deletes and s/2 inserts consecutively
	for (size_t p_counter = 0; p_counter < p; p_counter++) {
		
		if (chosenHeap() == BINARY_HEAP && printCounter == PRINT_COUNTER_MAX) {
			// stop timer for progress 
			counter_end = time(NULL);
			
			// calculate the begining/end iterations that are being profiled
			size_t begIterationRange = PRINT_COUNTER_MAX * numOfPrintsDone;
			size_t endIterationRange = begIterationRange + PRINT_COUNTER_MAX;
			cout << "\tIteration #: " << begIterationRange << " - " << endIterationRange << ", ";
			cout << "Time = " << difftime(counter_end, counter_start) << ", ";
			cout << "I/O Wait Time = " << stats->get_io_wait_time() << ", ";
			cout << "Reads = " << stats->get_reads() << ", ";
			cout << "Writes = " << stats->get_writes() << ", ";
			cout << "Bytes Read From Disk = " << stats->get_read_volume() << "(" <<  stats->get_read_volume()/(1024*1024) << " Mib), ";
			cout << "Bytes Written to Disk = " << stats->get_written_volume() << "(" << stats->get_written_volume()/(1024*1024) << " Mib)" << endl;
			
			// increment accumulators
			timeAcc += difftime(counter_end, counter_start);
			IOWaitAcc += stats->get_io_wait_time();
			readsAcc += stats->get_reads();
			writesAcc += stats->get_writes();
			diskReadsAcc += stats->get_read_volume();
			diskWritesAcc += stats->get_written_volume();
			
			
			stats->reset();
			printCounter = 0;	// reset the counter

			counter_start = time(NULL);			// restart the timer 
			

			numOfPrintsDone++;
			if (numOfPrintsDone == MAX_NUM_OF_PRINTS) {
				// calculate the averages of previous iterations
				double timeAvg = timeAcc/numOfPrintsDone;
				double IOWaitAvg = IOWaitAcc/numOfPrintsDone;
				size_t readsAvg = readsAcc/numOfPrintsDone, writesAvg = writesAcc/numOfPrintsDone;
				size_t diskReadsAvg = diskReadsAcc/numOfPrintsDone, diskWritesAvg = diskWritesAcc/numOfPrintsDone;
				
				cout << "Averages: " << endl;
				cout << "Time = " << timeAvg << ", ";
				cout << "I/O Wait Time = " << IOWaitAvg << ", ";
				cout << "Reads = " << readsAvg << ", ";
				cout << "Writes = " << writesAvg << ", ";
				cout << "Bytes Read From Disk = " << diskReadsAvg << "(" <<  diskReadsAvg/(1024*1024.0) << " Mib), ";
				cout << "Bytes Written to Disk = " << diskWritesAvg << "(" << diskWritesAvg/(1024*1024.0) << " Mib)" << endl;
				
				cout << endl;
				
				// predict total values
				double timeTotal = (timeAvg * (p/(double)PRINT_COUNTER_MAX)) /(3600*24.0);
				double IOWaitTotal = IOWaitAvg * (p/(double)PRINT_COUNTER_MAX);
				size_t readsTotal = readsAvg * (p/(double)PRINT_COUNTER_MAX), writesTotal = writesAvg * (p/(double)PRINT_COUNTER_MAX);
				size_t diskReadsTotal = diskReadsAvg * (p/(double)PRINT_COUNTER_MAX), diskWritesTotal = diskWritesAvg * (p/(double)PRINT_COUNTER_MAX);
				
				cout << "Total Prediction: " << endl;
				cout << "Time = " << timeTotal << " days, ";
				cout << "I/O Wait Time = " << IOWaitTotal << " (s), ";
				cout << "Reads = " << readsTotal << ", ";
				cout << "Writes = " << writesTotal << ", ";
				cout << "Bytes Read From Disk = " << diskReadsTotal << "(" <<  diskReadsTotal/(1024*1024.0) << " Mib), ";
				cout << "Bytes Written to Disk = " << diskWritesTotal << "(" << diskWritesTotal/(1024*1024.0) << " Mib)" << endl << endl;
				exit(0);
			}
			
		}
		
		printCounter++;
		
		size_t deletedElem = 0;
		// perform s/2 deletes
		for (size_t i = 0; i < (size_t) s/2; i++) {
			deletedElem = A->poll();
			
			if (checkHeap) {
				if (GET_DEGREE(deletedElem) != checker.top()) {
					working = false;
				}
				
				checker.pop();
			}
		}
		
		
		// perform s/2 inserts
		for (size_t i = 0; i < (size_t) s/2; i++) {
			num = (size_t) (rand() % UPPER_LIMIT) + 1;       // generate a number to insert
			
			A->insert(num, 1, 0);
			
			if (checkHeap) checker.push(num);
		}
	}
	
	size_t deletedElem = 0;
	// finally perform s deletes
	for (size_t i = 0; i < s; i++) {
		deletedElem = A->poll();
		
		if (checkHeap) {
			if (GET_DEGREE(deletedElem) != checker.top()) {
				working = false;
			}
			
			checker.pop();
		}
	}
}


/*
 * test2: Add 'n' elements. Then perform 'p' consecutive delete/insert operations, 
 *		  and finally perform 'n' deletes, leaving the heap empty.
 */
void test2 (size_t s, size_t p, Heap *A) {
	// setup the random generator
	const size_t UPPER_LIMIT = 1000;
	srand((size_t)time(NULL));
	size_t num = 0;
	
	// queue for checking heap correctness (if enabled)
	priority_queue<size_t> checker;
	
	
	// start by adding s elements
	for (size_t i = 0; i < s; i++) {
		num = (size_t) (rand() % UPPER_LIMIT) + 1;       // generate a number to insert
		
		A->insert(num, 1, 0);
		
		// insert in checker
		if (checkHeap) checker.push(num);
	}
	
	
	// perform p consecutive delete/insert operations
	size_t deletedElem = 0;
	
	for (size_t p_counter = 0; p_counter < p; p_counter++) {
		
		// perform delete 
		deletedElem = A->poll();
		
		if (checkHeap) {
			if (deletedElem != checker.top())
				working = false;
			
			checker.pop();
		}
		
		// perform insert
		num = (size_t) (rand() % UPPER_LIMIT) + 1;       // generate a number to insert
		
		A->insert(num, 1, 0);
		
		if (checkHeap) checker.push(num);
	}
	
	
	
	// perform s deletes
	for (size_t i = 0; i < s; i++) {
		deletedElem = A->poll();
		
		if (checkHeap) {
			if (deletedElem != checker.top())
				working = false;
			
			checker.pop();
		}
	}
}



/*
 * test3: SHOULD ONLY BE USED FOR FUNNEL HEAP WITH MERGING:
 *		  
 */
void test3 (size_t s, size_t p, Heap *A) {
	
	if (verboseLevel(VERBOSE_LOW)) {
		cout << "---------Performing Test 3----------" << endl;
		cout << "---------Heap Size = " << s << " - NOI = " << p << "------" << endl;
	}
	
	// setup the random generator
	const size_t UPPER_LIMIT = 100;
	srand((size_t)time(NULL));
	size_t num = 0;
	
	// queue for checking heap correctness (if enabled)
	priority_queue<size_t> checker;
	
	
	for (size_t i = 0; i < s; i++) {
		
		num = (size_t) (rand() % UPPER_LIMIT) + 1;       // generate a number to insert
		
		A->insert(num, 1, 0);
		

		if (checkHeap) checker.push(num);
	}
	
	size_t deletedElem = 0;
	
	
	while (!A->isEmpty()) {
		deletedElem = A->poll();
		
		if (checkHeap) {
			if (GET_DEGREE(deletedElem) != checker.top())
				working = false;
			
			checker.pop();
			while (!checker.empty() && GET_DEGREE(deletedElem) == checker.top())
				checker.pop();
		}
	}
	
	
	
}




void setupHeap(Heap *&A) {
	A = new FunnelHeapWithMerging (6);
	if (chosenHeap() == BINARY_HEAP_WITH_CHAINING) {	
		// for the binaryheap with chaining, we first need to perform some preprocessing
		// to determine some characteristics of the chains. Thie preprocessing should not
		// be included in the profiling sections.

		printError("Generic testing for Binary Heap with Chaining has not yet been implemented!\n Please choose another heap.");
		exit(0);

		if (verboseLevel(VERBOSE_LOW)) cout << "Using Binary Heap (with chaining)"  << endl;
		// allocate heap
		A = new BinaryHeapWithChaining (s_global);

	}
	else if (chosenHeap() == FUNNEL_HEAP_WITH_MERGING) {
	
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Funnel Heap (with merging)" << endl;
		
		A = new FunnelHeapWithMerging (6);	
	}
	else if (chosenHeap() == BINARY_HEAP) {

		if (verboseLevel(VERBOSE_LOW)) cout << "Using Binary Heap (capacity = " << s_global << ")" << endl;
		
		A = new BinaryHeap (s_global);		
	}
	else if (chosenHeap() == FUNNEL_HEAP) {
		
		if (verboseLevel(VERBOSE_LOW)) cout << "Using Funnel Heap" << endl;
		
		A = new FunnelHeap (6);
	}
	else {
		cout << "ERROR: No heap was chosen" << endl;
		usage();
	}
}