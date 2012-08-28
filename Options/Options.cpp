// Options.h

/* 
 *	This file defines options and modes (configurable by the user) that
 *	can be executed at runtime. These include: verbose, debugging, etc..
 * 	It also includes profiling tools: timer.
 */


#include <iostream>
#include <string>
#include <time.h>
#include "Options.h"
// #include <stxxl.h>

using namespace std;

// Define error and warning printing methods
void printError(string err) {
	cerr << "Error: " << err.c_str() << endl;
}

void printWarning(string warning) {
	cerr << "Warning: " << warning.c_str() << endl;
}


// HANDLING OUTPUT FILES (time and multiplication output)
string time_filename = "";
void getTimeFilename(string &filename) {
	filename = time_filename;
}

void setTimeFilename(const string filename) {
	time_filename = filename;
}

string output_filename = "";
bool writeOutput = false;
bool printMultOutput() {
	return writeOutput;
}

void printMultOutput(bool value) {
	writeOutput = value;
}

void getOutputFilename(string &filename) {
	filename = output_filename;
}

void setOutputFilename(const string filename) {
	output_filename = filename;
}




// CHOOSING WHICH HEAP TO USE
unsigned int HEAP = 0;		// defaults to binary heap

// chooses the heap to be used in computation
void chooseHeap(string heap_name) {
	
	if (strcmp(heap_name.c_str(), "binary") == 0) {
		HEAP = BINARY_HEAP;
	}
	else if (strcmp(heap_name.c_str(), "funnel") == 0) {
		HEAP = FUNNEL_HEAP;
	}
	else if (strcmp(heap_name.c_str(), "funnel-merging") == 0) {
		HEAP = FUNNEL_HEAP_WITH_MERGING;
	}
	else if (strcmp(heap_name.c_str(), "binary-chaining") == 0) {
		HEAP = BINARY_HEAP_WITH_CHAINING;
	}
	else {
		printWarning ("Heap id invalid! Defaulting to binary heap.");
		HEAP = BINARY_HEAP;
	}
}

// returns the ID of the previously chosen heap
int chosenHeap () {
	return HEAP;
}



// USING THE PROFILER

time_t start, end;
double elapsedTime;
double totalIOWaitTime;
double totalReads;
double totalWrites;


void startTimer() {
	elapsedTime = 0;
	start = time(NULL);	// record the time that the task begins
	
	// // init stxxl stats
	// stats = stxxl::stats::get_instance();
	// stats->reset();
	// 
	// totalIOWaitTime = 0;
	// totalReads = 0;
	// totalWrites = 0;
}

void continueTimer() {
	start = time(NULL);
	
//	stats->reset();
}

void pauseTimer() {
	end = time(NULL);
	elapsedTime += difftime(end, start);
	
	// totalIOWaitTime += stats->get_io_wait_time();
	// totalReads += stats->get_reads();
	// totalWrites += stats->get_writes();
}

void stopTimer() {
	end = time(NULL);
	elapsedTime += difftime(end, start);
	
	// totalIOWaitTime += stats->get_io_wait_time();
	// totalReads += stats->get_reads();
	// totalWrites += stats->get_writes();
}

double getElapsedTime() {
	return elapsedTime;
}

double getIOWaitTime() {
	return totalIOWaitTime;
}

double getReads() {
	return totalReads;
}

double getWrites() {
	return totalWrites;
}

void resetStats() {
	// stats->reset();
}


// DEFINE VERBOSE MODE

unsigned int VERBOSE_LEVEL = 0;		// defaults off

// Sets the verbose level
void setVerbose(const unsigned int vLevel) {
	if (vLevel > VERBOSE_MAX) 
		printWarning("Verbose level invalid!");
	else 
		VERBOSE_LEVEL = vLevel;
}

// Returns true if the verbose level is set at 'vLevel' or higher
bool verboseLevel(const unsigned int vLevel) {
	if (vLevel > VERBOSE_MAX) {
		printWarning("Verbose level invalid!");
		return false;
	}
	
	if (VERBOSE_LEVEL == VERBOSE_OFF) 
		return false;
	
	// if we're on low verbose mode
	if (VERBOSE_LEVEL == VERBOSE_LOW) {
		if (vLevel == VERBOSE_LOW) 
			return true;
		else
			return false;
	}
	
	// if we're on medium verbose mode
	if (VERBOSE_LEVEL == VERBOSE_MED) {
		if (vLevel == VERBOSE_LOW || vLevel == VERBOSE_MED)
			return true;
		else
			return false;
	}
	
	// if we're on high verbose mode
	if (VERBOSE_LEVEL == VERBOSE_HIGH) {
		if (vLevel == VERBOSE_LOW || vLevel == VERBOSE_MED || vLevel == VERBOSE_HIGH)
			return true;
		else 
			return false;
	}
	
	
	// otherwise
	return false;
}





// CHAINING INFORMATION (FOR BINARY HEAP WITH CHAINING)
string bhWithChainingFilename = "bhWithChainingInfo";

void getBHWithChainingFilename(string &filename) {
	filename = bhWithChainingFilename;
}



void pausePrompt() {
	cout << "Press Enter to Continue...";
	cin.get();
	
	cout << endl;
}



size_t DEBsize;
size_t get_deb_size() {
	return DEBsize;
}

void init_deb_size() {
	DEBsize = 0;
}

void inc_deb_size(size_t inc) {
	DEBsize += inc;
	cout << "DEBsize = " << DEBsize << endl;
}

void dec_deb_size(size_t dec) {
	cout << "DEBsize = " << DEBsize << ", dec = " << dec << endl;
	DEBsize -= dec;
}




bool usingOSI = false;

void setOptimizedSequenceOfInserts(bool value) {
	usingOSI = value;
}

bool usingOptimizedSequenceOfInserts() {
	return usingOSI;
}



bool global_boolean = false;

bool getGlobalBoolean() {
	return global_boolean;
}

void setGlobalBoolean(bool value) {
	global_boolean = value;
}