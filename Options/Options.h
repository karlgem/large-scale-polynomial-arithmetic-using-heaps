// Options.h

/* 
 *	This file defines options and modes (configurable by the user) that
 *	can be executed at runtime. These include: verbose, debugging, etc..
 * 	It also includes profiling tools: timer.
 */

#ifndef Options_h
#define Options_h

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <unistd.h>
#include <vector>
#include "../Polynomial/Monomial.h"

// check to see if compilation was done with/without STXXL


// UNCOMMENT FOR STXXL (COMMENT THE SECOND TYPEDEF)
// #include <stxxl.h>
// typedef stxxl::VECTOR_GENERATOR<w_type, 1, 16*256, 4096>::result heap_t;

typedef std::vector<w_type> heap_t;



typedef heap_t::iterator mem_t;

using namespace std;


// Define error and warning //printing methods
void printError(string err);

void printWarning(string warning);


// HANDLING OUTPUT FILES (time and multiplication output)
void getTimeFilename(string &filename);
void setTimeFilename(const string filename);

bool printMultOutput();
void printMultOutput(bool value);
void getOutputFilename(string &filename);
void setOutputFilename(const string filename);




// CHOOSING WHICH HEAP TO USE
#define BINARY_HEAP 0
#define FUNNEL_HEAP 1
#define FUNNEL_HEAP_WITH_MERGING 2
#define BINARY_HEAP_WITH_CHAINING 3

// chooses the heap to be used in computation
void chooseHeap(string heap_name);

// returns the ID of the previously chosen heap
int chosenHeap ();


// PROFILING OPTIONS
void startTimer();
void continueTimer();
void pauseTimer();
void stopTimer();

void resetStats();
double getElapsedTime();
double getIOWaitTime();
double getReads();
double getWrites();




// DEFINE VERBOSE MODE

#define VERBOSE_OFF 0	// verbose turned off
#define VERBOSE_LOW 1		// first verbose level
#define VERBOSE_MED 2		// second verbose level
#define VERBOSE_HIGH 3		// third verbose level

#define VERBOSE_MAX 3	// keep this updated to the max possible verbose level

// Sets the verbose level
void setVerbose(const unsigned int vLevel);

// Returns true if the verbose level is set at 'vLevel' or higher
bool verboseLevel(const unsigned int vLevel);



// CHAINING INFORMATION (FOR BINARY HEAP WITH CHAINING)
void getBHWithChainingFilename(string &filename);
void parseInputFileForBHWithChaining();


void pausePrompt();


size_t get_deb_size();
void init_deb_size();
void inc_deb_size(size_t inc);
void dec_deb_size(size_t dec);


void setOptimizedSequenceOfInserts(bool value);
bool usingOptimizedSequenceOfInserts();


// global boolean
bool getGlobalBoolean();
void setGlobalBoolean(bool value);

#endif