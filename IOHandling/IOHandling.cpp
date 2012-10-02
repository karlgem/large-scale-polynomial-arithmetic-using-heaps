#include "IOHandling.h"

#include <fstream>

/***** Function Prototypes *****/

void parsePolynomialString (string p, poly_t &result);


/***** Function Implementations *****/

/*
 *	Reads the two sets of polynomials from the input file. The file must be formatted
 *	according to the following guidelines:
 *		- The first line is a number (call it N) indicating the number of polynomial 
 *			pairs in the file.
 *		- The next N lines are polynomials constituting the first set.
 *		- The final N lines are polynomials constituting the second set.
 *
 *	NOTE: The format of the lines containing polynomials is:
 *				p = a1*x^n1 + a2*x^n2 + ... 
 *			where:
 *				- p can be anything WITHOUT the "=" sign (eg f_1, g_2, ...)
 *				- a1, a2, ... are positive integers
 *				- n1, n2, ... are positive integers (in decreasing order)
 *
 *	@param inputName the name of the input file
 *	@param set1 vector that will be filled with the first set of polynomials
 *	@param set2 vector that will be filled with the second set of polynomials
 *	@param numberOfPairs will be set to be the number of pairs
 */
void readInputFile (string inputName, vector<poly_t> &set1, vector<poly_t> &set2, int &numberOfPairs) {
	// open the input file in a filestream
	ifstream in(inputName.c_str());
    string line;

	// get the number of pairs
	getline(in, line);
    numberOfPairs = atoi(line.c_str());

    poly_t p;

	// Read the first set of polynomials
    for (int i = 0; i < numberOfPairs; i++) {

        // read the f polynomial of the i'th pair
        getline(in, line);
        line = line.substr(line.find("=") + 1);
        parsePolynomialString(line, p);


        // add f to its corresponding vectors
        //f_polynomials[i] = p;
		set1.push_back(p);

		string p_string;
		polyToString(p, p_string);

		// reset the polynomials to read new ones
		p.clear();
    }

	// Read the second set of polynomials
    for (int i = 0; i < numberOfPairs; i++) {

        // read the g polynomial of the i'th pair
        getline(in, line);
        line = line.substr(line.find("=") + 1);
        parsePolynomialString(line, p);


        // add g to its corresponding vectors
        //g_polynomials[i] = g;
		set2.push_back(p);

		string p_string;
		polyToString(p, p_string);

		// reset the polynomial to read a new one
		p.clear();
    }
}




/*
 *	Skips the spaces in a string until it reaches the next non-space character.
 *
 *	@param current the current character from which to start skipping spaces
 *	@param index the index of the current character in the string
 *	@param line the actual string
 */
inline void skipSpaces (char &current, size_t &index, string line) {
    while (index <= line.length() && current == ' ') {
        current = line[index++];
    }
}


/*
 *	Parses a string representing a polynomial into a variable of type poly_t.
 *
 *	NOTE: The format of the polynomial string is:
 *				a1*x^n1 + a2*x^n2 + ... 
 *			where:
 *				- a1, a2, ... are positive integers
 *				- n1, n2, ... are positive integers (in decreasing order)
 *
 *	@param p the string holding representing a polynomial
 *	@param result a variable that will hold the parsed polynomial
 */
void parsePolynomialString (string p, poly_t &result) {
    
    p.append("?");      // ? is a sentinel value
    	
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
        
		result.push_back(createMonomial(degree, coef, 0));
        
        
        if (current == '?') {       // done reading poly
            break;
        }
        else {      // current is '+'
            current = p[i++];
            skipSpaces(current, i, p);
        }
	}
	
}

