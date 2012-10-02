#include "Polynomial.h"

#include <iostream>
#include <string>
#include <sstream>

using namespace std;


void printPolynomial(const poly_t &p) {
    for (size_t i = 0; i < p.size(); i++) {
        
        if (i != 0) {
            cout << " + ";
        }
        

        if (GET_COEF(p[i]) != 0) {
            cout << GET_COEF(p[i]);
        }
        
        if (GET_DEGREE(p[i]) > 0) {
            cout << "*x";
        }
        
        if (GET_DEGREE(p[i]) > 1) {
            cout << "^" << GET_DEGREE(p[i]);
        }
    }
    
	cout << endl;
}


void polyToString(const poly_t &p, string &output) {
		
    stringstream ss (stringstream::in | stringstream::out);
        
    output = "";
    
    for (size_t i = 0; i < p.size(); i++) {
        
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



/*
 *	Returns the degree of the polynomial
 */
deg_t deg(const poly_t &p) {
	return GET_DEGREE(p[0]);
}