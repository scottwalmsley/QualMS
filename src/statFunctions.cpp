/*
 * statFunctions.cpp
 *
 * All generic statistical functions used by qualMS
 */

#include<cmath>
#include<cstdlib>
#include<iostream>

#include "statFunctions.hpp"

using namespace std;

const double PI = 4.0 * atan(1);


// Guassian density function
double dnorm(double x, double mean, double sigma2) {
	double ret = 0;

	double fx = -0.5 * pow( (x-mean), 2.0 ) / sigma2 - 0.5 * log(2.00 * PI * sigma2);

	ret = exp(fx);
	return ret;
}


// LOG Guassian density function
double log_dnorm(double x, double mean, double sigma2) {
	double ret = 0;

	ret = -0.5 * pow( (x-mean), 2.0 ) / sigma2 - 0.5 * log(2.00 * PI * sigma2);
	
	return ret;
}

// Function returns the maximum value between the 3 passed values
double get_max_3(double x, double y, double z) {
	double ret = x; // assume x is the maximum

	if(y > ret) ret = y;

	if(z > ret) ret = z;

	return ret;
}


