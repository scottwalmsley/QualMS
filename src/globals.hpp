/*
 * globals.hpp
 *
 * Global functions and variables are defined here.
 *
 */

#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_

#include <string>
#include <map>
#include <deque>

#include "boost/multi_array.hpp"

using namespace std;
using namespace boost;

// Global variables
extern int g_STATES_K; // number of states 'K' to consider
extern string g_inputDirPath;
extern int g_BIC_low;
extern int g_BIC_high;
extern bool g_runBIC;

void print_usage();

void parse_command_line_args(int argc, char *argv[]);

void parseInputFiles(int &Ncol, int &Nrow, deque<string> &fileNameDeq);

void getMaxIndexVec(vector<double> &v, double &maxVal, int &index);

double getMaxAbsVectorDiff(vector<double> &A, vector<double> &B);

double round_dbl(double r, int places);

#endif /* GLOBALS_HPP_ */
