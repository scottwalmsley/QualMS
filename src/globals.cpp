/*
 * globals.cpp
 * Copyright 2012 Scott Walmsley, Hyungwon Choi and Damian Fermin
 *
 */
 
#include <getopt.h>
#include <iostream>
#include <string>
#include <map>
#include <deque>
#include <cstdlib>
#include <fstream>

#include "boost/multi_array.hpp"
#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "boost/algorithm/string/regex.hpp" // for split_regex() function

using namespace std;
using namespace boost;

namespace fs = filesystem;

// Global variables
int g_STATES_K = -1;
string g_inputDirPath;
int g_BIC_low = 0;
int g_BIC_high = 0;
bool g_runBIC = false;


void print_usage() {
	cerr << "\nUSAGE: qualMS.exe -k -b -d\n\n"
	     << "Command flags:\n"
	     << "   -k           Assume 'k' states.\n"
	     << "   -d <path>    Path to folder containing input matrix files. ** Required input **\n\n"
	     << "   -b <#:#>     Provide an estimate for best number of 'k' states using BIC\n"
	     << "                Example usage: -b 3:10 means consider between 3 and 10 states\n"
	     << "                Using the -b flag, qualMS will only report the best value for K and exit"
	     << "\n\n"; 
}


void parse_command_line_args(int argc, char *argv[]) {
	
	if(argc < 3) {
		print_usage();
		exit(0);
	}
	
	g_inputDirPath = "empty";
	string bic_str;
	int c;
	while( (c = getopt(argc, argv, "k:d:b:")) != -1 ) {
		switch(c) {
		
			case 'k':
				g_STATES_K = atoi(optarg);
				break;
			case 'd':
				g_inputDirPath = optarg;
				break;
			case 'b':
				bic_str = optarg;
				g_runBIC = true;
				break;
			default:
				print_usage();
				exit(0);
		}		
	}
	
	
	string L = bic_str.substr(0, bic_str.find(":") );
	string H = bic_str.substr(bic_str.find(":")+1);
	
	g_BIC_low = atoi( L.c_str() );
	g_BIC_high = atoi( H.c_str() );
	
	if( g_runBIC && (g_BIC_low < 4) ) {
		cerr << endl << "WARNING: The lower bound for option -b is too small.\n"
		     << "Adjusting the lowest estimator for k to be 4\n\n";
		g_BIC_low = 4;
	}
	
	if( !g_runBIC && g_STATES_K < 4) {
		cerr << "\nERROR: k = " << g_STATES_K << " is an invalid value.\n"
		     << "You must provided a number >= 4 for -k\n\n";
		exit(0);
	}
	
	if( g_runBIC && g_STATES_K > 1 ) {
		cerr << "\nERROR: You can't use -k and -b options together.\n\n";
		exit(0);
	}
	
	if(g_inputDirPath == "empty") {
		cerr << "\nERROR: -d option must be provided.\n"
		     << "Please indicate the path to the matrix files\n\n";
		exit(0);
	}
}


// Function reads in each input file and stores it into a multi_array object
void parseInputFiles(int &Ncol, int &Nrow, deque<string> &fileNameDeq) {

    Ncol = 0;
    Nrow = 0;
    
    fs::path p = g_inputDirPath.c_str();
    fs::directory_iterator end_iter;

    if( !fs::exists(p) ) {
    	cerr << "\nERROR: Unable to find source directory " << g_inputDirPath << endl;
    	exit(EXIT_FAILURE);
    }

    if( fs::is_regular_file(p) ) {
    	cerr << "\nERROR: " << g_inputDirPath << " is a regular file. I need a folder\n";
    	exit(EXIT_FAILURE);
    }

    // Record the names of the files we are going to read in and store them
    // into the 'fileNameDeq'
    for(fs::directory_iterator dir_iter(p); dir_iter != end_iter; ++dir_iter) {
    	if(fs::is_regular_file(dir_iter->status())) {
    	
    		// skip files that don't end in .txt.
    		if(dir_iter->path().extension() != ".txt") continue;
    		
    		string s = dir_iter->path().generic_string();
    		
    		fileNameDeq.push_back(s);
    	}
    }

    // Determine the dimensions of the files we are going to read.
    // We are assuming that all of the files have the exact same dimensions.
    ifstream in;
    in.open(fileNameDeq.at(0).c_str());
    if(!in.is_open() ) {
    	cerr << "\nERROR: Unable to open '" << fileNameDeq.at(0) << "'\nExiting now.\n";
    	exit(EXIT_FAILURE);
    }


    vector<string> fields;
    string line;
    while( !in.eof() ) {
    	getline(in, line);
    	if(line.length() < 3) continue;

    	Nrow++;
    	if(Ncol == 0) { // get number of columns
    		split_regex(fields, line, regex( "[\\s|\\t]+" ) );
    		Ncol = (signed) fields.size();
    	}
    }
    in.close();
}



// Function returns the maximum value in a double vector along
// with the index of that value in the vector
void getMaxIndexVec(vector<double> &v, double &maxVal, int &index) {
	index = 0;
	maxVal = v.at(0);

	for(int i = 0; i < (signed) v.size(); i++) {
		if(v.at(i) > maxVal) {
			maxVal = v.at(i);
			index = i;
		}
	}
}


// Function returns the maximum value obtained from the difference of two
// equal length vectors
double getMaxAbsVectorDiff(vector<double> &A, vector<double> &B) {
	double ret = 0;
	double delta;
	double abs_delta;
	int N = (signed) A.size(); // A.size() must equal B.size()
	
	for(int i = 0; i < N; i++) {
		delta = A.at(i) - B.at(i);
		abs_delta = fabs(delta);
		
		if(abs_delta > ret) ret = abs_delta;
	}
	return ret;
}



// function to round doubles to N number of decimal places
double round_dbl(double r, int places) {
        double off = pow((double)10, places);
        return ( round(r * off) / off );
}





