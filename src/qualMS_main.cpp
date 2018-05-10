/*
 * qualMS_main.cpp
 *
 *  
 *     Copyright 2012 Scott Walmsley, Hyungwon Choi and Damian Fermin    
 */

#include <iostream>
#include <string>
#include <map>
#include <deque>

#include "globals.hpp"
#include "QualMRFClass.hpp"

#include "boost/multi_array.hpp"

using namespace std;
using namespace boost;



int main(int argc, char *argv[]) {

	parse_command_line_args(argc, argv);
	
	int Ncol = 0;
	int Nrow = 0;

	deque<string> F;

	QualMRFClass *Q = NULL;

	// get the names of the matrix files and the matrix dimensions
	parseInputFiles(Ncol, Nrow, F);
	
	if( !g_runBIC ) {
		cerr << endl << "Number of states to use: " << g_STATES_K << endl 
			 << "Assumed matrix dimensions\n# columns: " << Ncol
			 << "\n# rows: " << Nrow << endl << endl;
	}
	
	// Iterate through the files and load them into a boost multi_array
	for(deque<string>::iterator d = F.begin(); d != F.end(); d++) {

		Q = new QualMRFClass(*d, Nrow, Ncol);
		
		Q->logTransformMatrix();
		if(g_runBIC) {
			Q->chooseK();
		}
		else {
			Q->initializeModelParams();
			Q->fitEM();
			Q->rerankParams();
			Q->initializeMRF();
			Q->ICM();
			Q->writeFile();
		}
		
		delete(Q); Q = NULL;
	}




	return 0;
}

