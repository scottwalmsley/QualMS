/*
 * QualMRFClass.cpp  Copyright 2012 Scott Walmsley, Hyungwon Choi and Damian Fermin
 *
 * Functions of the QualMRFClass object
 */

#include<iostream>
#include<string>
#include<cstdlib>
#include<fstream>
#include<vector>
#include<functional>

#include "globals.hpp"
#include "QualMRFClass.hpp"
#include "statFunctions.hpp"

#include "nlopt_wrapper.hpp" // NLOPT optimization library

#include "boost/multi_array.hpp"
#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "boost/algorithm/string/regex.hpp" // for split_regex() function

using namespace std;
using namespace boost;

namespace fs = filesystem;



// Constructor
QualMRFClass::QualMRFClass(string &d, int R, int C) {
	ifstream in;
	vector<string> fields;
	string line, x;
	double val;
	int c = 0;
	int r = 0;
	
	// set these object variables
	Nrow = R;
	Ncol = C;

	// define the dimensions of the matrices
	rawMAT.resize(boost::extents[Nrow][Ncol]);


	// record the source file
	srcFile = d;

	in.open(d.c_str());
	if(!in.is_open() ) {
		cerr << "\nERROR: Unable to open '" << d << "'\nExiting now.\n";
		exit(EXIT_FAILURE);
	}

	if(g_runBIC) cout << "Estimating optimal states for ";
	else cout << endl << "Parsing " ;
		
	cout << srcFile.substr( (srcFile.find_last_of("/")+1) ) << endl;



	while( !in.eof() ) {
		getline(in, line);
		if(line.length() < 3) continue; // likely a useless line

		// split current line into a string vector
		fields.clear();
		split_regex(fields, line, regex( "[\\s|\\t]+" ) );

		for(c = 0; c < Ncol; c++) {
			x = fields.at(c);

			if(x == "NA") val = 0;
			else val = atof( x.c_str() );

			rawMAT[ r ][ c ] = val;

		}
		r++; // increment row index
	}
	in.close();
}



// Function to log transform the data in rawMAT and store it in a normalized matrix
void QualMRFClass::logTransformMatrix() {

	maxLogVal = 0; // class variable initialized

	logMAT.resize(boost::extents[Nrow][Ncol]);

	for(int r = 0; r < Nrow; r++) {
		for(int c = 0; c < Ncol; c++) {
			double x = rawMAT[ r ][ c ] + 1.0;
			logMAT[ r ][ c ] = log(x);

			if(log(x) > maxLogVal) maxLogVal = log(x);
		}
	}
}



// Function intitializes the model parameter vector with the most likely starting values
void QualMRFClass::initializeModelParams() {
	
	double k = (double) g_STATES_K;
	
	modelParamVec.resize(g_STATES_K);

	for(double i = 0; i < g_STATES_K; i++) {
		modelParamVec[i].pi = 1.0 / k;
		modelParamVec[i].sigma2 = 1.0;
		modelParamVec[i].mu = ( (i+1)/(g_STATES_K + 1) ) * maxLogVal;
	}

	//initalize Z-matrix object
	zMAT.resize(boost::extents[Nrow][Ncol][g_STATES_K]);

	for(int r = 0; r < Nrow; r++) {
		for(int c = 0; c < Ncol; c++) {
			for(int i = 0; i < g_STATES_K; i++) {
				zMAT[ r ][ c ][ i ] = 1.0 / k;
			}
		}
	}
}



// Function fits the mixture model using the data in the QualMRFClass
void QualMRFClass::fitEM() {

	double prob, X, state_sum;
	double tol = 100;
	int cur_iter = 0;

	vector<mixModelStruct> oldParamVec;

	while( (tol > 1e-5) && (cur_iter < 50) ) {
		
		cur_iter++;

		// record the paramters from the last iteration.
		// if this is the first iteration, just copy the initial values
		oldParamVec.clear();
		oldParamVec = modelParamVec;

		// E-STEP
		for(int r = 0; r < Nrow; r++) {
			for(int c = 0; c < Ncol; c++) {

				state_sum = 0;
				for(int k = 0; k < g_STATES_K; k++) {
					X = logMAT[ r ][  c ];
					prob = dnorm(X, modelParamVec.at(k).mu, modelParamVec.at(k).sigma2);
					zMAT[ r ][ c ] [ k ] = modelParamVec.at(k).pi * prob;
					
					// if the cell is zero, assume state vector of 1,0,0,0
					if(X == 0) zMAT[ r ][ c ][ k ] = (k == 0 ? 1 : 0);
					
					state_sum += zMAT[ r ][ c ][ k ];
				}

				// normalize values in zMAT so that the state values 'k' sum to one
				for(int k = 0; k < g_STATES_K; k++) {
					X = zMAT[ r ][ c ][ k ];
					zMAT[ r ][ c ][ k ] =  X / state_sum; 
				}
			}
		}

		// M-STEP
		double numer, z_sum, diff, sq_diff;
		for(int k = 0; k < g_STATES_K; k++) {

			// Get the z_sum for the current state. This value
			// is constant for each state iteration
			z_sum = 0;
			for(int r = 0; r < Nrow; r++) {
				for(int c = 0; c < Ncol; c++) {
					if(logMAT[ r ][ c ] > 0) z_sum += zMAT[ r ][ c ][ k ];
				}
			}


			// compute the new mean (mu) values for this state
			numer = 0;
			for(int r = 0; r < Nrow; r++) {
				for(int c = 0; c < Ncol; c++) {
					numer += logMAT[ r ][ c ] * zMAT[ r ][ c ][ k ];
				}
			}
			modelParamVec.at(k).mu = numer / z_sum;
			
			// compute the new variance (sigma2) values for this state
			numer   = 0;
			diff    = 0;
			sq_diff = 0;
			for(int r = 0; r < Nrow; r++) {
				for(int c = 0; c < Ncol; c++) {
					X = logMAT[ r ][ c ];
					if(X > 0) {
						diff = X - modelParamVec.at(k).mu;
						sq_diff = pow(diff, 2.0);
						numer += sq_diff * zMAT[ r ][ c ][ k ];
					}
				}
			}
			modelParamVec.at(k).sigma2 = numer / z_sum;

			// compute the new mixture proportions (pi) values for this state
			numer = 0;
			X = 0;
			for(int r = 0; r < Nrow; r++) {
				for(int c = 0; c < Ncol; c++) {
					if(logMAT[ r ][ c ] > 0) {
						numer += zMAT[ r ][ c ][ k ];
						X++; // use X here to normalize 'pi' at end of this double loop
					}
				}
			}
			modelParamVec.at(k).pi = numer / X;

		} // end loop over g_STATES_K


		// compute the difference betweent the values in oldParamVec and modelParamVec
		double delta_pi = 0;
		double delta_mu = 0;
		double delta_sigma2 = 0;
		double maxVal = -1;
		for(int k = 0; k < g_STATES_K; k++) {

			delta_pi = fabs( oldParamVec.at(k).pi - modelParamVec.at(k).pi );
			delta_mu = fabs( oldParamVec.at(k).mu - modelParamVec.at(k).mu );
			delta_sigma2 = fabs( oldParamVec.at(k).sigma2 - modelParamVec.at(k).sigma2 );
			
			X = get_max_3(delta_pi, delta_mu, delta_sigma2);
			
			/* //debugging code
			cout << k << "\t(pi, mu, sigma2) = ("
				 << modelParamVec.at(k).pi << ", "
				 << modelParamVec.at(k).mu << ", "
				 << modelParamVec.at(k).sigma2 << ") X = " << X << endl;
			*/

			if(X > maxVal) maxVal = X;

		}
		tol = maxVal;
	} // end while loop
}




// Function to select the optimal number of states 'k' for this data type
void QualMRFClass::chooseK() {
	
	int cur_K = g_BIC_low;
	map<int, double> BICmap;
	double lowestBIC;
	int bestK;
	
	while(cur_K <= g_BIC_high) {
		g_STATES_K = cur_K;
		initializeModelParams();
		fitEM();
		BICmap[ cur_K ] = calcBIC();
		cur_K++;
	}
	
	map<int, double>::iterator m;
	lowestBIC = 1e10;
	bestK = 0;
	for(m = BICmap.begin(); m != BICmap.end(); m++) {
		if(m->second < lowestBIC) {
			lowestBIC = m->second;
			bestK = m->first;
		}
		cout << "k = " << m->first << "\tBIC = " << m->second << endl;
	}
	cout << endl << "** Lowest K = " << bestK << "\tBIC = " << lowestBIC << endl << endl;
}


// Function to calculate Bayesian Information Criterion (BIC)
double QualMRFClass::calcBIC() {
	
	double tmp;
	double logLik = 0;
	double pi, mu, sigma2;
	double ret = 0;
	
	for(int r = 0; r < Nrow; r++) {
		for(int c = 0; c < Ncol; c++) {
			tmp = 0;
			for(int k = 0; k < g_STATES_K; k++) {
				pi = modelParamVec.at(k).pi;
				mu = modelParamVec.at(k).mu;
				sigma2 = modelParamVec.at(k).sigma2;
				tmp += pi * dnorm( logMAT[ r ][ c ], mu, sigma2 );
			}
			logLik += log(tmp);
			
		}
	}
	
	ret = ( (3*g_STATES_K) * log(Nrow * Ncol) ) -2*logLik;
	return ret;
}




// Function reorders the data in modelParamVec based upon the mu
void QualMRFClass::rerankParams() {
	
	list<double> rawMU_list;
	vector<mixModelStruct> newModelVec;
	
	vector<mixModelStruct>::iterator v;
	list<double>::iterator cur_mu;
	
	
	// record the mu's
	for(v = modelParamVec.begin(); v != modelParamVec.end(); v++) rawMU_list.push_back( v->mu );
	
	rawMU_list.unique(); // keep only 1 instance of each mu
	rawMU_list.sort(); // sorted low to high
	
	for(cur_mu = rawMU_list.begin(); cur_mu != rawMU_list.end(); cur_mu++) { // iterate over ordered mu's
		for(v = modelParamVec.begin(); v != modelParamVec.end(); v++) { // iterate over parameter structs
			if( v->mu == *cur_mu ) newModelVec.push_back( *v );
		}
	}
	
	/******************* // for debugging
	cout << "OLD(mu, pi, sigma2)\tNEW(mu, pi,sigma2)\n";
	for(int i = 0; i < g_STATES_K; i++) {
		cout << "(" << modelParamVec.at(i).mu << ", "
		     << modelParamVec.at(i).pi << ", "
		     << modelParamVec.at(i).sigma2 << ")\t("
		     << newModelVec.at(i).mu << ", "
		     << newModelVec.at(i).pi << ", "
		     << newModelVec.at(i).sigma2 << ")\n";
	}
	*****************************************************/
	
	modelParamVec = newModelVec;
}



// Function initializes MRF States model
void QualMRFClass::initializeMRF() {

	double X, maxState, maxVal;

	stateMAT.resize(boost::extents[Nrow][Ncol]);

	// initialize stateMAT assigning it the max. state
	for(int r = 0; r < Nrow; r++) {
		for(int c = 0; c < Ncol; c++) {
			X = logMAT[ r ][ c ];
			if(X == 0) stateMAT[ r ][ c ] = 0;
			else {

				maxState = 0;
				maxVal = 0;
				for(int k = 0; k < g_STATES_K; k++) {
					if(zMAT[ r ][ c ][ k ] > maxVal) {
						maxState = k;
						maxVal = zMAT[ r ][ c ][ k ];
					}
				}
				stateMAT[ r ][ c ] = maxState;
			}
		}
	}
}



// Function computes the ICM for this QualMRFClass object
// This function is the big boy
void QualMRFClass::ICM() {
	
	double threshold = 0.01;
	double diff = 1000;
	int iter = 0;
	
	// Use this vector to hold the values of theta from the previous iteration
	// this is initially set to 1 just like theta.
	vector<double> oldTheta( (g_STATES_K + 1), 1 );

	theta.resize( (g_STATES_K + 1) );
	for(int i = 0; i < (g_STATES_K + 1); i++) theta[i] = 1;


	while( (diff > threshold) && (iter < 100) ) {
		oldTheta = theta;


		// update theta
		optimizeTheta();

		// update Z
		updateStates();

		diff = getMaxAbsVectorDiff(theta, oldTheta);

	} // end while loop
	
}




// Function returns the counts of the states for the neighboring data points
// This is a private function and cannot be called outside of this class
vector<double> QualMRFClass::getStateCount(int r, int c) {

	vector<double> ret(g_STATES_K, 0); 
	vector<double> d(g_STATES_K, 0); 

	// matrix boundaries in 0-based coordinates
	int rB = Nrow - 1;
	int cB = Ncol - 1;

	// compute the coordinates for all of the neighbors for 'r' and 'c'
	int rr, cc;
	double ntotal = 0;
	for(rr = r-1; rr <= r+1; rr++) {
		for(cc = c-1; cc <= c+1; cc++) {
			if((rr >= 0) && (rr <= rB) && (cc >= 0) && (cc <= cB) && (rr != r && cc != c) ) {
				d.at( stateMAT[ rr ][ cc ] )++;
				ntotal++;
			}
		}
	}

	for(int k = 0; k < g_STATES_K; k++) {
		double x = (ntotal - d.at(k)) / ntotal;
		ret.at(k) = x;
	}
	
	return ret;
}



// Function computes the likelihood function
double QualMRFClass::likelihood(const vector<double> &x, vector<double> &grad) {

	double res = 0;
	double tmp = 0;
	double S = 0;
	vector<double> counts(g_STATES_K, 0); // holds counts for each state
	vector<double> d(g_STATES_K, 0);	

	for(int r = 0; r < Nrow; r++) {
		for(int c = 0; c < Ncol; c++) {
			
			counts = getStateCount(r,c);
			tmp = 0;
			for(int k = 0; k < g_STATES_K; k++) {
				tmp += exp(x[ k ] - (x[ g_STATES_K ]  * counts[ k ]) );
			}
			S = stateMAT[ r ][ c ];
			res += (x[ S ] - x[ g_STATES_K ] * counts[ S ]);
			res -= log(tmp);
		}
	}
	return res;
}




// Function to optimize the values stored in theta
void QualMRFClass::optimizeTheta() {

	vector<double> upperB(g_STATES_K+1, 10);
	vector<double> lowerB(g_STATES_K+1, -10);
	lowerB[ g_STATES_K ] = 0;
	
	//initialize nlopt::opt object
	nlopt::opt opt(nlopt::LN_COBYLA, (g_STATES_K + 1));
	
	opt.set_lower_bounds(lowerB);
	opt.set_upper_bounds(upperB);

	opt.set_ftol_abs(1e-4);
	opt.set_maxeval(1e3);

	vector<double> x = theta;
	
	vector<double> stateVec( (Nrow*Ncol), 0 );

	int i = 0;
	for(int r = 0; r < Nrow; r++) {
		for(int c = 0; c < Ncol; c++) {
			i = (r * Ncol) + c;
			stateVec[ i ] = stateMAT[ r ][ c ];
		}
	}

	nlopt_wrap(opt, x, &QualMRFClass::likelihood, this);
	
	theta = x;
	
	// Let's the user know something is happening 
	for(int k = 0; k < (g_STATES_K + 1); k++) cout << theta.at(k) << "\t";
	cout << endl;
}




// Function updates the states and probabilities after you optimize Theta
void QualMRFClass::updateStates() {

	double x, j, maxProb, sumProb;
	int idx;
	vector<double> prob(g_STATES_K, 0);
	vector<double> counts(g_STATES_K, 0); // holds counts for each state
	

	for(int r = 0; r < Nrow; r++) {
		for(int c = 0; c < Ncol; c++) {
			counts = getStateCount(r,c);
			
			for(int k = 0; k < g_STATES_K; k++) {
				x = theta[ k ] - (theta[ g_STATES_K ]  * counts[ k ]);
				j = log_dnorm( logMAT[ r ][ c ], modelParamVec.at(k).mu, modelParamVec.at(k).sigma2 );
				prob[ k ] = x + j;
			}

			getMaxIndexVec(prob, maxProb, idx); // get max. value in prob
			
			for(int k = 0; k < g_STATES_K; k++) {
				x = prob.at(k);
				prob.at(k) = exp( x - maxProb );
			}

			sumProb = 0;
			for(int k = 0; k < g_STATES_K; k++) sumProb += prob.at(k);
			
			for(int k = 0; k < g_STATES_K; k++) {
				x = prob.at(k);
				prob.at(k) = x / sumProb;
			}

			if(logMAT[ r ][ c ] == 0) {
				stateMAT[ r ][ c ] = 0;
				zMAT[ r ][ c ][ 0 ] = 1;
				for(int k = 1; k < g_STATES_K; k++) zMAT[ r ][ c ][ k ] = 0;
			}
			else {
				stateMAT[ r ][ c ] = idx;
				for(int k = 0; k < g_STATES_K; k++) {
					zMAT[ r ][ c ][ k ] = prob.at(k);
				}
			}
		}
	}

}


// Write output file
void QualMRFClass::writeFile() {
	
	// Create output file's name
	int idx = srcFile.find_last_of(".");
	string outFile = srcFile.substr(0,idx) + ".mrf";
	
	cerr << "Writing results to '" << outFile << "'\n";

	// open output file
	fstream outF;
	outF.open(outFile.c_str(), ios::out);
	if( !outF.is_open() ) {
		cerr << "\nERROR: Unable to create output file " << outFile << endl
			 << "Exiting now...\n\n";
		exit(EXIT_FAILURE);
	}

	// construct header line
	outF << "index\trow\tcol\trawCount\tlogCount\topt_state\tstate_probs\n";

	string line;
	double d;
	idx = 1; 
	for(int r = 0; r < Nrow; r++) {
		for(int c = 0; c < Ncol; c++) {
			outF << idx << "\t"
				 << r << "\t"
				 << c << "\t"
				 << rawMAT[ r ][ c ] << "\t"
				 << logMAT[ r ][ c ] << "\t"
				 << stateMAT[ r ][ c ] << "\t";
		
			for(int k = 0; k < g_STATES_K; k++) {
				d = round_dbl( zMAT[ r ][ c ][ k ], 4 );
				
				outF << d;
				//line += boost::lexical_cast<string>( d );

				if(k < (g_STATES_K - 1)) outF << ",";
			}
			
			outF << endl;
			idx++; // increment index
		}
	}

	outF.close();
}





