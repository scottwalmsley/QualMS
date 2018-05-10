/*
 * qualMRFClass.hpp  Copyright 2012 Scott Walmsley, Hyungwon Choi and Damian Fermin
 *
 * Definition of the Qual MRF object class
 */

#include<iostream>
#include<string>
#include<vector>

#include "structs.hpp"

#include "boost/multi_array.hpp"


using namespace std;
using namespace boost;

class QualMRFClass;


class QualMRFClass {
	private:
		int Nrow;
		int Ncol;
		string srcFile;
		double maxLogVal;  // maximum value observed in 'logMAT' multi_array
		vector<double> theta;
		multi_array<double, 2> rawMAT; // raw matrix
		multi_array<double, 2> logMAT; // log-transformed matrix
		multi_array<double, 2> stateMAT; // state assignment matrix
		multi_array<double, 3> zMAT;   // Z matrix of states
		
		vector<mixModelStruct> modelParamVec; // statistical parameters 

		vector<double> getStateCount(int r, int c); // private function

	public:
		QualMRFClass() {};
		QualMRFClass(string &d, int R, int C); // default contructor
		
		void logTransformMatrix();
		void initializeModelParams();
		void fitEM();
		void chooseK();
		void rerankParams();
		void initializeMRF();
		void ICM();
		void optimizeTheta();
		void updateStates();
		void writeFile();

		double calcBIC();
		double likelihood(const vector<double> &x, vector<double> &grad/*, vector<double> &stateVec*/);


};

