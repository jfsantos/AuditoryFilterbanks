#include "CochlearFilterbank.h"
#include "mex.h"
#include <iostream>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]){

	// Parsing input parameters (WARNING: no error checking!)
	double *matx;
	matx = mxGetPr(prhs[0]);
	size_t mrows = mxGetM(prhs[0]);
	size_t mcols = mxGetN(prhs[0]);
	double fs = mxGetScalar(prhs[1]);
	int num_channels = static_cast<int>(mxGetScalar(prhs[2]));
	double low_freq = mxGetScalar(prhs[3]);

	CochlearFilterbank cf = CochlearFilterbank(fs, num_channels, low_freq);
	
	// Mapping MATLAB array to Eigen::VectorXd
	Eigen::Map<Eigen::VectorXd> x(matx,mrows);

	// Mapping MATLAB output matrix to Eigen::MatrixXd
	plhs[0] = mxCreateDoubleMatrix(mrows, num_channels, mxREAL);
	Eigen::Map<Eigen::MatrixXd> y(mxGetPr(plhs[0]),mrows,num_channels);
	
	// Running filters
	y = cf.process(x, int(mrows));
}

