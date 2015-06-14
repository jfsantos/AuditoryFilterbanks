#include "../CochlearFilterbank/ModulationFilterBank.h"
#include "mex.h"
#include <iostream>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]){

	// Parsing input parameters (WARNING: no error checking!)
	double *matx, *mcf;
	matx = mxGetPr(prhs[0]);
	size_t x_mrows = mxGetM(prhs[0]);
	size_t x_mcols = mxGetN(prhs[0]);
	mcf = mxGetPr(prhs[1]);
	size_t cf_mrows = mxGetM(prhs[1]);
	size_t cf_mcols = mxGetN(prhs[1]);
	double fs = mxGetScalar(prhs[2]);
	double q = mxGetScalar(prhs[3]);
	
	// Mapping MATLAB array to Eigen::VectorXd
	Eigen::Map<Eigen::VectorXd> x(matx,x_mrows);
	Eigen::Map<Eigen::VectorXd> cf(mcf,cf_mrows);
	int num_channels = cf.rows();
	//mexPrintf("%d", num_channels);

	ModulationFilterBank mf = ModulationFilterBank(fs, num_channels, cf, q);
	
	// Mapping MATLAB output matrix to Eigen::MatrixXd
	plhs[0] = mxCreateDoubleMatrix(x_mrows, num_channels, mxREAL);
	Eigen::Map<Eigen::MatrixXd> y(mxGetPr(plhs[0]),x_mrows,num_channels);
	
	// Running filters
	y = mf.process(x, int(x_mrows));
	
}

