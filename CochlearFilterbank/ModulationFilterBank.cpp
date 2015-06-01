#include "ModulationFilterBank.h"

//#define DEBUG

#ifdef DEBUG
#include <iostream>
#endif

using namespace std;

ModulationFilterBank::ModulationFilterBank(double fs, int num_channels, Eigen::VectorXd mf, double q)
{
	filters = makeModulationFilters(fs, num_channels, mf, q);
}


ModulationFilterBank::~ModulationFilterBank(void)
{
}

Eigen::VectorXd ModulationFilterBank::compute_modulation_cfs(double min_cf, double max_cf, int n)
{
	double spacing_factor = pow(max_cf/min_cf,1.0/(n-1));
	Eigen::VectorXd cfs(n);
	cfs[0] = min_cf;
	for (int k = 1; k < n; k++)
		cfs[k] = cfs[k - 1]*spacing_factor;
	return cfs;
}

Eigen::MatrixXd ModulationFilterBank::process(const Eigen::VectorXd &input, int n)
{
	Eigen::MatrixXd y(input.rows(), filters.size());
	vector<future<Eigen::VectorXd>> futures;
	for (unsigned int ch=0; ch < filters.size(); ch++)
		futures.push_back(async(&ModulationFilterBank::process_channel, this, input, n, ch));
	for (unsigned int ch=0; ch < filters.size(); ch++)
	{
		y.col(ch) = futures[ch].get();
	}
	return y;

}

Eigen::VectorXd ModulationFilterBank::process_channel(const Eigen::VectorXd &input, int n, int ch)
{
	return filters[ch].process(input, n);
}

std::vector<Biquad> ModulationFilterBank::makeModulationFilters(double fs, int num_channels, Eigen::VectorXd mf, double q)
{
	vector<Biquad> filters = vector<Biquad>(num_channels);
	for (int k=0; k < num_channels; k++)
	{
    double w0 = 2*M_PI*mf(k)/fs;
		double W0 = tan(w0/2);
		double W02 = pow(W0,2);
		double B0 = W0/q;

		double a0 = 1 + B0 + W02;
		double a1 = 2*W02 - 2;
		double a2 = 1 - B0 + W02;

		Eigen::Matrix<double, 5, 1> f;

		f << B0, 0, -B0, a1, a2;
#ifdef DEBUG
		cout << "[" << k << "]: " << f << endl;
#endif
		filters[k] = Biquad(f/a0);
	}
	return filters;
}
