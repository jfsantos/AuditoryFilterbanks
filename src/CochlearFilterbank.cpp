#include "CochlearFilterbank.h"

//#define DEBUG
#define MULTITHREAD

#ifdef DEBUG
#include <iostream>
#endif

using namespace std;

static double const EAR_Q = 9.26449;				//  Glasberg and Moore Parameters
static double const MIN_BW = 24.7;
static int const ORDER = 1;

CochlearFilterbank::CochlearFilterbank(double fs, int num_channels, double low_freq)
{
	filters = makeERBFilters(fs, num_channels, low_freq);
}


CochlearFilterbank::~CochlearFilterbank(void)
{

}

vector< vector<Biquad> >CochlearFilterbank::makeERBFilters(double fs, int num_channels, double low_freq)
{
	double T = 1.0/fs;
	Eigen::VectorXd ERB(num_channels);
	Eigen::VectorXd B(num_channels);
	Eigen::VectorXd cf = ERBspace(low_freq, fs/2, num_channels);
	double A0 = T;
	double A2 = 0;
	double B0 = 1;
	Eigen::VectorXd B1(num_channels);
	Eigen::VectorXd B2(num_channels);
	Eigen::VectorXd A11(num_channels);
	Eigen::VectorXd A12(num_channels);
	Eigen::VectorXd A13(num_channels);
	Eigen::VectorXd A14(num_channels);
	Eigen::VectorXd gain(num_channels);

	vector< vector<Biquad> > filter_bank = vector< vector<Biquad> >(num_channels);

	complex<double> i (0,1);
	complex<double> aux1, aux2, aux3, aux4, aux5, aux6;

	for (int k=0; k < num_channels; k++)
	{
		ERB[k] = pow(pow((cf[k]/EAR_Q),ORDER) + pow(MIN_BW,ORDER),1/ORDER);
		B[k] = 1.019*2*M_PI*ERB[k];
		B1[k] = -2.0*cos(2*cf[k]*M_PI*T)/exp(B[k]*T);
		B2[k] = exp(-2*B[k]*T);
		A11[k] = -(-B1[k]*T + 2.0*sqrt(3.0+pow(2.0,1.5))*T*sin(2.0*cf[k]*M_PI*T)/exp(B[k]*T))/2.0;
		A12[k] = -(-B1[k]*T - 2.0*sqrt(3.0+pow(2.0,1.5))*T*sin(2.0*cf[k]*M_PI*T)/exp(B[k]*T))/2.0;
		A13[k] = -(-B1[k]*T + 2.0*sqrt(3.0-pow(2.0,1.5))*T*sin(2.0*cf[k]*M_PI*T)/exp(B[k]*T))/2.0;
		A14[k] = -(-B1[k]*T - 2.0*sqrt(3.0-pow(2.0,1.5))*T*sin(2.0*cf[k]*M_PI*T)/exp(B[k]*T))/2.0;

		aux1 = complex<double>(-2)*exp(complex<double>(4)*i*cf[k]*M_PI*T)*T;
		aux2 = complex<double>(2)*exp(-(B[k]*T) + complex<double>(2)*i*cf[k]*M_PI*T)*T*(cos(2*cf[k]*M_PI*T) - sqrt(3.0 - pow(2.0,1.5))*sin(2*cf[k]*M_PI*T));
		aux3 = complex<double>(2)*exp(-(B[k]*T) + complex<double>(2)*i*cf[k]*M_PI*T)*T*(cos(2*cf[k]*M_PI*T) + sqrt(3.0 - pow(2.0,1.5))*sin(2*cf[k]*M_PI*T));
		aux4 = complex<double>(2)*exp(-(B[k]*T) + complex<double>(2)*i*cf[k]*M_PI*T)*T*(cos(2*cf[k]*M_PI*T) - sqrt(3.0 + pow(2.0,1.5))*sin(2*cf[k]*M_PI*T));
		aux5 = complex<double>(2)*exp(-(B[k]*T) + complex<double>(2)*i*cf[k]*M_PI*T)*T*(cos(2*cf[k]*M_PI*T) + sqrt(3.0 + pow(2.0,1.5))*sin(2*cf[k]*M_PI*T));
		aux6 = pow(-2 / exp(2*B[k]*T) - complex<double>(2)*exp(complex<double>(4)*i*cf[k]*M_PI*T) + complex<double>(2)*(complex<double>(1) + exp(complex<double>(4)*i*cf[k]*M_PI*T))/exp(B[k]*T),4);
		gain[k] = abs(abs((aux1 + aux2)*(aux1 + aux3)*(aux1 + aux4)*(aux1 + aux5)/aux6));
		vector<Biquad> filters = vector<Biquad>(4);
		Eigen::Matrix<double, 5, 1> coeffs1, coeffs2, coeffs3, coeffs4;
		coeffs1 << A0*(1.0/gain[k]), A11[k]*(1.0/gain[k]), A2*(1.0/gain[k]), B1[k], B2[k];
		coeffs2 << A0, A12[k], A2, B1[k], B2[k];
		coeffs3 << A0, A13[k], A2, B1[k], B2[k];
		coeffs4 << A0, A14[k], A2, B1[k], B2[k];
#ifdef DEBUG
		cout << "Coeffs for channel 1:" << endl;
		cout << "Filter 1" << endl << coeffs1 << endl << endl;
		cout << "Filter 2" << endl << coeffs2 << endl << endl;
		cout << "Filter 3" << endl << coeffs3 << endl << endl;
		cout << "Filter 4" << endl << coeffs4 << endl << endl;
#endif
		filters[0] = Biquad(coeffs1);
		filters[1] = Biquad(coeffs2);
		filters[2] = Biquad(coeffs3);
		filters[3] = Biquad(coeffs4);
		filter_bank[k] = filters;
	}
	return filter_bank;
}

Eigen::VectorXd CochlearFilterbank::ERBspace(double low_freq, double high_freq, int num_channels)
{
	Eigen::VectorXd cf_array(num_channels);
	double aux = EAR_Q * MIN_BW;
	for (int i=1; i <= num_channels; i++)
	{
		cf_array[i-1] = -(aux) + exp((i)*(-log(high_freq + aux) + log(low_freq + aux))/num_channels) * (high_freq + aux);
	}
	return cf_array;
}

vector<Eigen::VectorXd> CochlearFilterbank::process(const Eigen::VectorXd &input, unsigned int n)
{
	vector<Eigen::VectorXd> out;
	Eigen::VectorXd y;
#ifdef MULTITHREAD
	vector<future<Eigen::VectorXd>> futures;
#endif
	for (unsigned int ch=0; ch < filters.size(); ch++)
	{
#ifdef MULTITHREAD
        futures.push_back(async(launch::async, &CochlearFilterbank::process_channel_sample, this, input, n, ch));
#else
		y = process_channel(input, n, ch);
		out.push_back(y);
#endif
	}
#ifdef MULTITHREAD
	for (unsigned int ch=0; ch < filters.size(); ch++)
	{
		y = futures[ch].get();
		out.push_back(y);
	}
#endif
	return out;
}

Eigen::VectorXd CochlearFilterbank::process_channel(const Eigen::VectorXd &input, int n, int ch)
{
		Eigen::VectorXd y1(input.rows());
		Eigen::VectorXd y2(input.rows());
		Eigen::VectorXd y3(input.rows());
		Eigen::VectorXd y4(input.rows());

		y1 = filters[ch][0].process(input, n);
		y2 = filters[ch][1].process(y1, n);
		y3 = filters[ch][2].process(y2, n);
		y4 = filters[ch][3].process(y3, n);
		return y4;
}

Eigen::VectorXd CochlearFilterbank::process_channel_sample(const Eigen::VectorXd &input, unsigned int n, unsigned int ch)
{
	Eigen::VectorXd y(Eigen::VectorXd::Zero(input.rows()));

	Eigen::VectorXd x1(1);
	Eigen::VectorXd y1(1);
	Eigen::VectorXd y2(1);
	Eigen::VectorXd y3(1);
	Eigen::VectorXd y4(1);

	for (unsigned int i=0; i < n; i++)
	{
		x1[0] = input[i];
		y1 = filters[ch][0].process(x1, 1);
		y2 = filters[ch][1].process(y1, 1);
		y3 = filters[ch][2].process(y2, 1);
		y4 = filters[ch][3].process(y3, 1);
		y[i] = y4[0];
	}
	return y;
}
