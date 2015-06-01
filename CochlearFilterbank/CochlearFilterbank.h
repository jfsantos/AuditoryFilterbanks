#pragma once

#include <cmath>
#include <complex>
#include <vector>
#include <eigen3/Eigen/Eigen>
#include <future>
#include <vector>
#include "Biquad.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class CochlearFilterbank
{
public:
	// Constructors and destructors
	CochlearFilterbank(double fs, int num_channels, double low_freq);
	~CochlearFilterbank(void);
	// Methods
	Eigen::MatrixXd process(const Eigen::VectorXd &input, int n);
	Eigen::VectorXd process_channel(const Eigen::VectorXd &input, int n, int ch);
	static Eigen::VectorXd ERBspace(double low_freq, double high_freq, int num_channels);
	// Constants are defined on source file
private:
	std::vector<std::vector<Biquad> > filters;
	std::vector<std::vector<Biquad> > makeERBFilters(double fs, int num_channels, double low_freq);
};
