#include "Biquad.h"

Biquad::Biquad(void)
{
	a0 = 0.0;
	a1 = 0.0;
	a2 = 0.0;
	b1 = 0.0;
	b2 = 0.0;
}


Biquad::Biquad(Eigen::Matrix<double, 5, 1> coeffs)
{
	a0 = coeffs[0];
	a1 = coeffs[1];
	a2 = coeffs[2];
	b1 = coeffs[3];
	b2 = coeffs[4];
}

Eigen::Matrix<double, 5, 1> Biquad::coeffs()
{
	Eigen::Matrix<double, 5, 1> coeffs;
	coeffs << a0, a1, a2, b1, b2;
	return coeffs;
}

Biquad::~Biquad(void)
{
}

// FIXME: store internal state into Biquad to allow streaming
Eigen::VectorXd Biquad::process(const Eigen::VectorXd &x, int n)
{
	Eigen::VectorXd y(Eigen::VectorXd::Zero(x.rows()));
	double xi, yi, z0 = 0, z1 = 0;
	for (int i = 0; i < n ; ++i)
	{
		xi = x[i];
		yi = a0*xi + z0;
		z0 = a1 * xi - b1 * yi + z1;
		z1 = a2 * xi - b2 * yi;
		y[i] = yi;
		//y[i] = a0 * x[i] + a1 * x[i-1] + a2 * x[i-2]
        //               - b1 * y[i-1] - b2 * y[i-2];
	}
	return y;
}
