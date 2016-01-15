#pragma once

#include <eigen3/Eigen/Eigen>

class Biquad {
 public:
  Biquad(void);
  Biquad(Eigen::Matrix<double, 5, 1> coeffs);
  ~Biquad(void);
  // Methods
  Eigen::VectorXd process(const Eigen::VectorXd &x, int n);
  Eigen::Matrix<double, 5, 1> coeffs();

 private:
  double a0, a1, a2;
  double b1, b2;
  double z0, z1;
};
