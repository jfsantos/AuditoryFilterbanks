#pragma once

#include <cmath>
#include <complex>
#include <vector>
#include <eigen3/Eigen/Eigen>
#include "Biquad.h"
#include <future>

class ModulationFilterBank {
 public:
  ModulationFilterBank(double fs, int num_channels, Eigen::VectorXd mf,
                       double q);
  ~ModulationFilterBank(void);
  Eigen::MatrixXd process(const Eigen::VectorXd &input, int n);
  Eigen::VectorXd process_channel(const Eigen::VectorXd &input, int n, int ch);
  static Eigen::VectorXd compute_modulation_cfs(double min_cf, double max_cf,
                                                int n);
  // Constants are defined on source file
 private:
  std::vector<Biquad> filters;
  std::vector<Biquad> makeModulationFilters(double fs, int num_channels,
                                            Eigen::VectorXd mf, double q);
};
