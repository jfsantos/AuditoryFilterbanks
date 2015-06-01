export CXXFLAGS="-g -O3 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security -fopenmp -std=c++11"
mkoctfile --mex ModulationFilterbankMex.cpp ../CochlearFilterbank/ModulationFilterBank.cpp ../CochlearFilterbank/Biquad.cpp -I../CochlearFilterbank/ -I/usr/include/eigen3/Eigen/ -lpthread -L../CochlearFilterbank/
