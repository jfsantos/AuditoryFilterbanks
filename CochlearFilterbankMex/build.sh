export CXXFLAGS="-g -O3 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security -fopenmp -std=c++11"
mkoctfile --mex CochlearFilterbankMex.cpp ../CochlearFilterbank/CochlearFilterbank.cpp ../CochlearFilterbank/Biquad.cpp -I../CochlearFilterbank/ -I/usr/include/eigen3/Eigen/ -lpthread -L../CochlearFilterbank/
