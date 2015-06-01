# AuditoryFilterbanks
Gammatone and modulation filterbanks in C++ (plus MEX files for MATLAB/Octave)

These are multithreaded implementations of gammatone and modulation filterbanks using C++11 and [Eigen](http://eigen.tuxfamily.org). This repository also includes MEX wrappers for MATLAB and Octave.

## Build instructions

Clone the repository and then do as follows:

```
cd AuditoryFilterbanks
mkdir build
cd build
cmake ../
make
```

For building the MEX files for Octave, use the `build.sh` scripts included in the wrapper folders. For building MEX files for MATLAB on Windows, you will need a recent version of Microsoft Visual Studio that supports C++11 (VS2013 should do the job). I needed to use the `Source.def` files to compile the MEX files directly from Visual Studio (you may need to do that if your MATLAB version does not support the Visual Studio version you have installed).
