#include "CochlearFilterbank.h"
#include "ModulationFilterBank.h"
#include <iostream>
#include <vector>
#include <eigen3/Eigen/Eigen>
#include "timer.h"

using namespace std;

#define N_SAMPLES 8000*60*2

int main()
{
	CochlearFilterbank cf = CochlearFilterbank(8000.0, 23, 150.0);

	// Testing ERBspace
	Eigen::VectorXd cf_array = cf.ERBspace(150, 8000/2, 23);
	cout << cf_array << endl;

	// Testing filters
	Eigen::VectorXd x(Eigen::VectorXd::Zero(N_SAMPLES));
	x[0] = 1.0; // unit impulse

	// Running filters
	vector<Eigen::VectorXd> y;
    Timer tt;
    y = cf.process(x, N_SAMPLES);
	double t = tt.elapsed();
    cout << "Total processing time for the cochlear filterbank: " << t << " s." << endl;
#ifdef DEBUG
	for (int k=0; k < 23; k++)
	{
		char fname[15];
		sprintf(fname,"channel%d.txt",k+1);
		FILE* f = fopen(fname,"w");
		for (int n=0; n < N_SAMPLES; n++)
			fprintf(f,"%e\n",y[k](n));
		fclose(f);
	}
#endif
    cout << "Tested CochlearFilterBank successfully." << endl;

// 	// Testing the ModulationFilterBank
 	Eigen::VectorXd central_freqs(8);
 	central_freqs = ModulationFilterBank::compute_modulation_cfs(4.0, 128.0, 8);
 	cout << central_freqs << endl;
 	ModulationFilterBank mf = ModulationFilterBank(16000, 8, central_freqs, 2);
 	// Running filters
 	Eigen::MatrixXd ym;
 	tt.reset();
 	ym = mf.process(x, N_SAMPLES);
    t = tt.elapsed();
    cout << "Total processing time for the modulation filterbank: " << t << " s." << endl;
#ifdef DEBUG
 	for (int k=0; k < 23; k++)
 	{
 		char fname[15];
 		sprintf(fname,"mf%d.txt",k+1);
 		FILE* f = fopen(fname,"w");
 		for (int n=0; n < N_SAMPLES; n++)
 			fprintf(f,"%e\n",ym(n,k));
 		fclose(f);
 	}
#endif
        cout << "Tested ModulationFilterBank successfully." << endl;
}
