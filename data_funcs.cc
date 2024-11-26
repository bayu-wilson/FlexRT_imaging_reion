#include <math.h>
#include <vector>
#include "global_variables.cc"
using namespace std;

//count the number of photons in the box
double count_photons(void)  {
	double photon_counter = 0.;
	#pragma omp parallel reduction (+ : photon_counter)
	{
	#pragma omp for
	for (long int x = 0; x < size_buffer; x++)  {
		if (ray_tags[x] == TRUE)  {
			photon_counter += (double) ray[x].nphot[0];
		}
	}
	}
	return photon_counter*phot_norm;
}

//ADDED 03/20: find the maximum value in an array
float maxarr(float arr[Nx][Ny][Nz])  {
	float out = 0.;
	
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				if (arr[i][j][k] > out)  {
					out = arr[i][j][k];
				}
			}
		}
	}
	return out;
}

//ADDED 03/20: find the minimum value in an array
float minarr(float arr[Nx][Ny][Nz])  {
	float out = 1e32;
	
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				if (arr[i][j][k] < out)  {
					out = arr[i][j][k];
				}
			}
		}
	}
	return out;
}

float frac_gtr_thresh(float f[Nx][Ny][Nz], float thresh)  {
	double tot    = 0.;
	 double weight = (double) Nx;
	
	#pragma omp parallel reduction (+ : tot)
	{
	#pragma omp for
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nx; k++)  {
				if (f[i][j][k] > thresh)  {
					tot += 1.;
				}
			}
		}
	}
	}
	float ans = (float) tot/pow(weight, 3.);
	return ans;
}

//Volume-weighted average
float calc_vol_avg(float f[Nx][Ny][Nz])  {
	 double tot    = 0.;
	 double weight = (double) Nx;
	
	#pragma omp parallel reduction (+ : tot)
	{
	#pragma omp for
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nx; k++)  {
				tot += (double) f[i][j][k];
			}
		}
	}
	}
	float ans = (float) tot/pow(weight, 3.);
	return ans;
}

//Mass-weighted average
float calc_mass_avg(float f[Nx][Ny][Nz], float den[Nx][Ny][Nz])  {
	double tot    = 0;
	double m_tot  = 0;

	#pragma omp parallel reduction (+ : tot, m_tot)
	{
	#pragma omp for
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				double m = (double) den[i][j][k];
				m_tot += m;
				tot   += (double) f[i][j][k]*m;
			}
		}
	}
	}
	float ans = (float) tot/m_tot;
	return ans;
}
