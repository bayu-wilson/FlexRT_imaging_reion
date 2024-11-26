#include <math.h>
#include <fstream>
#include "healpix_funcs.cc"

//cosmic time variables
const int N_cos = 100000;
double a_cos1[N_cos];
double t_cos1[N_cos];
double a_cos2[N_cos];
double t_cos2[N_cos];

//Hubble parameter at scale factor a
float H(float a)  {
	return hh*100./mpc_to_km*pow(Omega_m/pow(a, 3.) + Omega_l,0.5);
}

//Critical density of the universe
float rho_crit(float a)  {
	return 3.*pow(H(a),2.)/8./pi/G_grav;
}

//Radius of a halo, as defined by the mass enclosed within the radius at which 
//the mean density within the halo is 200x the critical density of the universe
double R200(double M200, double z)  {
	double a = 1./(1.+z);
	return pow(3./800./pi*M200/hh*Msun/rho_crit((float) a), 1./3.); //returns in cgs
}

void import_cosmic_time(void)  {
	const char cosmo_file[] = "input_files/cosmo_table";
	ifstream file (cosmo_file, ios::in | ios::binary);
	
	printf("importing cosmic time\n");
	
	for (int i = 0; i < N_cos; i++)  {
		file.read((char*)&a_cos1[i], sizeof(double));
		file.read((char*)&t_cos1[i], sizeof(double));
	}
}

//cosmic time as a function of scale factor
float cosmic_time(float a)  {
	
	int ind = (int) floor((a - a_cos1[0])/(a_cos1[N_cos-1] - a_cos1[0])*(N_cos-1));
	
	float time = (float) (t_cos1[ind] + (t_cos1[ind+1] - t_cos1[ind])/(a_cos1[ind+1] - a_cos1[ind])*(a - a_cos1[ind]));
	
	return time;
}

void make_scale_factor_table(void)  {

	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < N_cos; i++)  {
			int ind;
			t_cos2[i] = ((float) i/(N_cos-1))*(t_cos1[N_cos-1] - t_cos1[0]) + t_cos1[0];
			for (int j = 0; j < N_cos; j++)  {
					if ( (t_cos2[i] >= t_cos1[j]) & (t_cos2[i] <= t_cos1[j+1]) )  {
							ind = j;
							break;
					}
			}
			a_cos2[i] = (a_cos1[ind] + (a_cos1[ind+1] - a_cos1[ind])/(t_cos1[ind+1] - t_cos1[ind])*(t_cos2[i] - t_cos1[ind]));
	}
	}
}

//get scale factor from cosmic time
float scale_factor(float t)  {
	
	int ind = (int) floor((t - t_cos2[0])/(t_cos2[N_cos-1] - t_cos2[0])*(N_cos-1));
	
	float a = (float) (a_cos2[ind] + (a_cos2[ind+1] - a_cos2[ind])/(t_cos2[ind+1] - t_cos2[ind])*(t - t_cos2[ind]));
	
	return a;
}

//cosmic time as a function of scale factor
float cosmic_time_integrate(float a)  {
	int n = 10000;
	double tc; 
	float sf[n];
	float integrand[n];
	
	for (int i = 0; i < n; i++)  {
		float ii = (float) i;
		sf[i]     = ii / (n - 1) * a + 1e-10;
		integrand[i] = 1./sf[i]/H(sf[i]);
	}
	
	tc = trapz_int(integrand, sf, n);
	
	return tc;
}

/*using the Einstein desitter space approximation for high redshift*/
//ADDED 05/27/22 - Fahad's time/redshift functions
double toTime(double redshift){

    const double H0 = 3.24086e-18; //in 1/s.  Does not include h.

    double atime =  1.0/(1.0+redshift);
    return 1.0/(1.5*H0*hh*sqrt(1.0-Omega_m))*asinh(sqrt((1.0-Omega_m)/Omega_m)*sqrt(atime)*atime);
}

double toRedshift(double time){

    const double H0 = 3.24086e-18; //in 1/s.  Does not include h.
    double atime = pow(sqrt(Omega_m/(1.0-Omega_m))*sinh(1.5*H0*hh*sqrt(1.0-Omega_m)*time), 2.0/3.0);
    return (1.0/atime)-1.0;
}

