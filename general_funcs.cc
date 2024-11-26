//General functions for coding various physical processes not related to the rt algorithm itself
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <algorithm>
#include <random>
#include <time.h> 
#include <vector>
#include "global_constants.h"
using namespace std; 

//ADDED 03/07 min/max functions for short ints
//min of two ints
int mins(short int a, short int b)  {
	if (a >= b)  {
		return b;
	}
	else  {
		return a;
	}
}

//max of two ints
int maxs(short int a, short int b)  {
	if (b >= a)  {
		return b;
	}
	else  {
		return a;
	}
}

//min of two ints
int min(int a, int b)  {
	if (a >= b)  {
		return b;
	}
	else  {
		return a;
	}
}

//max of two ints
int max(int a, int b)  {
	if (b >= a)  {
		return b;
	}
	else  {
		return a;
	}
}

//Absolute value of a float
float absd(float x)  {

	if (x >= 0.)  {
		return x;
	}
	else  {
		return -1.*x;
	}
}

//min of two floats
float mind(float a, float b)  {
	if (a >= b)  {
		return b;
	}
	else  {
		return a;
	}
}

//max of two floats
float maxd(float a, float b)  {
	if (b >= a)  {
		return b;
	}
	else  {
		return a;
	}
}

//min of two doubles
double mindd(double a, double b)  {
	if (a >= b)  {
		return b;
	}
	else  {
		return a;
	}
}

//max of two doubles
double maxdd(double a, double b)  {
	if (b >= a)  {
		return b;
	}
	else  {
		return a;
	}
}

//sign of an int
int sign(int a)  {
	if (a > 0)  {
		return 1;
	}
	else if (a < 0)  {
		return -1;
	}
	else  {
		return 0;
	}
}

//sign of a float
float signd(float a)  {
	if (a > 0.)  {
		return 1.;
	}
	else if (a < 0.)  {
		return -1.;
	}
	else  {
		return 0.;
	}
}

//sign of a float
double signdd(double a)  {
	if (a > 0.)  {
		return 1.;
	}
	else if (a < 0.)  {
		return -1.;
	}
	else  {
		return 0.;
	}
}

//dot product of two vectors
float dot(vec3 v1, vec3 v2)  {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

//length of a vector
float length(vec3 v1)  {
	return pow(dot(v1, v1), 0.5);
}

//cross product of two vectors
vec3 cross(vec3 v1, vec3 v2)  {
	vec3 v_out;
	v_out.x = v1.y*v2.z - v1.z*v2.y;
	v_out.y = v1.z*v2.x - v1.x*v2.z;
	v_out.z = v1.x*v2.y - v1.y*v2.x;
	return v_out;
}

//initialize random number generator
void init_rand()  {
	srand(time(NULL)*time(NULL)/rand());
	printf("initializing random number generator\n");
	printf("time: %le\n", (double) (time(NULL)));
}

//random integer on the interval [0, a)
unsigned short int siRand(int a)  {
   int v1 = rand() % a + 0;
   return (unsigned short int) v1;
}

//random integer on the interval [0, a)
int iRand(int a)  {
   int v1 = rand() % a + 0;
   return v1;
}

//Random float on the interval (a, b)
float fRand(float a, float b)  {
	int v1 = rand() % 1001 + 0;
	float f1 = ((float) v1)/1000.;
	return a + f1*(b - a);
}

//compute a mod b of two ints
int mod(int a, int b) {
	while (a < 0) {
		a += b;
	}
	
	int m = a % b;
	
	return m;
}

//a mod b of two long ints
long int modlong(long int a, long int b) {
	
	if (a < 0) {
		a += b;
	}
	
	long int m = a % b;
	
	return m;
}

//a mod b for doubles
double moddouble(double a, double b)  {
	if (a < 0) {
		a += b;
	}
	
	int m = (int) a % (int) b;
	
	return (double) m;
}

//compute the remainder of a/b
double remainder(double a, double b)  {
	return a - floor(a/b)*b;
}

//Blackbody radiation function.  
double b_nu(double nu, double T)
{
	return 2*h*pow(nu,3.)/pow(cl,2.)*(1./(exp(h*nu/k_B/T) - 1.));
} 

//power law density function to test the rt algorithm.  
//Use alpha = 0 for constant density plane parallel case.  
double power_law(double r, double r0, double A, double alpha)
{
	if (r != 0)  {
		return A*pow(r/r0, alpha);
	}
	else  {
		return A;
	}
}

//compute int y(x)dx using the trapezoidal method (n = length of arrays)
float trapz_int(float y[], float x[], int n)  {
	float F  = 0;
	float dF = 0;

	for (int i = 0; i < n - 1; i++)  {
		dF  = (y[i] + y[i+1])*(x[i + 1] - x[i]);
		F  += dF/2.;
	}
	return F;
}

//compute the cumulative integral F(x) = int_x0^x y(x')dx' using the trapezoidal rule
float cum_trapz_int(float y[], float x[], int n)  {
	int i,j;
	float cumint[n];
	cumint[0] = 0.;
	for (i = 1; i < n; i++)  {
		float temp_x[i + 1];
		float temp_y[i + 1];
		for (j = 0; j < i + 1; j++)  {
			temp_x[j] = x[j];
			temp_y[j] = y[j];
		}
		cumint[i] = trapz_int(temp_y, temp_x, i + 1);
	}
	return *cumint;
}

//interpolate over domain x and range y to get y0 = y(x0).  
float interpolate(float x[], float y[], float x0, int n)  {
	float y0 = 0.;
    
	bool itp = FALSE;
	
	for (int i = 0; i < n - 1; i++)  {
		if ( ( (x[i] <= x0) && (x[i+1] >= x0) )  || ( (x[i] >= x0) && (x[i+1] <= x0) ) )  {
			y0 = y[i] + (y[i+1] - y[i])/(x[i+1] - x[i])*(x0 - x[i]);
			itp = TRUE;
			break;
		}
	}
	if (itp == FALSE)  {
		if (x[0] < x[n-1])  {
			if (x0 > x[n-1])  {
				float slope = (y[n-1] - y[n-2])/(x[n-1] - x[n-2]);
				y0 = y[n-1] + slope*(x0 - x[n-1]);
			}
			else  {
				float slope = (y[1] - y[0])/(x[1] - x[0]);
				y0 = y[0] + slope*(x0 - x[0]);
			}
		}
		else  {
			if (x0 < x[n-1])  {
				float slope = (y[n-1] - y[n-2])/(x[n-1] - x[n-2]);
				y0 = y[n-1] + slope*(x0 - x[n-1]);
			}
			else  {
				float slope = (y[1] - y[0])/(x[1] - x[0]);
				y0 = y[0] + slope*(x0 - x[0]);
			}
		}
	}
	return y0;
}

/*void cic_assignment(int i, int j, int k, float xmid, float ymid, float zmid, float f, float f_array[Nx][Ny][Nz], omp_lock_t locks[Nx][Ny][Nz])  {
	float weights[3][3][3];
	
	if (xmid >= 0.5)  {
		dxim = 0.;
		dxi  = 1.5 - xmid;
		dxip = xmid - 0.5;
	}
	else  {
		dxim = 0.5 - xmid;
		dxi  = xmid + 0.5;
		dxip = 0.;
	}
	
	if (ymid >= 0.5)  {
		dyjm = 0.;
		dyj  = 1.5 - ymid;
		dyjp = ymid - 0.5;
	}
	else  {
		dyjm = 0.5 - ymid;
		dyj  = ymid + 0.5;
		dyjp = 0.;
	}
	
	if (zmid >= 0.5)  {
		dzkm = 0.;
		dzk  = 1.5 - zmid;
		dzkp = zmid - 0.5;
	}
	else  {
		dzkm = 0.5 - zmid;
		dzk  = zmid + 0.5;
		dzkp = 0.;
	}
	
	weights[i-1][j-1][k-1] = dxim*dyjm*dzkm;
	weights[i][j-1][k-1]   = dxi*dyjm*dzkm;
	weights[i+1][j-1][k-1] = dxip*dyjm*dzkm;
	weights[i-1][j][k-1] = dxim*dyj*dzkm;
	weights[i][j][k-1]   = dxi*dyj*dzkm;
	weights[i+1][j][k-1] = dxip*dyj*dzkm;
	weights[i-1][j+1][k-1] = dxim*dyjp*dzkm;
	weights[i][j+1][k-1]   = dxi*dyjp*dzkm;
	weights[i+1][j+1][k-1] = dxip*dyjp*dzkm;
	
	weights[i-1][j-1][k] = dxim*dyjm*dzk;
	weights[i][j-1][k]   = dxi*dyjm*dzk;
	weights[i+1][j-1][k] = dxip*dyjm*dzk;
	weights[i-1][j][k] = dxim*dyj*dzk;
	weights[i][j][k]   = dxi*dyj*dzk;
	weights[i+1][j][k] = dxip*dyj*dzk;
	weights[i-1][j+1][k] = dxim*dyjp*dzk;
	weights[i][j+1][k]   = dxi*dyjp*dzk;
	weights[i+1][j+1][k] = dxip*dyjp*dzk;
	
	weights[i-1][j-1][k+1] = dxim*dyjm*dzkp;
	weights[i][j-1][k+1]   = dxi*dyjm*dzkp;
	weights[i+1][j-1][k+1] = dxip*dyjm*dzkp;
	weights[i-1][j][k+1] = dxim*dyj*dzkp;
	weights[i][j][k+1]   = dxi*dyj*dzkp;
	weights[i+1][j][k+1] = dxip*dyj*dzkp;
	weights[i-1][j+1][k+1] = dxim*dyjp*dzkp;
	weights[i][j+1][k+1]   = dxi*dyjp*dzkp;
	weights[i+1][j+1][k+1] = dxip*dyjp*dzkp;
	
	for (int m = i-1; m <= i+1; m++)  {
		for (int p = j-1; p <= j+1; p++)  {
			for (int q = k-1; q <= k+1; q++)  {
				omp_set_lock(&(locks[m][p][q]));
				f_array[m][p][q] += f*weights[m-i-1][p-j-1][q-k-1];
				omp_unset_lock(&(locks[m][p][q]));
			}
		}
	}
	
}*/