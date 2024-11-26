#include <math.h>
#include "cosmo_funcs.cc"

//cross-section parameters from Verner et. al. 1996, Table
double const sigma0_H1  = 5.475e-14;
double const sigma0_He1 = 9.492e-16;
double const sigma0_He2 = 1.369e-14;

double const Eth_H1     = 1.360e+1;
double const Eth_He1    = 2.459e+1;
double const Eth_He2    = 5.442e+1;

double const Emax_H1    = 5e+4;
double const Emax_He1   = 5e+4;
double const Emax_He2   = 5e+4;

double const E0_H1      = 4.298e-1;
double const E0_He1     = 1.361e+1;
double const E0_He2     = 1.720e0;

double const y0_H1      = 0;
double const y0_He1     = 4.434e-1;
double const y0_He2     = 0;

double const y1_H1      = 0;
double const y1_He1     = 2.136;
double const y1_He2     = 0;

double const yw_H1      = 0;
double const yw_He1     = 2.039;
double const yw_He2     = 0;

double const ya_H1      = 3.288e+1;
double const ya_He1     = 1.469;
double const ya_He2     = 3.288e+1;

double const P_H1       = 2.963;
double const P_He1      = 3.188;
double const P_He2      = 2.963;

double lambda(double T, double Eth)  {

    double K2eV = k_B/eV_cgs;
    double eV2K = 1.0/K2eV;
    double x    = 2*Eth*eV2K/T;
    return x;

}

double F_of_y(double nu, double Eth, double Emax, double E0, double y0, double y1, double yw, double ya, double P)  {

	double x, y;
	double E, F;    

	E = h_eV*nu;
	x = E/E0 - y0;
	y = pow(pow(x,2) + pow(y1,2),0.5);
	F = (pow(x-1,2) + pow(yw,2))*pow(y,0.5*P - 5.5)*pow(1 + pow(y/ya,0.5),-P);

	if (h_eV*nu >= Eth && h_eV*nu <= Emax)  {
		return F;
	}
	else  {
		return 0;
    } 
}

//photoionization cross-sections

//photo-ionization cross-section of HI from Verner et. al. 1996
double sigmapi_H1(double nu)  {
	return sigma0_H1*F_of_y(nu, Eth_H1, Emax_H1, E0_H1, y0_H1, y1_H1, yw_H1, ya_H1, P_H1);
}

//photo-ionization cross-section of HeI from Verner et. al. 1996
double sigmapi_He1(double nu)  {
	return sigma0_He1*F_of_y(nu, Eth_He1, Emax_He1, E0_He1, y0_He1, y1_He1, yw_He1, ya_He1, P_He1);
}

//photo-ionization cross-section of HeII from Verner et. al. 1996
double sigmapi_He2(double nu)  {
	return sigma0_He2*F_of_y(nu, Eth_He2, Emax_He2, E0_He2, y0_He2, y1_He2, yw_He2, ya_He2, P_He2);
}

//collisional ionization rates

//collisional ionization rate of HI from Hui & Gnedin 1996
double cic_H1(double T)  {

	double x = lambda(T, Eth_H1);
	return 21.11*pow(T,-1.5)*exp(-x/2.)*pow(x,-1.089)/(pow(1. + pow(x/0.354,0.874),1.101));
}

//collisional ionization rate of HeI from Hui & Gnedin 1996
double cic_He1(double T)  {

	double x = lambda(T, Eth_He1);
	return 32.38*pow(T,-1.5)*exp(-x/2.)*pow(x,-1.146)/(pow(1. + pow(x/0.416,0.987),1.056));
}

//collisional ionization rate of HeII from Hui & Gnedin 1996
double cic_He2(double T)  {

	double x = lambda(T, Eth_He2);
	return 19.95*pow(T,-1.5)*exp(-x/2.)*pow(x,-1.089)/(pow(1. + pow(x/0.553,0.735),1.275));
}

//collisional ionization cooling rates

//collisional ionization cooling rate of HI from Hui & Gnedin 1996
double cicr_H1(double T)  {
    return Eth_H1*ev_to_erg*cic_H1(T);
}

//collisional ionization cooling rate of HeI from Hui & Gnedin 1996
double cicr_He1(double T)  {
	return Eth_He1*ev_to_erg*cic_He1(T);
}

//collisional ionization cooling rate of HeII from Hui & Gnedin 1996
double cicr_He2(double T)  {
	return Eth_He2*ev_to_erg*cic_He2(T);
}

//recombination rates

//Case A recombination coefficient of HII from Hui & Gnedin 1996
double alphaA_H2(double T)  {

	double x = lambda(T, Eth_H1);
	return 1.269e-13*pow(x,1.503)/pow((1.0+pow((x/0.522),0.470)),1.923);
}

//Case B recombination coefficient of HII from Hui & Gnedin 1996
double alphaB_H2(double T)  {
	double x = lambda(T, Eth_H1);
	return 2.753e-14*pow(x,1.500)/pow((1.0+pow((x/2.740),0.407)),2.2420);
}

//Case A recombination coefficient of HeII from Hui & Gnedin 1996
double alphaA_He2(double T)  {
	double x = lambda(T, Eth_He1);
	return 3.0e-14*pow(x,0.654);
}

//Case B recombination coefficient of HeII from Hui & Gnedin 1996
double alphaB_He2(double T)  {

	double x    = lambda(T, Eth_He1);
	return 1.26e-14*pow(x,0.750);
}

//Case A recombination coefficient of HeIII from Hui & Gnedin 1996
double alphaA_He3(double T)  {

	double x    = lambda(T, Eth_He2);
	return 2.*1.269e-13*pow(x,1.503)/pow((1.0+pow((x/0.522),0.470)),1.923);
}

//Case B recombination coefficient of HeIII from Hui & Gnedin 1996
double alphaB_He3(double T)  {

	double x    = lambda(T, Eth_He2);
	return 2*2.753e-14*pow(x,1.500)/pow((1.0+pow((x/2.740),0.407)),2.2420);
}

//Dielectric recombination coefficient of He II from Hui & Gnedin 1996
double Dalpha_He2(double T)  {
	double x    = lambda(T, Eth_He2);
	return 1.90e-3 * pow(T, -1.5) * exp(-0.75*x/2)*(1 + 0.3*exp(-0.15*x/2));
}

//recombination cooling rates

//Case A recombination cooling rate for HII from Hui & Gnedin 1996
double rcrA_H2(double T)  {

	double x    = lambda(T, Eth_H1);
	return 1.778e-29*T*pow(x,1.965)/(pow(1. + pow(x/0.541,0.502),2.697));
}

//Case B recombination cooling rate for HII from Hui & Gnedin 1996
double rcrB_H2(double T)  {
    
	double x = lambda(T, Eth_H1);
	return 3.435e-30*T*pow(x,1.970)/(pow(1. + pow(x/2.250,0.376),3.720));
}

//Case A recombination cooling rate for HeII from Hui & Gnedin 1996
double rcrA_He2(double T)  {
	return k_B*T*alphaA_He2(T);
}

//Case B recombination cooling rate for HeII from Hui & Gnedin 1996
double rcrB_He2(double T)  {
	return k_B*T*alphaB_He2(T);
}

//Case A recombination cooling rate for HeIII from Hui & Gnedin 1996
double rcrA_He3(double T)  {
	double x    = lambda(T, Eth_He2);
	return 1.778e-29*T*pow(x,1.965)/(pow(1. + pow(x/0.541,0.502),2.697));
}

//Case B recombination cooling rate for HeIII from Hui & Gnedin 1996
double rcrB_He3(double T)  {
	double x = lambda(T, Eth_He2);
	return 3.435e-30*T*pow(x,1.970)/(pow(1. + pow(x/2.250,0.376),3.720));
}

//dielectric recombination cooling (He II).  
double drcr_He2(double T)  {
	return 0.75*Eth_He2*ev_to_erg*Dalpha_He2(T);
}

//collisional excitation cooling rates from Cen 1992
double coll_ex_rate_H1(double T)  {
	return 7.5e-19*pow((1 + pow(T/1e5, 0.5)), -1)*exp(-118348./T);
}

double coll_ex_rate_He1(double T)  {
	return 9.10e-27*pow(T, -0.1687)*pow((1 + pow(T/1e5, 0.5)), -1)*exp(-13179./T);
}

double coll_ex_rate_He2(double T)  {
	return 5.54e-17*pow(T, -0.397)*pow((1 + pow(T/1e5, 0.5)), -1)*exp(-473638./T);
}

double free_free_rate(double T)  {
	return 1.42e-27*pow(T, 0.5)*1.30;
}

//compton heating/cooling rate from Hui & Gnedin 1996 (CORRECTED 07/09/21)
double compton_rate(double z)  {
	double chi = 0.08;
	//return 6.35e-41*Omega_b*pow(hh, 2.)*(1. + chi)/(1.+4.*chi)*pow(1 + z, 7.);
	return 5.6072e-36*pow(1 + z, 4.);
}
//H I photoionization heating
double phothtHI(double alpha, double T, double z, double ne)
{ 
  return Eth_H1/h_eV*h*((alpha+2.8)/(1.8+alpha)*( 1.- pow(4.,-1.8-alpha ) )/( 1. - pow(4.,-2.8-alpha) ) - 1.)*alphaA_H2(T)*ne; 
}

//He I photoionization heating
double phothtHeI(double alpha, double T , double z, double ne)
{
  return Eth_He1/h_eV*h*((alpha+1.7)/(0.7+alpha)*( 1. - pow(2.2,-0.7-alpha) )/( 1. - pow(2.2,-1.7-alpha)) - 1.)*alphaA_He2(T)*ne;
}

double Photoheat_HIHeI(double alpha, double T, double z, double ne, double nH, double nHe)
{
  return phothtHI(alpha,T,z, ne)*nH + phothtHeI(alpha,T,z, ne)*nHe;
}