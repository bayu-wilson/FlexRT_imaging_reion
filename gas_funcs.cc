#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <string>
#include "subgrid.cc"
using namespace std;

//MOVED to gas file on 03/04 for use in the energy calculation
//Fraction of the total photons in each frequency bin.  This has to be set by hand by the user.  

void read_spectrum_table(void)  {
	double temp_double;
	FILE *file = NULL;

	const char * spectrum_table_file = spectrum_table.c_str();

	printf("opening spectrum interpolation table\n");
	file = fopen(spectrum_table_file, "r");

	if (file == NULL)  {
		printf("Spectrum table not found\n");
	}
	else  {
		printf("Reading spectrum table\n");
	}

	//get rid of the frequency header
	for (int i = 0; i < Nfreq-1; i++)  {
		fscanf(file, "%le", &temp_double);
	}
	fscanf(file, "%le\n", &temp_double);

	for (int a = 0; a < Nalph; a++)  {
		fscanf(file, "%le", &alpha[a]);
		for (int nu = 0; nu < Nfreq-1; nu++)  {
			fscanf(file, "%le", &spectrum_input[a][nu]);
		}
		fscanf(file, "%le\n", &spectrum_input[a][Nfreq-1]);
	}
}

//MOVED to gas file on 03/04 for use in the energy calculation
//Fraction of the total photons in each frequency bin.  This has to be set by hand by the user.
void set_spectrum(void)  {
	spectrum[0] = 1.0; //monochromatic case

	//interpolate spectrum from table if doing multi-frequency RT
	if (multi_freq == TRUE)  {
		double alpha_spec = 1.5;

		//locate alpha indices for interpolation
		int alpha_ind_m = floor((alpha_spec - alpha[0])/(alpha[Nalph-1] - alpha[0])*(Nalph-1));
		alpha_ind_m = min(alpha_ind_m, Nalph-2);
		int alpha_ind_p = alpha_ind_m + 1;

		//interpolate the spectrum in each bin over alpha
		float spec_sum = 0.;
		for (int nu = 0; nu < Nfreq; nu++)  {
			spectrum[nu] = (float) (spectrum_input[alpha_ind_m][nu] + (spectrum_input[alpha_ind_p][nu] - spectrum_input[alpha_ind_m][nu])/(alpha[alpha_ind_p] - alpha[alpha_ind_m])*(alpha_spec - alpha[alpha_ind_m]));
			spec_sum += spectrum[nu];
		}

		//make sure the spectrum is normalized
		for (int nu = 0; nu < Nfreq; nu++)  {
			spectrum[nu] /= spec_sum;
			printf("spectrum[%d] = %le\n", nu, spectrum[nu]);
		}

	}
}

//Compute the photo-ionization rate from the radiation energy density, if not using the C2Ray scheme.  
//If collisonal ionization is on, add this to the total.  
void update_gamma(void)  {
	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				//photoionization
				gamma_H1_tot[i][j][k]  = 0.;
				gamma_He1_tot[i][j][k] = 0.;
				gamma_He2_tot[i][j][k] = 0.;
				for (int nu = 0; nu < Nfreq; nu++)  {
					gamma_H1_tot[i][j][k]  += sigmapi_H1(nu_phot[nu])*clight[i][j][k]*u_nu[i][j][k][nu]/h/nu_phot[nu];
					gamma_He1_tot[i][j][k] += sigmapi_He1(nu_phot[nu])*clight[i][j][k]*u_nu[i][j][k][nu]/h/nu_phot[nu];
					gamma_He2_tot[i][j][k] += sigmapi_He2(nu_phot[nu])*clight[i][j][k]*u_nu[i][j][k][nu]/h/nu_phot[nu];
				}
				//collisional ionization counts as an effective photoionization rate with coeff
				//cic x ne
				if ( (coll_ion == TRUE) && (temp[i][j][k] >= 1e4) && (temp[i][j][k] <= 1e9) )  {
					gamma_H1_tot[i][j][k]  += cic_H1(temp[i][j][k])*ne[i][j][k];
					gamma_He1_tot[i][j][k] += cic_He1(temp[i][j][k])*ne[i][j][k];
					gamma_He2_tot[i][j][k] += cic_He2(temp[i][j][k])*ne[i][j][k];
				}
			}
		}
	}
	}
}

//set the initial temperature of the gas
void set_Treion_default(int i, int j, int k)  {
	float fion_new = (nH[i][j][k] - nH1[i][j][k])/nH[i][j][k];
	float fion_old = (nH[i][j][k] - nH1_prev[i][j][k])/nH[i][j][k];
	float deltax = dx[i]*(fion_new - fion_old);
	float vIF = deltax/dt;
	float Cn[5] = {9.5432, -1.6441e-2, -9.8010e-3, 4.1664e-3, -2.2710e-4};
	float lnTreion = 0.;
	for (int t = 0; t < 5; t++)  {
		lnTreion += Cn[t]*pow(log(vIF/1e5), t); //Treion from Daloisio+18, in K
	}
	//Modified 12/16/21: only set Treion during I-front crossing (up to 95% ionized)
	if ((fion_new > fion_old) && (fion_new < 0.98))  { 
		temp[i][j][k] = (fion_old*temp[i][j][k] + (fion_new - fion_old)*exp(lnTreion))/fion_new;
	}
	else  {
		temp[i][j][k] = temp[i][j][k];
	}
	
	//Added 08/31/22 - keep track of IF speed and incident flux
	IF_speed[i][j][k] = maxd(0., vIF); //0 if negative IF speed
	for (int nu=0;nu<Nfreq;nu++){
		Inc_Flux[i][j][k] = (1. + nhe_to_nh)*nH[i][j][k]*maxd(0., vIF); //use approximate equation assuming non-relativistic IF.
		//Inc_Flux[i][j][k] += ngam_nu[i][j][k][nu]; 
	}
		//(1. + nhe_to_nh)*nH[i][j][k]*maxd(0., vIF); //use approximate equation assuming non-relativistic IF.  
}

//MODIFIED on 03/04 to contain heating/cooling terms for sub-grid model temperature evolution
//Compute the heating and cooling rates, if tracking temperature evolution.  Might wanna update this to use 
//gamma to compute photo-heating rates
void update_heat_cool(int i, int j, int k)  {
	heat_rate[i][j][k] = 0;
	cool_rate[i][j][k] = 0;

	//float local_spectrum[Nfreq];
	//float ngam_sum = 0.; // 23/1/6 BAYU: moved outside if statement so it stays in this scope
	//if (multi_freq == TRUE)  { //heating correction for spectral hardening
	//	//float ngam_sum = 0.;
	//	for (int nu = 0; nu < Nfreq; nu++)  {
	//		ngam_sum += ngam_nu[i][j][k][nu];
	//	}
	//}
	
	//photoheating
	if (subgrid == FALSE)  {
		if (f_H2[i][j][k] > 0.98)  { //turn off photo-heating during IF passage
			for (int nu = 0; nu < Nfreq; nu++)  {
				//if (multi_freq == TRUE)  {
				//	local_spectrum[nu] = ngam_nu[i][j][k][nu]/ngam_sum;
				//}
				//else  {
				//	local_spectrum[nu] = spectrum[nu];
				//}
				heat_rate[i][j][k] += avg_nh1[i][j][k]*gamma_H1_nu[i][j][k][nu]*h*(h_eV*nu_phot[nu] - Eth_H1)/h_eV;
				heat_rate[i][j][k] += avg_nhe1[i][j][k]*gamma_He1_nu[i][j][k][nu]*h*(h_eV*nu_phot[nu] - Eth_He1)/h_eV;
	 		}
		}
		else  {
			heat_rate[i][j][k] = 0.;
		}
	}
	else  { //MODIFIED 03/15 to use Anson's analyical expression for photoheating
		if ( (gamma_H1_tot[i][j][k] > 3e-14) )  {
			heat_rate[i][j][k] = (float) Photoheat_HIHeI(1.5, (double) maxd(temp_min, temp[i][j][k]), (double) zz, (double) (nH[i][j][k]+nHe[i][j][k]), (double) nH[i][j][k], (double) nHe[i][j][k]);
		}
		else  {
			heat_rate[i][j][k] = 0.;
		}
	}
	
	//collisional ionization cooling
	if ( (coll_ion == TRUE) && (temp[i][j][k] >= 1e4) && (temp[i][j][k] <= 1e9) )  {
		//collisional ionization cooling contributions from each ion
		if (subgrid == FALSE)  {
			cool_rate[i][j][k] += cicr_H1(maxd(temp_min, temp[i][j][k]))*ne[i][j][k]*nH1[i][j][k];
			cool_rate[i][j][k] += cicr_He1(maxd(temp_min, temp[i][j][k]))*ne[i][j][k]*nHe1[i][j][k];
		}
		//we will explicitly neglect collisional ionization cooling for now, as we have no reliable way to estimate what nHI should be
	}

	//case A or B recombination cooling
	if (recomb_cool == TRUE)  {
		if (strcmp(rec_case, "A") == 0)  {
			if (subgrid == FALSE)  {
				cool_rate[i][j][k] += rcrA_H2(maxd(temp_min, temp[i][j][k]))*ne[i][j][k]*nH2[i][j][k];
				cool_rate[i][j][k] += rcrA_He2(maxd(temp_min, temp[i][j][k]))*ne[i][j][k]*nHe2[i][j][k];
			}
			else  {
				cool_rate[i][j][k] += rcrA_H2(maxd(temp_min, temp[i][j][k]))*(nH[i][j][k] + nHe[i][j][k])*nH[i][j][k];
				cool_rate[i][j][k] += rcrA_He2(maxd(temp_min, temp[i][j][k]))*(nH[i][j][k] + nHe[i][j][k])*nHe[i][j][k];
			}
		}
		else if (strcmp(rec_case, "B") == 0) {
			if (subgrid == FALSE)  {
				cool_rate[i][j][k] += rcrB_H2(maxd(temp_min, temp[i][j][k]))*ne[i][j][k]*nH2[i][j][k];
				cool_rate[i][j][k] += rcrB_He2(maxd(temp_min, temp[i][j][k]))*ne[i][j][k]*nHe2[i][j][k];
			}
			else  {
				cool_rate[i][j][k] += rcrB_H2(maxd(temp_min, temp[i][j][k]))*(nH[i][j][k] + nHe[i][j][k])*nH[i][j][k];
				cool_rate[i][j][k] += rcrB_He2(maxd(temp_min, temp[i][j][k]))*(nH[i][j][k] + nHe[i][j][k])*nHe[i][j][k];
			}
		}
		else  {
			printf("Warning: Ignoring recombinations\n");
		}
	}
	//collisional excitation cooling contributions from HI, HeI, and HeII
	if (coll_exc_cool == TRUE)  {
		if (subgrid == FALSE)  {
			cool_rate[i][j][k] += coll_ex_rate_H1(maxd(temp_min, temp[i][j][k]))*ne[i][j][k]*nH1[i][j][k];
			cool_rate[i][j][k] += coll_ex_rate_He1(maxd(temp_min, temp[i][j][k]))*ne[i][j][k]*nHe1[i][j][k];
		}
		//Also neglect collisional excitation cooling in the subgrid model for now, since we do not know how to get nHI
	}
	
	if (free_free == TRUE)  {
		if (subgrid == FALSE)  {
			float nff = nH2[i][j][k] + nHe2[i][j][k] + 4.*nHe3[i][j][k];
			cool_rate[i][j][k] += free_free_rate(maxd(temp_min, temp[i][j][k]))*ne[i][j][k]*nff;
		}
		else  { //we can approximately include free-free cooling
			float nff = fion[i][j][k]*(nH[i][j][k] + nHe[i][j][k]);
			cool_rate[i][j][k] += free_free_rate(maxd(temp_min, temp[i][j][k]))*pow(nff, 2.);
		}
	}
	
	//CORRECTED on 03/04
	//CORRECTED on 07/09/21
	//compton heating/cooling.  Counts as a positive cooling rate if T > T_cmb.  
	if (compton == TRUE)  {
		//cool_rate[i][j][k] -= compton_rate(zz)*(2.726*(1 + zz) - maxd(temp_min, temp[i][j][k]));
		if (subgrid == FALSE)  {
			cool_rate[i][j][k] -= compton_rate(zz)*ne[i][j][k]*(2.726*(1 + zz) - maxd(temp_min, temp[i][j][k]));
		}
		else  {
			cool_rate[i][j][k] -= compton_rate(zz)*(nH[i][j][k] + nHe[i][j][k])*(2.726*(1 + zz) - maxd(temp_min, temp[i][j][k]));
		}
	}
}

//1st order backwards differece scheme for solving the ionization balance equations, if the C2Ray scheme is not
//turned on.  
void solve_ion(int i, int j, int k)  {
	
	//solve the ionizing balance equation for H using a 1st order implicit backwards-difference scheme 
	if (nH1_prev[i][j][k] >= nH2_prev[i][j][k])  {
		//Solve for nH2 when there is less ionized gas
		nH2[i][j][k]	= (nH2_prev[i][j][k] + gamma_H1_tot[i][j][k]*nH[i][j][k]*dt)
						/(1. + gamma_H1_tot[i][j][k]*dt + clump[i][j][k]*recomb_H2[i][j][k]*ne[i][j][k]*dt);
		nH1[i][j][k] 	= nH[i][j][k] - nH2[i][j][k];
	}
	else  {
		//solve for nH1 when there is less neutral gas
		nH1[i][j][k]	= (nH1_prev[i][j][k] + clump[i][j][k]*recomb_H2[i][j][k]*ne[i][j][k]*nH[i][j][k]*dt)
						/(1. + gamma_H1_tot[i][j][k]*dt + clump[i][j][k]*recomb_H2[i][j][k]*ne[i][j][k]*dt);
		nH2[i][j][k]  = nH[i][j][k] - nH1[i][j][k];
	}
	
	//solve the ionizing balance equations for He.  Solve the backwards difference equation for the
	//species that have the smaller abundances, and solve closing equation for the remaining species.  
	if ( (nHe1_prev[i][j][k] <= nHe3_prev[i][j][k]) && (nHe2_prev[i][j][k] <= nHe3_prev[i][j][k]) )  {
		
		nHe1[i][j][k] = (nHe1_prev[i][j][k] + clump[i][j][k]*recomb_He2[i][j][k]*ne[i][j][k]*(nHe[i][j][k] - nHe3[i][j][k])*dt)
				      /(1. + gamma_He1_tot[i][j][k]*dt + clump[i][j][k]*recomb_He2[i][j][k]*ne[i][j][k]*dt);
		nHe2[i][j][k] = (nHe2_prev[i][j][k] + gamma_He1_tot[i][j][k]*(nHe[i][j][k] - nHe3[i][j][k])*dt
					  + clump[i][j][k]*recomb_He3[i][j][k]*ne[i][j][k]*(nHe[i] - nHe1[i])*dt)
				      /(1. + (gamma_He1_tot[i][j][k] + gamma_He2_tot[i][j][k])*dt
					  + clump[i][j][k]*(recomb_He2[i][j][k] + recomb_He3[i][j][k])*ne[i][j][k]*dt);
		nHe3[i][j][k] = nHe[i][j][k] - nHe1[i][j][k] - nHe2[i][j][k];
	}
	else if( (nHe1_prev[i][j][k] <= nHe2_prev[i][j][k]) && (nHe3_prev[i][j][k] <= nHe2_prev[i][j][k]) )  {
		nHe1[i][j][k] = (nHe1_prev[i][j][k] + clump[i][j][k]*recomb_He2[i][j][k]*ne[i][j][k]*(nHe[i][j][k] - nHe3[i][j][k])*dt)
					  /(1. + gamma_He1_tot[i][j][k]*dt + clump[i][j][k]*recomb_He2[i][j][k]*ne[i][j][k]*dt);
		nHe3[i][j][k] = (nHe3_prev[i][j][k]+ gamma_He2_tot[i][j][k]*(nHe[i][j][k] - nHe1[i][j][k])*dt)
					  /(1. + gamma_He2_tot[i][j][k]*dt + clump[i][j][k]*recomb_He3[i][j][k]*ne[i][j][k]*dt);
		nHe2[i][j][k] = nHe[i][j][k] - nHe1[i][j][k] - nHe3[i][j][k];
	}
	else {
		nHe2[i][j][k] = (nHe2_prev[i][j][k] + gamma_He1_tot[i][j][k]*(nHe[i][j][k] - nHe3[i][j][k])*dt
					  + clump[i][j][k]*recomb_He3[i][j][k]*ne[i][j][k]*(nHe[i][j][k] - nHe1[i][j][k])*dt)
					  /(1. + (gamma_He1_tot[i][j][k] + gamma_He2_tot[i][j][k])*dt + clump[i][j][k]*(recomb_He2[i][j][k]
					  + recomb_He3[i][j][k])*ne[i][j][k]*dt);
		nHe3[i][j][k] = (nHe3_prev[i][j][k] + gamma_He2_tot[i][j][k]*(nHe[i][j][k] - nHe1[i][j][k])*dt)
					  /(1.+ gamma_He2_tot[i][j][k]*dt + clump[i][j][k]*recomb_He3[i][j][k]*ne[i][j][k]*dt);
		nHe1[i][j][k] = nHe[i][j][k] - nHe2[i][j][k] - nHe3[i][j][k];
	}
	
	//update free electron and total number densities
	ne[i][j][k]    = nH2[i][j][k] + nHe2[i][j][k] + 2.*nHe3[i][j][k];
}

//Update the chemical abundances
void update_chem(void)  {
	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				
				//electron density derivative
				dne_dt[i][j][k]   = -ne[i][j][k];
				
				//recombination coefficients
				if (strcmp(rec_case, "A") == 0)  {
					recomb_H2[i][j][k]  = (float) alphaA_H2(temp[i][j][k]);
					recomb_He2[i][j][k] = (float) alphaA_He2(temp[i][j][k]);
					recomb_He3[i][j][k] = (float) alphaA_He3(temp[i][j][k]);
				}
				else  {
					recomb_H2[i][j][k]  = (float) alphaB_H2(temp[i][j][k]);
					recomb_He2[i][j][k] = (float) alphaB_He2(temp[i][j][k]);
					recomb_He3[i][j][k] = (float) alphaB_He3(temp[i][j][k]);
				}
				
				//diaelectric recombination of He II.  
				if ( (temp[i][j][k] >= 3e4) && (temp[i][j][k] <= 1e6) )  {
					recomb_He2[i][j][k] += (float) Dalpha_He2(temp[i][j][k]);
				} 
				
				//solve the backwards difference scheme iteratively to get converged electron density
				for (int m = 0; m < 5; m++)  {
					solve_ion(i, j, k);
				}
				
				//total number density & electron derivative
				n_tot[i][j][k] = nH[i][j][k] + nHe[i][j][k] + ne[i][j][k];
				dne_dt [i][j][k] += ne[i][j][k];
				dne_dt [i][j][k] /= dt;
				
				//update previous time step
				nH1_prev[i][j][k]  = nH1[i][j][k];
				nH2_prev[i][j][k]  = nH2[i][j][k];
				nHe1_prev[i][j][k] = nHe1[i][j][k];
				nHe2_prev[i][j][k] = nHe2[i][j][k];
				nHe3_prev[i][j][k] = nHe3[i][j][k];
				
				//update abundance fractions
				f_H1[i][j][k]  = nH1[i][j][k]/nH[i][j][k];
				f_H2[i][j][k]  = nH2[i][j][k]/nH[i][j][k];
				f_He1[i][j][k] = nHe1[i][j][k]/nHe[i][j][k];
				f_He2[i][j][k] = nHe2[i][j][k]/nHe[i][j][k];
				f_He3[i][j][k] = nHe3[i][j][k]/nHe[i][j][k];
			}
		}
	}
	}
}

float get_thermal_tstep(int i, int j, int k)  {
	float dT_dt =  2./3./k_B/n_tot[i][j][k]*(heat_rate[i][j][k] - cool_rate[i][j][k])
				- (2.*H(1./(1. + zz)) + dne_dt[i][j][k]/n_tot[i][j][k])*temp[i][j][k];
	float min_dt = 0.1*dt;
	return maxd(min_dt, mind(0.1*temp[i][j][k]/absd(dT_dt), dt));
}

//Update the temperature.  Either solve the temperature equation if thermal evolution is on, or set neutral
//cells to temp0 and ionized cells to 10^4 K.  
void update_thermal(void)  {
	float a = 1/(1 + zz);
	float t_temp[Nx][Ny][Nz];
	
	//if temperature evolution is on, solve the equation for T.  
	if (temp_ev == TRUE)  {
		#pragma omp parallel
		{
		#pragma omp for
		for (int i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  {
					//Commented 09/21/23. this prevents the gas from getting Treion heated during first 50% of a cells ion history
					//if (treion_flag[i][j][k] == FALSE)  { //set to the neutral temp before ionization
					//	temp[i][j][k] = temp_0;
					//}

					//Bug fixed  09/21/23
					if ( (treion_flag[i][j][k] == FALSE) && (subgrid == TRUE) )  { //set to the neutral temp before ionization
						temp[i][j][k] = temp_0;
					}

					//Added 09/25/21: turn off photo-heating and cooling terms (in regular RT mode) in cells with neutral fractions > 0.05 to avoid
					//incorrect heating/cooling in cells with low average temperatures.  
					//else if ((treion_flag[i][j][k] == TRUE) && (subgrid == FALSE) && (f_H1[i][j][k] > 0.05)){
					//	temp[i][j][k] = temp[i][j][k];
					//}
					else  { //evolve the temperature in the ionized part of the cell
						t_temp[i][j][k] = 0.;
						//Added back in 12/16/21: set post I-front temperature if sub-grid is off
						if (subgrid == FALSE)  {
							set_Treion_default(i, j, k);
						}
						while (t_temp[i][j][k] < dt)  {
							float dt_temp = mind(get_thermal_tstep(i, j, k), dt - t_temp[i][j][k]);
							update_heat_cool(i, j, k);
							temp[i][j][k] 	= (temp[i][j][k] + 2./3./k_B/n_tot[i][j][k]*(heat_rate[i][j][k] - cool_rate[i][j][k])*dt_temp)
											/(1. + 2.*H(1./(1. + zz))*dt_temp + 1./n_tot[i][j][k]*dne_dt[i][j][k]*dt_temp - 2.*ddelta_dt[i][j][k]/3./(rho[i][j][k]/avg_rho)*dt_temp);
							t_temp[i][j][k] += dt_temp;
						}
						//FIXED typo bug on 06/20/22
						temp[i][j][k] = maxd(temp[i][j][k], temp_min);
						
						//Added 08/31/22 - update previous values here if not in sub-grid mode (for vIF calc)
						if (subgrid == FALSE)  {
							nH1_prev[i][j][k]  = nH1[i][j][k];
							nH2_prev[i][j][k]  = nH2[i][j][k];
							nHe1_prev[i][j][k] = nHe1[i][j][k];
							nHe2_prev[i][j][k] = nHe2[i][j][k];
							nHe3_prev[i][j][k] = nHe3[i][j][k];
						}
						
					}
				}
			}
		}
		}
	}
	//if temp_ev is off and subgrid is on, get the temp from the subgrid file in ionized cells
	else if ( (temp_ev == FALSE) && (subgrid == TRUE) )  {
		#pragma omp parallel
		{
		#pragma omp for
		for (int i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  {
					if (fion[i][j][k] > 0.)  {
						interp_temp(i, j, k, avg_gamma[i][j][k]);
					}
					else  {
						temp[i][j][k] = temp_0;
					}
				}
			}
		}
		}
	}
	else  { //otherwise, neutral cells = temp0 and ionized cells = 10^4K.  
		#pragma omp parallel
		{
		#pragma omp for
		for (int i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  {
					if (nH1[i][j][k]/nH[i][j][k] > 0.5)  {
						temp[i][j][k] = temp_0;
					}
					else  {
						temp[i][j][k] = 1e4;
					} 
				}
			}
		}
		}
	}
}

//ADDED 05/26/22 to update the density field between hydro steps assuming rho ~ (1+z)^3 (ignores structure formation between hydro steps)
void update_cosmo_expansion(void)  {
	float cosmic_time_current  = t + cosmic_time(1./(1.+zinit));
	float cosmic_time_previous = t + cosmic_time(1./(1.+zinit)) - dt;
	
	float a_current  = scale_factor(cosmic_time_current);
	float a_prev     = scale_factor(cosmic_time_previous);
	
	printf("checking scale factors: %le\n", cosmic_time_current);
	printf("%le\n", cosmic_time_previous);
	printf("%le\n", a_current);
	printf("%le\n", a_prev);
	
	Lx *= a_current/a_prev;
	Ly *= a_current/a_prev;
	Lz *= a_current/a_prev;
	
	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < Nx; i++)  {
		dx[i] *= a_current/a_prev;
		dy[i] *= a_current/a_prev;
		dz[i] *= a_current/a_prev;
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				rho[i][j][k] *= pow(a_prev/a_current, 3.);
				nH[i][j][k]    = (1. - Y)*rho[i][j][k]/m_H;
				nHe[i][j][k]   = Y*rho[i][j][k]/(4.*m_H);
				if (fion[i][j][k] > 0.)  {
					nH2[i][j][k]   = f_H2[i][j][k]*nH[i][j][k]/fion[i][j][k];
					nHe2[i][j][k]  = f_He2[i][j][k]*nHe[i][j][k]/fion[i][j][k];
					nHe3[i][j][k]  = f_He3[i][j][k]*nHe[i][j][k]/fion[i][j][k];
				}
				else  {
					nH2[i][j][k] = nH[i][j][k]*(1. - fH1_0);
					nHe2[i][j][k] = nHe[i][j][k]*(1. - fHe1_0);
					nHe3[i][j][k] = 0.;
				}
				nH1[i][j][k] = nH[i][j][k] - nH2[i][j][k];
				avg_nh1[i][j][k] = nH1[i][j][k];
				nHe1[i][j][k] = nHe[i][j][k] - nHe2[i][j][k] - nHe3[i][j][k];
				avg_nhe1[i][j][k] = nHe1[i][j][k];
				ne[i][j][k]    = nH2[i][j][k] + nHe2[i][j][k] + 2.*nHe3[i][j][k];
				n_tot[i][j][k] = nH[i][j][k] + nHe[i][j][k] + ne[i][j][k];  
				
				//Added 08/31/22 - need to update previous values here too or C2Ray mode won't work.  
				nH1_prev[i][j][k]  = nH1[i][j][k];
				nH2_prev[i][j][k]  = nH2[i][j][k];
				nHe1_prev[i][j][k] = nHe1[i][j][k];
				nHe2_prev[i][j][k] = nHe2[i][j][k];
				nHe3_prev[i][j][k] = nHe3[i][j][k];
			}
		}
	}
	}
	avg_rho *= pow(a_prev/a_current, 3.);
	
	//Rescale the positions of the rays within their cells, if this is not the first hydro step.  
	#pragma omp parallel
	{
	#pragma omp for
	for (long int x = 0; x < size_buffer; x++)  {
		//check if ray exists and moved on the last time step
		if (ray_tags[x] == TRUE)  {
			ray[x].x  *= a_current/a_prev;
			ray[x].y  *= a_current/a_prev;
			ray[x].z  *= a_current/a_prev;
			
			//CHRIS 07/05/22: corner correction
			if (ray[x].x <= 0.)  {
					ray[x].x = min_dist*dx[0];
			}
			if (ray[x].x >= dx[0]*(1. - min_dist/2.)) {
					ray[x].x = dx[0]*(1. - min_dist);
			}
			if (ray[x].y <= 0.)  {
					ray[x].y = min_dist*dy[0];
			}
			if (ray[x].y >= dy[0]*(1. - min_dist/2.)) {
					ray[x].y = dy[0]*(1. - min_dist);
			}
			if (ray[x].z <= 0.)  {
					ray[x].z = min_dist*dz[0];
			}
			if (ray[x].z >= dz[0]*(1. - min_dist/2.)) {
					ray[x].z = dz[0]*(1. - min_dist);
			}
		}
	}
	}
}

//Move to the next hydro step by re-scaling quantities and updating/re-reading the input files.  Also output
//the results from the end of the previous hydro step. 
void update_hydro_step(void)  {
	string s1, gas, source, num_source, gas_out, ray_out, tre_out;
	
	float zz_prev = (float) zz;
	zz = hydro_steps[ihydro];
	string zz_string = to_string(zz);
	printf("zz_string: %s\n", zz_string.c_str());
	
	float zz_next = hydro_steps[ihydro+1];
	float aa = 1./(1. + zz);
	float aa_next = 1./(1. + zz_next);
	
	float a_rt = scale_factor(t + cosmic_time(1./(1.+zinit))); //ADDED 05/26/22
	
	printf("checking scale factor in update hydro: %le\n", a_rt);
	printf("%le\n", 1./a_rt - 1.);
	printf("%le\n", t);
	printf("%le\n", cosmic_time(1./(1.+zinit)));
	
	//If we are not only the last hydro step, update.  
	if (ihydro < num_hydro_steps - 1)  {
		
		printf("zz: %le\n", zz);
		printf("zz_next: %le\n", zz_next);
		
		printf("t: %le\n", cosmic_time(aa));
		printf("t_next: %le\n", cosmic_time(aa_next));
		
		//Append the next time step to the total simulation runtime and rescale the lengths. 
		//FIXED 07/26/22 to work correctly with new restart scheme.  
		if (step == 0)  {
			t_tot += cosmic_time(aa_next) - cosmic_time(a_rt);
		}
		else  {
			t_tot += cosmic_time(aa_next) - cosmic_time(aa);
		}

		
		Lx = Lx0*a_rt; //convert to proper units
		Ly = Ly0*a_rt;
		Lz = Lz0*a_rt;
	}
	
	//If this is the first time step, grab the initial source and density fields.  
	if (zz == zinit)  {
		get_ugas(ugas_file_init);
		
		#pragma omp parallel
		{
		#pragma omp for
		for (int i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  {
					rho[i][j][k] *= pow(a_rt, -3.);
				}
			}
		}
		}
		avg_rho *= pow(a_rt, -3.);
		
		read_source_catalog(source_field_init, num_src_file_init);
	}
		
	else  { //Otherwise, write the gas, ray, and treion output files. 
		if (ihydro > 0)  {
			if (zz < 10.)  {
				s1      = std::string(output_dir)+"gas_z=0";
				gas_out = s1 + zz_string;
				s1      = std::string(output_dir)+"ray_z=0";
				ray_out = s1 + zz_string;
				s1      = std::string(output_dir)+"tre_z=0";
				tre_out = s1 + zz_string;
			}
			else  {
				s1      = std::string(output_dir)+"gas_z=";
				gas_out = s1 + zz_string;
				s1      = std::string(output_dir)+"ray_z=";
				ray_out = s1 + zz_string;
				s1      = std::string(output_dir)+"tre_z=";
				tre_out = s1 + zz_string;
			}

			printf("gas_out: %s\n", gas_out.c_str());
			printf("ray_out: %s\n", ray_out.c_str());
			printf("tre_out: %s\n", tre_out.c_str());

			gas_out.resize(((int) gas_out.length()) - 2);
			ray_out.resize(((int) ray_out.length()) - 2);
			tre_out.resize(((int) tre_out.length()) - 2);

			const char *gas_output = gas_out.c_str();
			const char *ray_output = ray_out.c_str();
			const char *tre_output = tre_out.c_str();

			write_gas_binary(gas_output);
			//write_reduced_rays_binary(ray_output);
			//if ( (ihydro + snapinit > 54) && (ihydro + snapinit < 63) )  {
			if ( (ihydro + snapinit > 50) && (ihydro + snapinit < 59) )  {
				write_rays_binary(ray_output);
			}
			write_treion(tre_output);
		}
		
		if (ihydro < num_hydro_steps - 1)  {
			//Generate the names of the input files for the next step.
			if (ihydro + snapinit < 10)  {
					//gas    = "/nobackup/ccain5/pipeline_for_chris/real_75_mod/density.00"+to_string(ihydro + snapinit)+".75";
					gas        = std::string(hydro_base)+std::string(density_dir)+std::string(density_base)+"00"+to_string(ihydro + snapinit)+"."+to_string(Nx);
					source     = std::string(hydro_base)+std::string(source_dir)+std::string(source_base)+"00"+to_string(ihydro + snapinit);
					num_source = std::string(hydro_base)+std::string(source_dir)+std::string(sourcenum_base)+"00"+to_string(ihydro + snapinit);
			}
			else if ((ihydro + snapinit < 100) && (ihydro + snapinit >= 10)){
					gas        = std::string(hydro_base)+std::string(density_dir)+std::string(density_base)+"0"+to_string(ihydro + snapinit)+"."+to_string(Nx);
					source     = std::string(hydro_base)+std::string(source_dir)+std::string(source_base)+"0"+to_string(ihydro + snapinit);
					num_source = std::string(hydro_base)+std::string(source_dir)+std::string(sourcenum_base)+"0"+to_string(ihydro + snapinit);
			}
			else  {
					gas        = std::string(hydro_base)+std::string(density_dir)+std::string(density_base)+to_string(ihydro + snapinit)+"."+to_string(Nx);
					source     = std::string(hydro_base)+std::string(source_dir)+std::string(source_base)+to_string(ihydro + snapinit);
					num_source = std::string(hydro_base)+std::string(source_dir)+std::string(sourcenum_base)+to_string(ihydro + snapinit);
			}
			
			const char * ugas_file = gas.c_str();
			const char * source_field = source.c_str();
			const char * num_src_file = num_source.c_str();
			
			//read density field and sources.  
			get_ugas(ugas_file);
			printf("read ugas: %s\n", ugas_file);
                        printf("read source: %s\n", source_field);
                        printf("read num_source: %s\n", num_src_file);
	
			#pragma omp parallel
			{
			#pragma omp for
			for (int i = 0; i < Nx; i++)  {
				for (int j = 0; j < Ny; j++)  {
					for (int k = 0; k < Nz; k++)  {
						rho[i][j][k] *= pow(a_rt, -3.);
					}
				}
			}
			}
			avg_rho *= pow(a_rt, -3.);
			
			read_source_catalog(source_field, num_src_file);
			float avg_rho_prev = calc_vol_avg(rho_prev);
			
			//update hydro variables, negelcting advection
			//MODIFIED 04/05/21 to include update to drho_dt for temperature calculation
			#pragma omp parallel
			{
			#pragma omp for
			for (int i = 0; i < Nx; i++)  {
				for (int j = 0; j < Ny; j++)  {
					for (int k = 0; k < Nz; k++)  {
						dx[i]          = Lx/Nx*kpc_to_cm;
						dy[j]          = Ly/Ny*kpc_to_cm;
						dz[k]          = Lz/Nz*kpc_to_cm;
						nH[i][j][k]    = (1. - Y)*rho[i][j][k]/m_H;
						nHe[i][j][k]   = Y*rho[i][j][k]/(4.*m_H);
						if (fion[i][j][k] > 0.)  {
							nH2[i][j][k]   = f_H2[i][j][k]*nH[i][j][k]/fion[i][j][k];
							nHe2[i][j][k]  = f_He2[i][j][k]*nHe[i][j][k]/fion[i][j][k];
							nHe3[i][j][k]  = f_He3[i][j][k]*nHe[i][j][k]/fion[i][j][k];
						}
						else  {
							nH2[i][j][k] = nH[i][j][k]*(1. - fH1_0);
							nHe2[i][j][k] = nHe[i][j][k]*(1. - fHe1_0);
							nHe3[i][j][k] = 0.;
						}
						nH1[i][j][k] = nH[i][j][k] - nH2[i][j][k];
						avg_nh1[i][j][k] = nH1[i][j][k];
						nHe1[i][j][k] = nHe[i][j][k] - nHe2[i][j][k] - nHe3[i][j][k];
						avg_nhe1[i][j][k] = nHe1[i][j][k];
						ne[i][j][k]    = nH2[i][j][k] + nHe2[i][j][k] + 2.*nHe3[i][j][k];
						n_tot[i][j][k] = nH[i][j][k] + nHe[i][j][k] + ne[i][j][k];    
						
						 if (avg_rho_prev > 0.)  {
								ddelta_dt[i][j][k]  = (rho[i][j][k]/avg_rho - rho_prev[i][j][k]/avg_rho_prev)/(cosmic_time(aa_next) - cosmic_time(aa));
						 }
						
						//Update previous values of ionization states
						rho_prev[i][j][k]  = rho[i][j][k];
						nH1_prev[i][j][k]  = nH1[i][j][k];
						nH2_prev[i][j][k]  = nH2[i][j][k];
						nHe1_prev[i][j][k] = nHe1[i][j][k];
						nHe2_prev[i][j][k] = nHe2[i][j][k];
						nHe3_prev[i][j][k] = nHe3[i][j][k];
					}
				}
			}
			}
		
		}
	}
	printf("Beginning hydro step: %d\n", ihydro);
	ihydro += 1;
}
