#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "io_funcs.cc"

//Initialize the gas grid
void init(void)  {
	
	tinit     = cosmic_time(1./(1 + zinit));
	
	//if restart, read the grid and treion restart files and the gas file
	if (input_grid == TRUE)  {
		read_grid_binary(start_file_init);
		read_treion(treion_file_init);
	}
	else  { //otherwise, read the gas data from the input file
		
		//spatial grid   
		for (int i = 0; i < Nx; i++)  {
			dx[i] = Lx/Nx*kpc_to_cm;
		}
		for (int j = 0; j < Ny; j++)  {
			dy[j] = Ly/Ny*kpc_to_cm;
		}
		for (int k = 0; k < Nz; k++)  {
			dz[k] = Lz/Nz*kpc_to_cm;
		}
		#pragma omp parallel
		{
		#pragma omp for 
		for (int i = 0; i < Nx; i++)  { //ionization fractions and initial temp
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  {
					f_H1[i][j][k]  = fH1_0;
					f_H2[i][j][k]  = 1. - f_H1[i][j][k];
					f_He1[i][j][k] = fHe1_0;
					f_He2[i][j][k] = fHe2_0;
					f_He3[i][j][k] = 1. - f_He1[i][j][k] - f_He2[i][j][k];
					
					temp[i][j][k] = temp_0;
					
					//ADDED 05/27/22 - initialize Fahad's MP-Gadget heating variables
					if (MP_GADGET_HEAT == TRUE)  {
						gadget_ini_heat[i][j][k] = 0.0;
						gadget_heat[i][j][k] = 0.0;
						gadget_gamma[i][j][k] = 0.0;
					}
					
				}
			}
		}
		}
	}
	
	//MODIFIED 04/08/21 to update number densities properly
	//compute other important quantities
	#pragma omp parallel
	{
	#pragma omp for 
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				//number densities
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
				
				//IF stuff - added 08/31/22
				IF_speed[i][j][k] = 0.;
				Inc_Flux[i][j][k] = 0.;
				
				//previous time step
				nH1_prev[i][j][k]  = nH1[i][j][k];
				nH2_prev[i][j][k]  = nH2[i][j][k];
				nHe1_prev[i][j][k] = nHe1[i][j][k];
				nHe2_prev[i][j][k] = nHe2[i][j][k];
				nHe3_prev[i][j][k] = nHe3[i][j][k];
				
				//on the fly stepping
				f_H1_step[i][j][k]  = f_H1[i][j][k];
				f_H2_step[i][j][k]  = f_H2[i][j][k];
				f_He1_step[i][j][k] = f_He1[i][j][k];
				f_He2_step[i][j][k] = f_He2[i][j][k];
				f_He3_step[i][j][k] = f_He3[i][j][k];
				
				//recombination coefficients
				if (strcmp(rec_case, "A") == 0)  {
					recomb_H2[i][j][k]  = alphaA_H2(temp[i][j][k]);
					recomb_He2[i][j][k] = alphaA_He2(temp[i][j][k]);
					recomb_He3[i][j][k] = alphaA_He3(temp[i][j][k]);
				}
				else  {
					recomb_H2[i][j][k]  = alphaB_H2(temp[i][j][k]);
					recomb_He2[i][j][k] = alphaB_He2(temp[i][j][k]);
					recomb_He3[i][j][k] = alphaB_He3(temp[i][j][k]);
				}
				
				if (input_grid == TRUE)  {
					avg_gamma[i][j][k]     = gamma_H1_tot[i][j][k];
					avg_gamma_prev[i][j][k] = avg_gamma[i][j][k];
					avg_gamma_he1[i][j][k] = gamma_He1_tot[i][j][k];
				}
				else  {
					treion[i][j][k] = 0.;
					treion_out[i][j][k] = 0.;
					//initialize MFP
					mfp[i][j][k]        = 1./sigmapi_H1(nu_phot[0])/avg_nh1[i][j][k];
					mfp_prev[i][j][k]   = 1./sigmapi_H1(nu_phot[0])/avg_nh1[i][j][k];
					//initialize MFP
					mfp_912[i][j][k]        = 1./sigmapi_H1(13.6/h_eV)/avg_nh1[i][j][k];
					mfp_912_prev[i][j][k]   = 1./sigmapi_H1(13.6/h_eV)/avg_nh1[i][j][k];
				}
				//compute the time derivative of the electron density
				float dnH2_dt  = gamma_H1_tot[i][j][k]*nH1[i][j][k] - recomb_H2[i][j][k]*ne[i][j][k]*nH2[i][j][k];
				float dnHe2_dt = nHe1[i][j][k]*gamma_He1_tot[i][j][k] - recomb_He2[i][j][k]*ne[i][j][k]*nHe2[i][j][k]
								  + recomb_He3[i][j][k]*ne[i][j][k]*nHe3[i][j][k] - nHe2[i][j][k]*gamma_He2_tot[i][j][k];
				float dnHe3_dt = nHe2[i][j][k]*gamma_He2_tot[i][j][k] - recomb_He3[i][j][k]*ne[i][j][k]*nHe3[i][j][k];
				dne_dt[i][j][k]   = dnH2_dt + dnHe2_dt + 2*dnHe3_dt;
				ddelta_dt[i][j][k]  = 0.;
			}
		}
	}
	}
}

