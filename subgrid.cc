#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <string>
#include "init_funcs.cc"
using namespace std;

//ADDED 03/14 to calculate the correct factor by which to boost the MFP over the value at 912 ang in the RH sims.  Assumes a column density slope beta_col
//and average over 1/mfp_nu to get the appropriate boost factor
void calc_mfp_boost(void)  {
	for (int i = 0; i < 5; i++)  {
		mfp_boost += (1./5.)*pow(nu_rh[i]/(13.6/h_eV), -3.*(beta_col - 1.));
	}
	mfp_boost = 1./mfp_boost;
}

//MODIFIED 03/14 to correct the boost to the 912 angstrom MFP
//MODIFIED 03/15 to include the 912 angstrom MFP and the frequency averaged one from the new subgrid files
//Function for importing the data from a RadHydro sub-grid model file. I have code elsewhere that makes these.  
void import_subgrid_file(FILE *file, int i, int j, int k)  { //i = zre, j = delta, k = gamma
	double temp_double;
	int gam_offset = 1;
	
	printf("Importing file: %d %d %d\n", i, j, k);
	
	if (file == NULL)  {
		printf("Subgrid file not found\n");
	}
	for (int m = 0; m < Nt; m++)  {
		fscanf(file, "%le", &temp_double);
		time_RH [m][i][j][k+gam_offset] = (float) temp_double; //time since ionization
		fscanf(file, "%le", &temp_double);
		mfp_RH  [0][m][i][j][k+gam_offset] = (float) temp_double; //Mean free path frequency averaged
		fscanf(file, "%le", &temp_double);
		mfp_RH  [1][m][i][j][k+gam_offset] = (float) temp_double; //Mean free path at 912 angstroms
		fscanf(file, "%le", &temp_double);
		fh1_RH  [m][i][j][k+gam_offset] = (float) temp_double; //Total HI fraction
		fscanf(file, "%le\n", &temp_double);
		temp_RH [m][i][j][k+gam_offset] = (float) temp_double; //Temperature
		time_RH[m][i][j][k+gam_offset] *= yr_to_s; //convert to seconds
		mfp_RH[0][m][i][j][k+gam_offset]  *= (float) (1e3*kpc_to_cm/hh); //convert to comoving cm
		mfp_RH[1][m][i][j][k+gam_offset]  *= (float) (1e3*kpc_to_cm/hh); //convert to comoving cm
		mfp_RH[0][m][i][j][k+gam_offset] /= const_clump;
		mfp_RH[1][m][i][j][k+gam_offset] /= const_clump; //scale down the MFP by 1/clumping factor
	}
	fclose(file);
}

void import_delta_over_sigma_table(void)  {
	double temp_double;
	FILE *file = NULL;
	//const char delta_over_sigma_table[] = "subgrid_files/delta_over_sigma_table.txt";
	string delta_over_sigma_table = std::string(subgrid_dir) + "delta_over_sigma_table.txt";
	const char * delta_over_sigma_table_file = delta_over_sigma_table.c_str();
	
	printf("opening delta_over_sigma interpolation table\n");
	file = fopen(delta_over_sigma_table_file, "r");
	
	if (file == NULL)  {
		printf("Delta over sigma table not found\n");
	}
	
	//read in redshift and delta/sigma
	for (int i = 0; i < n_dlin-1; i++)  {
		fscanf(file, "%le", &temp_double);
		t_deltaNL[i] = (float) temp_double;
		if (i == 50)  {
			printf("t_deltaNL[%d]: %le\n", i, t_deltaNL[i]);
		}
	}
	fscanf(file, "%le\n", &temp_double);
	t_deltaNL[n_dlin-1] = (float) temp_double;
	
	for (int i = 0; i < n_dlin-1; i++)  {
		fscanf(file, "%le", &temp_double);
		delta_over_sigma_interp[i] = (float) temp_double;
		if (i == 50)  {
			printf("delta_over_sigma_interp[%d]: %le\n", i, delta_over_sigma_interp[i]);
		}
	}
	fscanf(file, "%le\n", &temp_double);
	delta_over_sigma_interp[n_dlin-1] = (float) temp_double;
	
	for (int i = 0; i < n_dlin; i++)  {
		for (int j = 0; j < n_dlin-1; j++)  {
			fscanf(file, "%le", &temp_double);
			deltaNL_interp[i][j] = (float) temp_double;
			if ( (i == 50) && (j < 10) )  {
				printf("deltaNL_interp[%d][%d]: %le\n", i, j, deltaNL_interp[i][j]);
			}
		}
		fscanf(file, "%le\n", &temp_double);
		deltaNL_interp[i][n_dlin-1] = (float) temp_double;
	}
}

//MODIFIED on 03/09 to include a cap deltaNL on the over-dense case, to avoid turnaround.  
void delta_dom_to_deltaNL_dom(float time)  {
	deltaNL_dom[1] = 0.;
	float zcap = 5.5;
	float tcap = cosmic_time(1./(1. + zcap));
	
	int ti_m, ti_p;
	int tcap_m, tcap_p;
	
	for (int i = 0; i < n_dlin-1; i++)  {
		if ( (t_deltaNL[i] <= time) && (t_deltaNL[i+1] >= time) )  {
			ti_m = i;
			ti_p = i+1;
		}
		if ( (t_deltaNL[i] <= tcap) && (t_deltaNL[i+1] >= tcap) )  {
			tcap_m = i;
			tcap_p = i+1;
		}
	}
	
	float deltaNL_dom0_m = interpolate(delta_over_sigma_interp, deltaNL_interp[ti_m], delta_dom[0], n_dlin);
	float deltaNL_dom0_p = interpolate(delta_over_sigma_interp, deltaNL_interp[ti_p], delta_dom[0], n_dlin);
	
	float deltaNL_dom2_m = interpolate(delta_over_sigma_interp, deltaNL_interp[ti_m], delta_dom[2], n_dlin);
	float deltaNL_dom2_p = interpolate(delta_over_sigma_interp, deltaNL_interp[ti_p], delta_dom[2], n_dlin);
	
	float deltaNL_cap_m  = interpolate(delta_over_sigma_interp, deltaNL_interp[tcap_m], delta_dom[2], n_dlin);
	float deltaNL_cap_p  = interpolate(delta_over_sigma_interp, deltaNL_interp[tcap_p], delta_dom[2], n_dlin);
	
	float deltaNL_cap = deltaNL_cap_m  + (deltaNL_cap_p - deltaNL_cap_m)  /(t_deltaNL[tcap_p] - t_deltaNL[tcap_m])*(tcap - t_deltaNL[tcap_m]);
	deltaNL_dom[0]    = deltaNL_dom0_m + (deltaNL_dom0_p - deltaNL_dom0_m)/(t_deltaNL[ti_p] - t_deltaNL[ti_m])*(time - t_deltaNL[ti_m]);
	deltaNL_dom[2]    = deltaNL_dom2_m + (deltaNL_dom2_p - deltaNL_dom2_m)/(t_deltaNL[ti_p] - t_deltaNL[ti_m])*(time - t_deltaNL[ti_m]);
	deltaNL_dom[2]    = mind(deltaNL_dom[2], deltaNL_cap);
}

void init_radhydro_data(void)  {
	FILE *file = NULL;
	string model_dir, ext;
	
	//parameter domain available from the RadHydro runs
	zre_dom[0]   = 12.;
	zre_dom[1]   = 8.;
	zre_dom[2]   = 6.;
	
	delta_dom[0] = -pow(3., 0.5);
	delta_dom[1] = 0.;
	delta_dom[2] = pow(3., 0.5);
	
	gamma_dom[0] = 3e-15;
	gamma_dom[1] = 3e-14;
	gamma_dom[2] = 3e-13;
	gamma_dom[3] = 3e-12;
	gamma_dom[4] = 3e-11;
	
	string gamma_tags[Ngamma]   = {"G14", "G13", "G12"};
	string zre_tags[Nzre]       = {"zre12", "zre8", "zre6"};
	string density_tags[Ndelta] = {"dm3", "d0", "dp3"};
	
	//Specify the subfolder to use for the sub-grid model
	if (relaxed == TRUE)  {
		model_dir = "relaxed/";
		ext 	  = "_relaxed";
	}
	else if (equilibrium == TRUE)  {
		model_dir = "equilibrium/";
		ext		  = "_noclump";
	}
	else  {
		model_dir = "full_model/";
		ext  	  = "";
	}
	
	for (int q = 0; q < Ngamma; q++)  {
		for (int m = 0; m < Nzre; m++)  {
			for (int p = 0; p < Ndelta; p++)  {
				string subgrid_file_name = gamma_tags[q] + "_" + zre_tags[m] + "_" + density_tags[p] + ext + ".txt";
				string subgrid_file_address = std::string(subgrid_dir) + model_dir + subgrid_file_name;
				
				const char * subgrid_file = subgrid_file_address.c_str();
				file = fopen(subgrid_file, "r");
				import_subgrid_file(file, m, p, q);
			}
		}
	}
	
	//convert zreion to treion
	for (int i = 0; i < Nzre; i++)  {
		tre_dom[i] = cosmic_time(1. / (1. + zre_dom[i]));
	} 
}

//Get time interpolation index
int get_tradhydro_index(float tsince)  {
	int trh_index = (int) floor(tsince/time_RH[Nt-1][0][0][0]*(Nt-1));
	return max(0, min(trh_index, Nt-2));
}

//convert non-linear overdensity to linearly extrapolated overdensity / sigma
float deltaNL_to_deltasig(float deltaNL, float redshift)  {
	float dcmode = 5.10655;
	float sigma0 = dcmode/pow(3.,0.5);
	float d_of_z = 1./(1 + redshift)/0.78205; //assume matter domainated scaling + lambda correction
	
	//Mo & White 1996
	float delta = -1.35*pow(1. + deltaNL, -2./3.) + 0.78785*pow(1. + deltaNL, -0.58661) - 1.12431*pow(1. + deltaNL, -0.5) + 1.68647;
	
	return delta/d_of_z/sigma0; 
}

//Interpolate in log-log over the treion and delta/sigma variables to get the appropriate values of a variable
//at the three gamma gridpoints.  
//NOTE 03/13: using cosmic time as the interpolation variable here in lieu of treion means that the tcross correction does not explicitly affect
//the time variable here.  However, since the tsince of the input MFP is lower and the cosmic time remains fixed, treion is implicitly changed to the 
//correct value.  Even though tre + ts doesn't change the other terms are all modified appropriately.  
void get_interpolate_gamma_f(int i, int j, int k, float tsince, float tre, float deltaNL, int ind_tre[2], int ind_delta[2], int ind_gamma[Ngamma], float f_tsince[2][2][Ngamma], float f_gamma_interp[Ngamma+2], bool fractype)  {
		
	float f_tre_del0_gam0 = log10(f_tsince[0][0][0]) + (log10(f_tsince[1][0][0]) - log10(f_tsince[0][0][0]))
							   *(log10(tre + tsince) - log10(tre_dom[ind_tre[0]] + tsince))
								/(log10(tre_dom[ind_tre[1]] + tsince) - log10(tre_dom[ind_tre[0]] + tsince));
	float f_tre_del1_gam0 = log10(f_tsince[0][1][0]) + (log10(f_tsince[1][1][0]) - log10(f_tsince[0][1][0]))
								*(log10(tre + tsince) - log10(tre_dom[ind_tre[0]] + tsince))
								/(log10(tre_dom[ind_tre[1]] + tsince) - log10(tre_dom[ind_tre[0]] + tsince));
	float f_tre_del0_gam1 = log10(f_tsince[0][0][1]) + (log10(f_tsince[1][0][1]) - log10(f_tsince[0][0][1]))
							   *(log10(tre + tsince) - log10(tre_dom[ind_tre[0]] + tsince))
								/(log10(tre_dom[ind_tre[1]] + tsince) - log10(tre_dom[ind_tre[0]] + tsince));
	float f_tre_del1_gam1 = log10(f_tsince[0][1][1]) + (log10(f_tsince[1][1][1]) - log10(f_tsince[0][1][1]))
							   *(log10(tre + tsince) - log10(tre_dom[ind_tre[0]] + tsince))
								/(log10(tre_dom[ind_tre[1]] + tsince) - log10(tre_dom[ind_tre[0]] + tsince));
	float f_tre_del0_gam2 = log10(f_tsince[0][0][2]) + (log10(f_tsince[1][0][2]) - log10(f_tsince[0][0][2]))
							   *(log10(tre + tsince) - log10(tre_dom[ind_tre[0]] + tsince))
								/(log10(tre_dom[ind_tre[1]] + tsince) - log10(tre_dom[ind_tre[0]] + tsince));
	float f_tre_del1_gam2 = log10(f_tsince[0][1][2]) + (log10(f_tsince[1][1][2]) - log10(f_tsince[0][1][2]))
							   *(log10(tre + tsince) - log10(tre_dom[ind_tre[0]] + tsince))
							   /(log10(tre_dom[ind_tre[1]] + tsince) - log10(tre_dom[ind_tre[0]] + tsince));    
	
	float f_tre_del_gam0  = f_tre_del0_gam0 + (f_tre_del1_gam0 - f_tre_del0_gam0)
								*(log10(1. + deltaNL) - log10(1. + deltaNL_dom[ind_delta[0]]))
								/(log10(1. + deltaNL_dom[ind_delta[1]]) - log10(1. + deltaNL_dom[ind_delta[0]]));
	float f_tre_del_gam1  = f_tre_del0_gam1 + (f_tre_del1_gam1 - f_tre_del0_gam1)
								*(log10(1. + deltaNL) - log10(1. + deltaNL_dom[ind_delta[0]]))
								/(log10(1. + deltaNL_dom[ind_delta[1]]) - log10(1. + deltaNL_dom[ind_delta[0]]));
	float f_tre_del_gam2  = f_tre_del0_gam2 + (f_tre_del1_gam2 - f_tre_del0_gam2)
								*(log10(1. + deltaNL) - log10(1. + deltaNL_dom[ind_delta[0]]))
								/(log10(1. + deltaNL_dom[ind_delta[1]]) - log10(1. + deltaNL_dom[ind_delta[0]]));
	
	f_gamma_interp[1] = pow(10., f_tre_del_gam0);
	f_gamma_interp[2] = pow(10., f_tre_del_gam1);
	f_gamma_interp[3] = pow(10., f_tre_del_gam2);
	
	if (fractype == TRUE)  {
		f_gamma_interp[1] = mind(f_gamma_interp[1], max_subgrid_frac);
		f_gamma_interp[2] = mind(f_gamma_interp[2], max_subgrid_frac);
		f_gamma_interp[3] = mind(f_gamma_interp[3], max_subgrid_frac);
	}
}

//given a treion, delta, and t - treion, choose the relevant part of the grid and interpolate over those three variables
//to get f(gamma).  This operation only needs to be done once per time step for each cell that is ionized.  The resultant 
//f(gamma) will be interpolated during each iteration.  Also interpolate the slopes here for the i-front correction if necessary.
//MODIFIED 03/15 to interpolate over both MFP's.   
void get_gamma_interp(int i, int j, int k, float tre, float deltaNL, float tsince)  {
	int ind_tre[2], ind_delta[2], ind_gamma[Ngamma];
	int ind_tre_tcross[2];
	
	//Added 08/19/22
	int ind_tre_tcross_new[Ntcross][2];
	
	//interpolation grids
	float mfp_tsince_b   [2][2][2][Ngamma];
	float mfp_tsince_f   [2][2][2][Ngamma];
	float fh1_tsince     [2][2][Ngamma];
	float temp_tsince    [2][2][Ngamma];
	
	//correction for crossing time
	float mfp_tsince_tcross_b   [2][2][2][Ngamma];
	float mfp_tsince_tcross_f   [2][2][2][Ngamma];
	float fh1_tsince_tcross     [2][2][Ngamma];
	float temp_tsince_tcross    [2][2][Ngamma];
	
	//Added 08/19/22
	float mfp_tsince_tcross_new_b   [Ntcross][2][2][2][Ngamma];
	float mfp_tsince_tcross_new_f   [Ntcross][2][2][2][Ngamma];
	float fh1_tsince_tcross_new     [Ntcross][2][2][Ngamma];
	float temp_tsince_tcross_new    [Ntcross][2][2][Ngamma];
	
	//Added 08/19/22: use tcross bins
	float mfp_t_tcross_new_gamma_interp_b    [Ntcross][2][Ngamma+2];
	float mfp_t_tcross_new_gamma_interp_f    [Ntcross][2][Ngamma+2];
	float fh1_t_tcross_new_gamma_interp      [Ntcross][Ngamma+2];
	float temp_t_tcross_new_gamma_interp     [Ntcross][Ngamma+2];
	
	//require zreion to be 12 or less, so that no extrapolation happens
	tre = maxd(cosmic_time(1./(1.+12)), tre);
	
	//select the relevant parts of the domain
	if (tre > cosmic_time(1. / (1. + 8.))) {
		ind_tre[1] = 2;
		ind_tre[0] = 1;
	}
	else  {
		ind_tre[1] = 1;
		ind_tre[0] = 0;
	}
	
	if (tre + tcross[i][j][k] > cosmic_time(1. / (1. + 8.))) {
		ind_tre_tcross[1] = 2;
		ind_tre_tcross[0] = 1;
	}
	else  {
		ind_tre_tcross[1] = 1;
		ind_tre_tcross[0] = 0;
	}
	
	//Added 08/19/22
	for (int tt = 0; tt < Ntcross; tt++)  {
		if (tre + tcross_new[i][j][k][tt] > cosmic_time(1. / (1. + 8.))) {
			ind_tre_tcross_new[tt][1] = 2;
			ind_tre_tcross_new[tt][0] = 1;
		}
		else  {
			ind_tre_tcross_new[tt][1] = 1;
			ind_tre_tcross_new[tt][0] = 0;
		}
	}
	
	if (deltaNL_to_deltasig(deltaNL, zz) < 0.)  {
		ind_delta[1] = 1;
		ind_delta[0] = 0;
	}
	else  {
		ind_delta[1] = 2;
		ind_delta[0] = 1;
	}
    
	ind_gamma[2] = 3;
	ind_gamma[1] = 2;
	ind_gamma[0] = 1;
	
	//Modified 08/19/22
	int itb = get_tradhydro_index(tsince - dt + 1e5);
	int itf = get_tradhydro_index(tsince);
	
	//Added 08/19/22
	int itbtc_new[Ntcross];
	int itftc_new[Ntcross];
	
	for (int tt = 0; tt < Ntcross; tt++)  {
		itbtc_new[tt] = 0;
		itftc_new[tt] = 0;
		
		float tc_new_tt = 0.; 
		if (tt == 0)  {
			tc_new_tt = tcross_new[i][j][k][0]/2.;
		}
		else  {
			tc_new_tt = tcross_new[i][j][k][tt - 1] + (tcross_new[i][j][k][tt] - tcross_new[i][j][k][tt - 1])/2.;
		}
		
		if (tcross_flag[i][j][k] == TRUE)  {
			itbtc_new[tt] = get_tradhydro_index(tsince - tc_new_tt - dt + 1e5);
			itftc_new[tt] = get_tradhydro_index(tsince - tc_new_tt);
		}
	}
	
	
	//interpolate over the time domain to get each of the parameters at the different tre, delta/sigma and gamma values.  
	//Also interpolate the slope corrections if the ifront correction is turned on.  Use linear interpolation in time.  
	for (int m = 0; m < 2; m++)  {
		for (int p = 0; p < 2; p++)  {
			for (int q = 0; q < Ngamma; q++)  {
				
				//Modified 08/19/22
				int im = ind_tre[m];
				int ip = ind_delta[p];
				int iq = ind_gamma[q];
				
				for (int n = 0; n < 2; n++)  { //perform the same operation on both MFP's
					mfp_tsince_b[n][m][p][q] = mfp_RH[n][itb][im][ip][iq] + (mfp_RH[n][itb+1][im][ip][iq] - mfp_RH[n][itb][im][ip][iq])
										/ (time_RH[itb+1][im][ip][iq] - time_RH[itb][im][ip][iq]) * (tsince - dt - time_RH[itb][im][ip][iq]);
					mfp_tsince_f[n][m][p][q] = mfp_RH[n][itf][im][ip][iq] + (mfp_RH[n][itf+1][im][ip][iq] - mfp_RH[n][itf][im][ip][iq])
										/ (time_RH[itf+1][im][ip][iq] - time_RH[itf][im][ip][iq]) * (tsince - time_RH[itf][im][ip][iq]);
				}
				fh1_tsince[m][p][q]  = fh1_RH[itf][im][ip][iq] + (fh1_RH[itf+1][im][ip][iq] - fh1_RH[itf][im][ip][iq])
									 / (time_RH[itf+1][im][ip][iq] - time_RH[itf][im][ip][iq]) * (tsince - time_RH[itf][im][ip][iq]);
				temp_tsince[m][p][q] = temp_RH[itf][im][ip][iq] + (temp_RH[itf+1][im][ip][iq] - temp_RH[itf][im][ip][iq])
									 / (time_RH[itf+1][im][ip][iq] - time_RH[itf][im][ip][iq]) * (tsince - time_RH[itf][im][ip][iq]);
				
				//min/max critertia
				for (int n = 0; n < 2; n++)  {
					mfp_tsince_b[n][m][p][q]	= maxd(mfp_min, mfp_tsince_b[n][m][p][q]);
					mfp_tsince_f[n][m][p][q]	= maxd(mfp_min, mfp_tsince_f[n][m][p][q]);
				}
				fh1_tsince[m][p][q] 	= mind(maxd(fh1_tsince[m][p][q], 0.), 1.);
				temp_tsince[m][p][q] 	= maxd(temp_min, temp_tsince[m][p][q]);
				
				//if the i-front has started to cross a cell, interpolate the i-front correction slopes
				//worry about changing this later
				if (tcross_flag[i][j][k] == TRUE)  {
				
				//Added 08/19/22: get interp values over all tcross bins.  
					for (int tt = 0; tt < Ntcross; tt++)  {
						
						//Modified 08/19/22
						int imtc_new = ind_tre_tcross_new[tt][m];
						
						float tc_new_tt = 0.; 
						if (tt == 0)  {
							tc_new_tt = tcross_new[i][j][k][0]/2.;
						}
						else  {
							tc_new_tt = tcross_new[i][j][k][tt - 1] + (tcross_new[i][j][k][tt] - tcross_new[i][j][k][tt - 1])/2.;
						}
						
						for (int n = 0; n < 2; n++)  { 
							float mfp_tsince_tc_new_b  = mfp_RH[n][itbtc_new[tt]][imtc_new][ip][iq] + (mfp_RH[n][itbtc_new[tt]+1][imtc_new][ip][iq] - mfp_RH[n][itbtc_new[tt]][imtc_new][ip][iq])
											/ (time_RH[itbtc_new[tt]+1][imtc_new][ip][iq] - time_RH[itbtc_new[tt]][imtc_new][ip][iq]) * ((tsince - dt - tc_new_tt) - time_RH[itbtc_new[tt]][imtc_new][ip][iq]);
							float mfp_tsince_tc_new_f  = mfp_RH[n][itftc_new[tt]][imtc_new][ip][iq] + (mfp_RH[n][itftc_new[tt]+1][imtc_new][ip][iq] - mfp_RH[n][itftc_new[tt]][imtc_new][ip][iq])
											/ (time_RH[itftc_new[tt]+1][imtc_new][ip][iq] - time_RH[itftc_new[tt]][imtc_new][ip][iq]) * ((tsince - tc_new_tt) - time_RH[itftc_new[tt]][imtc_new][ip][iq]);
							mfp_tsince_tcross_new_b[tt][n][m][p][q]   = maxd(mfp_min, mfp_tsince_tc_new_b);
							mfp_tsince_tcross_new_f[tt][n][m][p][q]   = maxd(mfp_min, mfp_tsince_tc_new_f);
						}
						float fh1_tsince_tc_new  = fh1_RH[itftc_new[tt]][imtc_new][ip][iq] + (fh1_RH[itftc_new[tt]+1][imtc_new][ip][iq] - fh1_RH[itftc_new[tt]][imtc_new][ip][iq])
										 / (time_RH[itftc_new[tt]+1][imtc_new][ip][iq] - time_RH[itftc_new[tt]][imtc_new][ip][iq]) * ((tsince - tc_new_tt) - time_RH[itftc_new[tt]][imtc_new][ip][iq]);
						float temp_tsince_tc_new  = temp_RH[itftc_new[tt]][imtc_new][ip][iq] + (temp_RH[itftc_new[tt]+1][imtc_new][ip][iq] - temp_RH[itftc_new[tt]][imtc_new][ip][iq])
										 / (time_RH[itftc_new[tt]+1][imtc_new][ip][iq] - time_RH[itftc_new[tt]][imtc_new][ip][iq]) * ((tsince - tc_new_tt) - time_RH[itftc_new[tt]][imtc_new][ip][iq]);
						
						fh1_tsince_tcross_new[tt][m][p][q]     = mind(maxd(fh1_tsince_tc_new, 0.), 1.);
						temp_tsince_tcross_new[tt][m][p][q]    = maxd(temp_min, temp_tsince_tc_new);
					}
				}
			}
		}
	}
	
	//Interpolate over the treion and delta/sigma variables to get each variable as a function of gamma only.  
	for (int n = 0; n < 2; n++)  {
		get_interpolate_gamma_f(i, j, k, tsince, tre, deltaNL, ind_tre, ind_delta, ind_gamma, mfp_tsince_b[n], mfp_gamma_interp_b[n][i][j][k], FALSE);
		get_interpolate_gamma_f(i, j, k, tsince, tre, deltaNL, ind_tre, ind_delta, ind_gamma, mfp_tsince_f[n], mfp_gamma_interp_f[n][i][j][k], FALSE);
	}
	get_interpolate_gamma_f(i, j, k, tsince, tre, deltaNL, ind_tre, ind_delta, ind_gamma, fh1_tsince, fh1_gamma_interp[i][j][k], TRUE);
	get_interpolate_gamma_f(i, j, k, tsince, tre, deltaNL, ind_tre, ind_delta, ind_gamma, temp_tsince, temp_gamma_interp[i][j][k], FALSE);
	
	//enforce desired scaling above and below limits
	for (int n = 0; n < 2; n++)  {
		mfp_gamma_interp_b[n][i][j][k][0] = pow(gamma_dom[0]/gamma_dom[1], 1.)*mfp_gamma_interp_b[n][i][j][k][1];
		mfp_gamma_interp_f[n][i][j][k][0] = pow(gamma_dom[0]/gamma_dom[1], 1.)*mfp_gamma_interp_f[n][i][j][k][1];
		mfp_gamma_interp_b[n][i][j][k][4] = pow(gamma_dom[4]/gamma_dom[3], 1.)*mfp_gamma_interp_b[n][i][j][k][3];
		mfp_gamma_interp_f[n][i][j][k][4] = pow(gamma_dom[4]/gamma_dom[3], 1.)*mfp_gamma_interp_f[n][i][j][k][3];
	}
	
	fh1_gamma_interp[i][j][k][0] = pow(gamma_dom[0]/gamma_dom[1],  0.)*fh1_gamma_interp[i][j][k][1];
	fh1_gamma_interp[i][j][k][4] = pow(gamma_dom[4]/gamma_dom[3], -1.)*fh1_gamma_interp[i][j][k][3];
	
	temp_gamma_interp[i][j][k][0] = pow(gamma_dom[0]/gamma_dom[1], 0.)*temp_gamma_interp[i][j][k][1];
	temp_gamma_interp[i][j][k][4] = pow(gamma_dom[4]/gamma_dom[3], 0.)*temp_gamma_interp[i][j][k][3];
	
	//If the i-front has started to cross, interpolate corrections.  
	if (tcross_flag[i][j][k] == TRUE)  {
		
		//Added 08/19/22: average over tcross bins
		
		for (int tt = 0; tt < Ntcross; tt++)  {
			
			float tc_new_tt = 0.; 
			if (tt == 0)  {
				tc_new_tt = tcross_new[i][j][k][0]/2.;
			}
			else  {
				tc_new_tt = tcross_new[i][j][k][tt - 1] + (tcross_new[i][j][k][tt] - tcross_new[i][j][k][tt - 1])/2.;
			}
		
			for (int n = 0; n < 2; n++)  {
				get_interpolate_gamma_f(i, j, k, tsince - tc_new_tt, tre + tc_new_tt, deltaNL, ind_tre_tcross_new[tt], ind_delta, ind_gamma, mfp_tsince_tcross_new_b[tt][n], mfp_t_tcross_new_gamma_interp_b[tt][n], FALSE);
				get_interpolate_gamma_f(i, j, k, tsince - tc_new_tt, tre + tc_new_tt, deltaNL, ind_tre_tcross_new[tt], ind_delta, ind_gamma, mfp_tsince_tcross_new_f[tt][n], mfp_t_tcross_new_gamma_interp_f[tt][n], FALSE);
			}
		
			get_interpolate_gamma_f(i, j, k, tsince - tc_new_tt, tre + tc_new_tt, deltaNL, ind_tre_tcross_new[tt], ind_delta, ind_gamma, fh1_tsince_tcross_new[tt], fh1_t_tcross_new_gamma_interp[tt], TRUE);
			get_interpolate_gamma_f(i, j, k, tsince - tc_new_tt, tre + tc_new_tt, deltaNL, ind_tre_tcross_new[tt], ind_delta, ind_gamma, temp_tsince_tcross_new[tt], temp_t_tcross_new_gamma_interp[tt], FALSE);
			
			for (int n = 0; n < 2; n++)  {
				mfp_t_tcross_new_gamma_interp_b[tt][n][0] = pow(gamma_dom[0]/gamma_dom[1], 1.)*mfp_t_tcross_new_gamma_interp_b[tt][n][1];
				mfp_t_tcross_new_gamma_interp_f[tt][n][0] = pow(gamma_dom[0]/gamma_dom[1], 1.)*mfp_t_tcross_new_gamma_interp_f[tt][n][1];
				mfp_t_tcross_new_gamma_interp_b[tt][n][4] = pow(gamma_dom[4]/gamma_dom[3], 1.)*mfp_t_tcross_new_gamma_interp_b[tt][n][3];
				mfp_t_tcross_new_gamma_interp_f[tt][n][4] = pow(gamma_dom[4]/gamma_dom[3], 1.)*mfp_t_tcross_new_gamma_interp_f[tt][n][3];
			}
			
			fh1_t_tcross_new_gamma_interp[tt][0] = pow(gamma_dom[0]/gamma_dom[1],  0.)*fh1_t_tcross_new_gamma_interp[tt][1];
			fh1_t_tcross_new_gamma_interp[tt][4] = pow(gamma_dom[4]/gamma_dom[3], -1.)*fh1_t_tcross_new_gamma_interp[tt][3];
	
			temp_t_tcross_new_gamma_interp[tt][0] = pow(gamma_dom[0]/gamma_dom[1], 0.)*temp_t_tcross_new_gamma_interp[tt][1];
			temp_t_tcross_new_gamma_interp[tt][4] = pow(gamma_dom[4]/gamma_dom[3], 0.)*temp_t_tcross_new_gamma_interp[tt][3];
			
		}
		
		for (int q = 0; q < Ngamma+2; q++)  {
			
			//initialize
			for (int n = 0; n < 2; n++)  {
				mfp_gamma_interp_b[n][i][j][k][q]  = 0.; // 1./mfp_gamma_interp_b[n][i][j][k][q];
				mfp_gamma_interp_f[n][i][j][k][q]  = 0.;  // 1./mfp_gamma_interp_f[n][i][j][k][q];
			}
			fh1_gamma_interp[i][j][k][q]       = 0.;
			temp_gamma_interp[i][j][k][q]      = 0.;
			
			//case where the cell is already fully ionized
			if (fion_max[i][j][k] == 1.)  {
				
				for (int tt = 0; tt < Ntcross; tt++)  {
					float tc_frac = 1./((float) Ntcross);
					//sum opacities over all bins, equally weighted
					for (int n = 0; n < 2; n++)  {
						mfp_gamma_interp_b[n][i][j][k][q] += tc_frac/mfp_t_tcross_new_gamma_interp_b[tt][n][q]; //some kind of average
						mfp_gamma_interp_f[n][i][j][k][q] += tc_frac/mfp_t_tcross_new_gamma_interp_f[tt][n][q];
					}
					fh1_gamma_interp[i][j][k][q]  += tc_frac*fh1_t_tcross_new_gamma_interp[tt][q];
					temp_gamma_interp[i][j][k][q] += tc_frac*temp_t_tcross_new_gamma_interp[tt][q];
				}
			}
			//conditions for different stages of partial ionization
			else  {
				for (int tt = 0; tt < Ntcross; tt++)  {
					float tc_frac = 1./((float) Ntcross);
					if ( (fion_max[i][j][k] >= tt*tc_frac) && (fion_max[i][j][k] < (tt+1.)*tc_frac) )  {
						for (int ttt = 0; ttt < tt; ttt++)  {
							for (int n = 0; n < 2; n++)  {
								mfp_gamma_interp_b[n][i][j][k][q] += tc_frac/mfp_t_tcross_new_gamma_interp_b[ttt][n][q];
								mfp_gamma_interp_f[n][i][j][k][q] += tc_frac/mfp_t_tcross_new_gamma_interp_f[ttt][n][q];
							}
							fh1_gamma_interp[i][j][k][q]      += tc_frac*fh1_t_tcross_new_gamma_interp[ttt][q];
							temp_gamma_interp[i][j][k][q]     += tc_frac*temp_t_tcross_new_gamma_interp[ttt][q];
						}
						for (int n = 0; n < 2; n++)  {
							mfp_gamma_interp_b[n][i][j][k][q] += (maxd(0.01,fion_max[i][j][k]) - tt*tc_frac)/mfp_t_tcross_new_gamma_interp_b[tt][n][q];
							mfp_gamma_interp_f[n][i][j][k][q] += (maxd(0.01,fion_max[i][j][k]) - tt*tc_frac)/mfp_t_tcross_new_gamma_interp_f[tt][n][q];
						}
						fh1_gamma_interp[i][j][k][q]      += (maxd(0.01,fion_max[i][j][k]) - tt*tc_frac)*fh1_t_tcross_new_gamma_interp[tt][q];
						temp_gamma_interp[i][j][k][q]     += (maxd(0.01,fion_max[i][j][k]) - tt*tc_frac)*temp_t_tcross_new_gamma_interp[tt][q];
					}
				}
			}
			for (int n = 0; n < 2; n++)  {
				mfp_gamma_interp_b[n][i][j][k][q] = 1./(mfp_gamma_interp_b[n][i][j][k][q]/maxd(0.01,fion_max[i][j][k]));
				mfp_gamma_interp_f[n][i][j][k][q] = 1./(mfp_gamma_interp_f[n][i][j][k][q]/maxd(0.01,fion_max[i][j][k]));
			}
			fh1_gamma_interp[i][j][k][q]      = fh1_gamma_interp[i][j][k][q]/maxd(0.01,fion_max[i][j][k]);
			temp_gamma_interp[i][j][k][q]     = temp_gamma_interp[i][j][k][q]/maxd(0.01,fion_max[i][j][k]);
			fh1_gamma_interp[i][j][k][q]      = mind(fh1_gamma_interp[i][j][k][q], max_subgrid_frac);
		}
	}
}

//Evolution of the MFP. 
//MODIFIED 03/15 to include integration over both frequency averaged MFP and the 912 one
void interp_mfp(int i, int j, int k, float gammaH1_old, float gammaH1_new)  {
	
	float mfp_min  = 1e-2*dx[0];
	float gamma_min_intg = 3e-14;
	float gamma_old = maxd(gammaH1_old, gamma_min);
	float gamma_new = maxd(gammaH1_new, gamma_min);
	
	float ts = 1./(gamma_new + clump[i][j][k]*recomb_H2[i][j][k]*(1. + nhe_to_nh)*pow(nH[i][j][k], 2.));
	
	int gb = (int) floor(((float) Ngamma+2-1.)*(log10(gamma_old) - log10(gamma_dom[0]))/(log10(gamma_dom[Ngamma+2-1]) - log10(gamma_dom[0])));
	int gf = (int) floor(((float) Ngamma+2-1.)*(log10(gamma_new) - log10(gamma_dom[0]))/(log10(gamma_dom[Ngamma+2-1]) - log10(gamma_dom[0])));
	gb = max(0, min(gb, Ngamma+2-2));
	gf = max(0, min(gf, Ngamma+2-2));
	
	//convert to proper units
	for (int n = 0; n < 2; n++)  {
		
		//MODIFIED 05/26/22 to follow the new system for keeping track of redshift
		float log_mfp_tb_gbb_prop = log10(mfp_gamma_interp_b[n][i][j][k][gb]/(Lx0/Lx));
		float log_mfp_tb_gbf_prop = log10(mfp_gamma_interp_b[n][i][j][k][gb+1]/(Lx0/Lx));
		float log_mfp_tf_gbb_prop = log10(mfp_gamma_interp_f[n][i][j][k][gb]/(Lx0/Lx));
		float log_mfp_tf_gbf_prop = log10(mfp_gamma_interp_f[n][i][j][k][gb+1]/(Lx0/Lx));
		
		float log_mfp_tf_gfb_prop = log10(mfp_gamma_interp_f[n][i][j][k][gf]/(Lx0/Lx));
		float log_mfp_tf_gff_prop = log10(mfp_gamma_interp_f[n][i][j][k][gf+1]/(Lx0/Lx));
		
		float log_mfp_tb_gb = log_mfp_tb_gbb_prop + (log_mfp_tb_gbf_prop - log_mfp_tb_gbb_prop)/(log10(gamma_dom[gb+1]) - log10(gamma_dom[gb]))
					 *(log10(gamma_old) - log10(gamma_dom[gb]));
		float log_mfp_tf_gb = log_mfp_tf_gbb_prop + (log_mfp_tf_gbf_prop - log_mfp_tf_gbb_prop)/(log10(gamma_dom[gb+1]) - log10(gamma_dom[gb]))
					 *(log10(gamma_old) - log10(gamma_dom[gb]));
		float log_mfp_tf_gf = log_mfp_tf_gfb_prop + (log_mfp_tf_gff_prop - log_mfp_tf_gfb_prop)/(log10(gamma_dom[gf+1]) - log10(gamma_dom[gf]))
					 *(log10(gamma_new) - log10(gamma_dom[gf]));
		
		float delta_mfp_t = pow(10., log_mfp_tf_gb) - pow(10., log_mfp_tb_gb);
		
		float delta_log_mfp_gamma = alpha_mfp*(log10(gamma_new) - log10(gamma_old));
		
		float delta_log_mfp = delta_log_mfp_gamma;
		
		//MODIFIED 03/06 to ensure that cells that are ionized in one time step do not end up with pathological initial conditions.    Also
		//modified on 03/06 to ensure that the integration is not carried out over the noisy part of the solution.  
		//MODIFIED 03/15 to do both MFP's
		if (n == 0)  { //start with the frequency averaged MFP
			
			float delta_mfp_relax = -1./(t_relax*1e6*yr_to_s)*(mfp_prev[i][j][k] - pow(10., log_mfp_tf_gf))*dt;
			
			if ( (fion[i][j][k] < 1.) || (gamma_old < gamma_min_intg) || (gamma_new < gamma_min_intg) )  { //don't integrate
				mfp[i][j][k]      = pow(10., log_mfp_tf_gf);
			}
			else  { //integrate
				mfp[i][j][k] = log10(mfp_prev[i][j][k]) + delta_log_mfp;
				mfp[i][j][k] = maxd(pow(10., mfp[i][j][k]) + delta_mfp_relax + delta_mfp_t, mfp_min); 
			}
			//average kappa over the photo-ionization + recombination timescale
			
			mfp[i][j][k] = maxd(mfp[i][j][k], mfp_min);
		}
		else  { //now do the same thing to the 912 one
			
			float delta_mfp_relax = -1./(t_relax*1e6*yr_to_s)*(mfp_912_prev[i][j][k] - pow(10., log_mfp_tf_gf))*dt;
			
			if ( (fion[i][j][k] < 1.) || (gamma_old < gamma_min_intg) || (gamma_new < gamma_min_intg) )  {
				mfp_912[i][j][k]      = pow(10., log_mfp_tf_gf);
			}
			else  {
				mfp_912[i][j][k] = log10(mfp_912_prev[i][j][k]) + delta_log_mfp;
				mfp_912[i][j][k] = maxd(pow(10., mfp_912[i][j][k]) + delta_mfp_relax + delta_mfp_t, mfp_min); 
			}
			//average kappa over the photo-ionization + recombination timescale
			mfp_912[i][j][k] = maxd(mfp_912[i][j][k], mfp_min);
		}
	}
}

//Interpolate over the total HI fraction, same procedure as ionized HI.  
void interp_nh1(int i, int j, int k, float gammaH1)  {
	float gamma = maxd(gammaH1, gamma_min);
	float ts = 1./(gamma + clump[i][j][k]*recomb_H2[i][j][k]*ne[i][j][k]);
	
	int g = (int) floor(((float) Ngamma+2-1.)*(log10(gamma) - log10(gamma_dom[0]))/(log10(gamma_dom[Ngamma+2-1]) - log10(gamma_dom[0])));
	g = max(0, min(g, Ngamma+2-2));

	nH1[i][j][k] = log10(fh1_gamma_interp[i][j][k][g]) + (log10(fh1_gamma_interp[i][j][k][g+1]) - log10(fh1_gamma_interp[i][j][k][g]))/(log10(gamma_dom[g+1]) - log10(gamma_dom[g]))
				 *(log10(gamma) - log10(gamma_dom[g]));
	nH1[i][j][k] = maxd(mind(pow(10., nH1[i][j][k]), fH1_0), 0.);
	nH1[i][j][k] *= nH[i][j][k];
}

//Interpolate xHeI.  Assume it is equal to xHI.  
void interp_nhe1(int i, int j, int k, float gammaH1)  {
	nHe1[i][j][k]     = nH1[i][j][k]*nHe[i][j][k]/nH[i][j][k]; //Helium = Hydrogen
	avg_nhe1[i][j][k] = nHe1[i][j][k];
}

//Log-log interpolation to get temperature.  
//Rewrite this interpolation using the other approach above
void interp_temp(int i, int j, int k, float gammaH1)  {
	float Tmin = 1e2;
	float gamma = maxd(gammaH1, gamma_min);
	float log_gamma_dom[Ngamma+2];
	float log_temp_gamma_interp[Ngamma+2];
	for (int g = 0; g < Ngamma+2; g++)  {
		log_gamma_dom[g] = log10(gamma_dom[g]);
		log_temp_gamma_interp[g] = log10(temp_gamma_interp[i][j][k][g]);
	}
	temp[i][j][k] = interpolate(log_gamma_dom, log_temp_gamma_interp, log10(gamma), Ngamma+2);
	temp[i][j][k] = maxd(pow(10., temp[i][j][k]), Tmin);
}

//ADDED 03/04: calculate the Reion temperature Treion as a function of I-front velocity (D'Aloisio+18)
void set_Treion(float vIF, float fi_old, float fi_new, int i, int j, int k)  { //set Treion (vIF should be input in cm/s)
	float Cn[5] = {9.5432, -1.6441e-2, -9.8010e-3, 4.1664e-3, -2.2710e-4};
	float lnTreion = 0.;
	for (int t = 0; t < 5; t++)  {
		lnTreion += Cn[t]*pow(log(vIF/1e5), t); //Treion from Daloisio+18, in K
	}
	if (fi_new > 0.)  {
		temp[i][j][k] = (fi_old*temp[i][j][k] + (fi_new - fi_old)*exp(lnTreion))/fi_new;
	}
	else  {
		temp[i][j][k] = temp[i][j][k];
	}
	
	//ADDED 05/27/22 - Fahad's update to the initial heating term for MP-Gadget
	if (MP_GADGET_HEAT == TRUE)  {
		if ( fi_new > fi_old ){
			//(mu_ion)-1 = SUM (fraction_of_ion/Atomic_Mass_of_ion)
			//(mu_electron)-1 = SUM (fraction_of_ion ionized_fraction_of_ion Atomic_charge /Atomic_Mass_of_ion)

			//completely ionized hydrogen and heium
			//(find mu for ions and electrons separately) - 2/(3(1-Y)+1) -   0.6131207847946045
			float mu = 2/( 3*(1-Y)+1 );

			// erg/grams
			//PLEASE DONT ADD HERE ACCUMULATIVE TOTAL HEATING IS ADDING AT THE END OF EACH STEP
			gadget_ini_heat[i][j][k] = 1.5 * (fi_new-fi_old) * k_B*exp(lnTreion) / (mu * m_H);
		}
	}
}

//MODIFIED 03/04 to include the Treion update and 03/06 to fix potential issue with variable types
//update the fraction of a cell that has been passed by i-fronts and return the left-over time if there is any.  
//Assume I-front speed from D'Aloisio+18 and assume all the number of photons absorbed = number of HI atoms in the 
//ionized region.  
void update_fion(int i, int j, int k)  {
	float tdiff;
	float vIF_min = 5e1*1e5;
	double F = ngam[i][j][k][2] / pow(dx[0], 2.) / dt;
	float vIF = (float) (F / nH[i][j][k] / (1. + nhe_to_nh)); 
	
	float fion_old     = fion[i][j][k];
	float fion_max_old = fion_max[i][j][k];
	fion[i][j][k] = mind(fion[i][j][k] + vIF*dt/dx[0], 1.); //update fion
	
	//Added 08/23/22: keep track of maximum ionization fraction
	if (fion[i][j][k] > fion_max[i][j][k])  {
		fion_max[i][j][k] = fion[i][j][k];
	}
	
	if ( (vIF > vIF_min) || (fion_old == 0.) )  {
		set_Treion(maxd(vIF, vIF_min), fion_old, fion[i][j][k], i, j, k); //Append the temperature of the newly ionized part to the running average
	}
	
	if ( (fion_old < 0.5) && (fion[i][j][k] >= 0.5) )  {
		treion_out[i][j][k] = t - dt;
	}
	
	if (fion[i][j][k] >= 1.)  { //leftover time if fion > 1
		tdiff = 0.5*dt; //true on average
	}
	else  {
		tdiff = 0.;
	}
	
	//Added 08/31/22 - keep track of IF speed and incident flux
	IF_speed[i][j][k] = vIF;
	Inc_Flux[i][j][k] = F;
	Delta_fion[i][j][k] = fion[i][j][k] - fion_old;

	//update crossing time (will add this back in shortly)
	tcross[i][j][k] += dt - tdiff;
	
	//Added 08/19/22: keep track of the I-front crossing times at equally spaced bins in fion.  This will correctly keep track
	//of the crossing times if the I-front crosses multiple fion bins during a time step.  
	for (int tt = 0; tt < Ntcross; tt++)  {
		float tc_frac = 1./((float) Ntcross);
		if ( (fion_max[i][j][k] < tc_frac*(tt + 1.)) && (tcross_check[i][j][k][tt] == FALSE) )  {
			tcross_new[i][j][k][tt] += dt; 
		}
		else if ( (fion_max[i][j][k] >= tc_frac*(tt + 1.)) && (fion_max_old < tc_frac*(tt + 1.)) && (tcross_check[i][j][k][tt] == FALSE) )  {
			tcross_new[i][j][k][tt] += dt*(tc_frac*(tt + 1.) - fion_max_old)/(fion_max[i][j][k] - fion_max_old);
			tcross_check[i][j][k][tt] = TRUE;
		} 
		else  {
			tcross_new[i][j][k][tt] = tcross_new[i][j][k][tt];
		}
	}
	
	if (tcross_flag[i][j][k] == FALSE)  {
		tcross_flag[i][j][k] = TRUE;
	}
}

//Self-shielding threshold for the photo-ionization rate. Motivated by Nasir+21 threshold.  
float gamma_ss(int i, int j, int k)  {
	double chi = nhe_to_nh;
	return (float) (1./xh1_thresh*(1. + chi)*alphaB_H2(1e4)*nH[i][j][k]*pow(1. - xh1_thresh, 2.));
}

//Import the Mao+19 clumping model.  This consists of a set of redshift-dependent fit parameters and a redshift grid.  
void get_mao_2019_clumping(void)  {
	double temp_double;
	const char mao_clumping_file[] = "subgrid_files/mao_2019_clumping.txt";
	
	FILE *file = NULL;
	file = fopen(mao_clumping_file, "r");
	if (file == NULL)  {
		printf("Mao 2019 clumping file not found\n");
	}
	else  {
		for (int i = 0; i < n_mao; i++)  {
			fscanf(file, "%le\n", &temp_double);
			mao_z[i]  = (float) temp_double;
			fscanf(file, "%le\n", &temp_double);
			mao_a0[i] = (float) temp_double;
			fscanf(file, "%le\n", &temp_double);
			mao_a1[i] = (float) temp_double;
			fscanf(file, "%le\n", &temp_double);
			mao_a2[i] = (float) temp_double;
		}
	}
}

//Interpolate in redshift over the Mao+19 fit parameter asn construct the clumping factor as a function over overdensity.  
float interp_mao_clumping(float redshift, float delta)  {
	float x = log10(delta);
	float a0 = interpolate(mao_z, mao_a0, redshift, n_mao);
	float a1 = interpolate(mao_z, mao_a1, redshift, n_mao);
	float a2 = interpolate(mao_z, mao_a2, redshift, n_mao);
	if (redshift <= mao_z[4*(n_mao - 1)])  {
		a0 = mao_a0[4*(n_mao - 1) + 1];
		a1 = mao_a1[4*(n_mao - 1) + 2];
		a2 = mao_a2[4*(n_mao - 1) + 3];
	}
	float y = a0 + a1*x + a2*pow(x, 2.);
	return pow(10., y);
}

//Set a constant clumping factor if desired.  
void set_clump(void)  {
	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				clump[i][j][k] = const_clump;
			}
		}
	}
	}
}
