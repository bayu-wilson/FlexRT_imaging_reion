//Input/output functions
#include <omp.h>
#include <fstream>
#include <stdio.h>
#include <dirent.h> 
#include <string.h>
#include "data_funcs.cc"
using namespace std;

float sum_nphot(long int x)  {
	float nphot = 0.;
	for (int nu = 0; nu < Nfreq; nu++)  {
		nphot += ray[x].nphot[nu];
	}
	return nphot;
}

double count_photons_io(void)  {
	double photon_counter = 0.;
	#pragma omp parallel reduction (+ : photon_counter)
	{
	#pragma omp for
	for (long int x = 0; x < size_buffer; x++)  {
		if (ray_tags[x] == TRUE)  {
			photon_counter += (double) sum_nphot(x);
		}
	}
	}
	return photon_counter*phot_norm;
}

//Read the gas density from one of Hy's gas files
void get_ugas(const char *ugas_file)  {
	ifstream file (ugas_file, ios::in | ios::binary);

	float  temp_float;
	float om = 0.305147;
	float ol = 1 - om;
   	float ob = 0.0482266;
    float h0 = 0.68;
	float rho_c_cgs = 1.87890e-29;
	avg_rho = 0.;
	//get the conversion from Hy's grid units to the 
	//MODIFIED 05/26/22 to remove conversion to proper units from here (since it is now in update_hydro)s
	double rho_unit = rho_c_cgs*ob*pow(h0, 2.); // * pow(1. + zz, 3.);
	
	//Convert from fortran order.  Don't change loop order!!! 
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				file.read((char*)&temp_float, sizeof(float));
				rho[i][j][k] = temp_float;
				rho[i][j][k] *= (float) rho_unit; //convert to physical units
				//file.read((char*)&temp_float, sizeof(float)); //skip past the other four quantities in Hy's output file
				//file.read((char*)&temp_float, sizeof(float));
				//file.read((char*)&temp_float, sizeof(float));
				//file.read((char*)&temp_float, sizeof(float)); //I could eventually get the temperture from here if I wanted to
			}
		}
	}
	avg_rho = calc_vol_avg(rho);
}

//Read the treion file for a restart.  
void read_treion(const char *treion_file)  {
	ifstream file (treion_file, ios::in | ios::binary);
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				file.read((char*)&treion[i][j][k], sizeof(float));
				file.read((char*)&treion_flag[i][j][k], sizeof(bool));
				file.read((char*)&tcross_flag[i][j][k], sizeof(bool));
				file.read((char*)&neutral_flag[i][j][k], sizeof(bool));
			}
		}
	}
}

//MODIFIED 04/08/21 to initialize previous MFP values on restart
//Read the gas grid for a restart.  
void read_grid_binary(const char *start_file)  {
	ifstream file (start_file, ios::in | ios::binary);
	float temp_float;
	file.read((char*)&t, sizeof(float));
	t_tot = t; //total simulation runtime
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				file.read((char*)&rho[i][j][k], sizeof(float));
				file.read((char*)&temp[i][j][k], sizeof(float));
				file.read((char*)&nH1[i][j][k], sizeof(float));
				file.read((char*)&f_H1[i][j][k], sizeof(float));
				file.read((char*)&f_He1[i][j][k], sizeof(float));
				file.read((char*)&f_H2[i][j][k], sizeof(float));
				file.read((char*)&f_He2[i][j][k], sizeof(float));
				file.read((char*)&f_He3[i][j][k], sizeof(float));
				file.read((char*)&gamma_H1_tot[i][j][k], sizeof(float));
				file.read((char*)&gamma_He1_tot[i][j][k], sizeof(float));
				file.read((char*)&mfp[i][j][k], sizeof(float));
				file.read((char*)&mfp_912[i][j][k], sizeof(float));
				file.read((char*)&heat_rate[i][j][k], sizeof(float));
				file.read((char*)&cool_rate[i][j][k], sizeof(float));
				file.read((char*)&avg_nh1[i][j][k], sizeof(float));
				file.read((char*)&fion[i][j][k], sizeof(float));
				file.read((char*)&tcross[i][j][k], sizeof(float));
				file.read((char*)&treion_out[i][j][k], sizeof(float));
				for (int nu = 0; nu < Nfreq; nu++) {
          				file.read((char*)&ngam_nu[i][j][k][nu], sizeof(float));
				}
				mfp_prev[i][j][k] = mfp[i][j][k];
				mfp_912_prev[i][j][k] = mfp_912[i][j][k];
			}
		}
	}
	file.close();
}

//Read the rays from a save file for a restart.  
void read_rays_binary(const char *ray_file)  {
	long int x = 0;
	ifstream file (ray_file, ios::in | ios::binary);
	file.read((char*)&x, sizeof(long int));
	while (x != -1)  {
		ray_tags[x] = TRUE;
		file.read((char*)&ray[x].active, sizeof(bool));
		file.read((char*)&ray[x].order, sizeof(short int));
		file.read((char*)&ray[x].pixel, sizeof(int));
		file.read((char*)&ray[x].itheta, sizeof(unsigned short int));
		file.read((char*)&ray[x].xhat, sizeof(float));
		file.read((char*)&ray[x].yhat, sizeof(float));
		file.read((char*)&ray[x].zhat, sizeof(float));
		for (int nu = 0; nu < Nfreq; nu++)  {
			file.read((char*)&ray[x].nphot[nu], sizeof(float));
		}
		file.read((char*)&ray[x].r, sizeof(float));
		file.read((char*)&ray[x].dr, sizeof(float));
		file.read((char*)&ray[x].x, sizeof(float));
		file.read((char*)&ray[x].y, sizeof(float));
		file.read((char*)&ray[x].z, sizeof(float));
		file.read((char*)&ray[x].time, sizeof(float));
		file.read((char*)&ray[x].path[0].ind, sizeof(int));
		file.read((char*)&x, sizeof(long int));
	}
	file.close();
}

//calculate the fractional cross-section of a source within it's cell
void update_source_opacity(int i, int j, int k, double M200)  {
	double r200 = R200(M200, zz);
	double exp_cc = 3./2.;
	double cell_area = pow(Lx/Nx*kpc_to_cm, 2.);
	source_tau[i][j][k] += 2.*pi*pow(r200, 2.)/3./cell_area;
}

//Read the sources from the catalog and convert to luminosity
//Read the sources from the catalog and convert to luminosity
void read_source_catalog(const char *source_field, const char *num_src_file)  {
        ifstream file (source_field, ios::in | ios::binary);
        ifstream file2;
        int CATDIM  = 5; //for reading ndot with mass in there
        long int temp_int;
        double x,y,z,massH,LumGAL;

        file2.open(num_src_file, ios::in | ios::binary);
        file2.read((char*)&temp_int, sizeof(long int));
        num_src = (int) (temp_int);
        num_inc = (int) (temp_int);

        printf("%d\n", num_src);
        printf("%d\n", temp_int);

        double source_grid[Nx][Ny][Nz];
        //initialize source grid
        #pragma omp parallel
        {
        #pragma omp for 
        for (int i = 0; i < Nx; i++)  {
                for (int j = 0; j < Nx; j++)   {
                        for (int k = 0; k < Nz; k++)   {
                                source_grid[i][j][k] = 0.;
                                source_tau[i][j][k] = 0.;
                        }
                }
        }
        }

        //resize the source luminosity vector
        source_lum.resize(min(num_src, Nx*Ny*Nz));

		for (int si = 0; si < num_inc; si++)  {
				file.read((char*)&x, sizeof(double));
				file.read((char*)&y, sizeof(double));
				file.read((char*)&z, sizeof(double));
				file.read((char*)&massH, sizeof(double));
				file.read((char*)&LumGAL, sizeof(double));
				int ix   = (int) floor(x/hh/Lx0*1e3*Nx);
				int jx   = (int) floor(y/hh/Ly0*1e3*Ny);
				int kx   = (int) floor(z/hh/Lz0*1e3*Nz);
				ix = mod(ix, Nx); //max(0, min(ix, Nx-1));
				jx = mod(jx, Ny); //max(0, min(jx, Ny-1));
				kx = mod(kx, Nz); //max(0, min(kx, Nz-1));
				//use formula from Scorch II to get ionizing photon counts
				source_grid[ix][jx][kx] += fesc*LumGAL; //WHAT UNITS?
				
				//update_source_opacity(ix, jx, kx, massH); //update source opacity
		}

        //assign the gridded sources to the source luminosity vector
        int counter = 0;
        for (int i = 0; i < Nx; i++)  {
                for (int j = 0; j < Nx; j++)   {
                        for (int k = 0; k < Nz; k++)   {
                                if (source_grid[i][j][k] != 0.)  {
                                        source_lum[counter].i = i;
                                        source_lum[counter].j = j;
                                        source_lum[counter].k = k;
                                        source_lum[counter].lum = source_grid[i][j][k];
                                        counter += 1;
                                        ray_dt[i][j][k] = 0.;
                                }
                                else  {
								//MODIFIED 07/26/22 to avoid double counting of ray times in casting routine
                                        ray_dt[i][j][k] = 0.;
                                }
								source_tau[i][j][k] = -1.*log(1. - source_tau[i][j][k]); //source opacity
                        }
                }
        }

        num_src = counter;
        num_inc = counter;
        source_lum.resize(num_src);
}

//Read the hydro step redshifts from the input file.  
void read_hydro_steps(void)  {
        int temp_int;
        double temp_double;
        FILE *file = NULL;
        file = fopen(hydro_steps_file, "r");
        if (file == NULL)  {
                printf("Hydro steps file not found\n");
        }
        else  {
                printf("Reading hydro steps\n");
                fscanf(file, "%d\n", &temp_int);
                num_hydro_steps = (int) temp_int;
                printf("Number of hydro steps: %d\n", num_hydro_steps - 1);
                hydro_steps.resize(num_hydro_steps, 0.);
                for (int i = 0; i < num_hydro_steps; i++)  {
                        fscanf(file, "%le\n", &temp_double);
                        fscanf(file, "%le\n", &temp_double);
                        fscanf(file, "%le\n", &temp_double);
                        hydro_steps[i] = (float) temp_double;
                        printf("hydro step %d : %le\n", i, hydro_steps[i]);
                }
                zz_end = hydro_steps[num_hydro_steps - 1];
                printf("Hydro steps read\n");
        }
}

//Create the output files for the run.  
void make_output(void)  {
	FILE *file = NULL;
	//printf("please work 3\n");
        //fflush(stdout);	
	file = fopen(otf_output, "w");
	fprintf(file, "On-the-fly output\n");
	fprintf(file, "\n");
	fclose(file);
}

//Write the rays to save file
void write_rays_binary(const char *ray_output)  {
	ofstream file;
	long int end = -1;
	printf("Writing ray file\n");
	file.open(ray_output, ios::out | ios::binary);
	for (long int x = 0; x < size_buffer; x++)  {
		if (ray_tags[x] == TRUE)  {
			file.write((char*)&x, sizeof(long int));
			file.write((char*)&ray[x].active, sizeof(bool));
			file.write((char*)&ray[x].order, sizeof(short int));
			file.write((char*)&ray[x].pixel, sizeof(int));
			file.write((char*)&ray[x].itheta, sizeof(unsigned short int));
			file.write((char*)&ray[x].xhat, sizeof(float));
			file.write((char*)&ray[x].yhat, sizeof(float));
			file.write((char*)&ray[x].zhat, sizeof(float));
			for (int nu = 0; nu < Nfreq; nu++)  {
				file.write((char*)&ray[x].nphot[nu], sizeof(float));
			}
			file.write((char*)&ray[x].r, sizeof(float));
			file.write((char*)&ray[x].dr, sizeof(float));
			file.write((char*)&ray[x].x, sizeof(float));
			file.write((char*)&ray[x].y, sizeof(float));
			file.write((char*)&ray[x].z, sizeof(float));
			file.write((char*)&ray[x].time, sizeof(float));
			file.write((char*)&ray[x].path[0].ind, sizeof(int));
		}
	}
	file.write((char*)&end, sizeof(long int));
	file.close();
}

//ADDED 06/02/22 to write a reduced-size ray file for testing purposes
void write_reduced_rays_binary(const char *ray_output)  {
	ofstream file;
	long int end = -1;
	printf("Writing ray file\n");
	file.open(ray_output, ios::out | ios::binary);
	for (long int x = 0; x < size_buffer; x++)  {
		if (ray_tags[x] == TRUE)  {
			//file.write((char*)&x, sizeof(long int));
			//file.write((char*)&ray[x].active, sizeof(bool));
			//file.write((char*)&ray[x].order, sizeof(short int));
			//file.write((char*)&ray[x].pixel, sizeof(int));
			//file.write((char*)&ray[x].itheta, sizeof(unsigned short int));
			file.write((char*)&ray[x].xhat, sizeof(float));
			file.write((char*)&ray[x].yhat, sizeof(float));
			file.write((char*)&ray[x].zhat, sizeof(float));
			//for (int nu = 0; nu < Nfreq; nu++)  {
			//	file.write((char*)&ray[x].nphot[nu], sizeof(float));
			//}
			//file.write((char*)&ray[x].r, sizeof(float));
			//file.write((char*)&ray[x].dr, sizeof(float));
			file.write((char*)&ray[x].x, sizeof(float));
			file.write((char*)&ray[x].y, sizeof(float));
			file.write((char*)&ray[x].z, sizeof(float));
			//file.write((char*)&ray[x].time, sizeof(float));
			file.write((char*)&ray[x].path[0].ind, sizeof(float));
		}
	}
	file.write((char*)&end, sizeof(long int));
	file.close();
}

//Write the treion grid to save file
void write_treion(const char *tre_output)  {
	ofstream file;
	
	printf("Writing treion file\n");
	file.open(tre_output, ios::out | ios::binary);
	
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				file.write((char*)&treion[i][j][k], sizeof(float));
				file.write((char*)&treion_flag[i][j][k], sizeof(bool));
				file.write((char*)&tcross_flag[i][j][k], sizeof(bool));
				file.write((char*)&neutral_flag[i][j][k], sizeof(bool));
			}
		}
	}
	file.close();
}

//Write the gas grid to save file
void write_gas_binary(const char *gas_file)  {
	
	ofstream file;
	printf("Writing gas file\n");
	file.open(gas_file, ios::out | ios::binary);
	file.write((char*)&t, sizeof(float));
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				file.write((char*)&rho  [i][j][k],         sizeof(float)); //0
				file.write((char*)&temp [i][j][k],         sizeof(float)); //1
				file.write((char*)&nH1  [i][j][k],         sizeof(float)); //2
				file.write((char*)&f_H1 [i][j][k],         sizeof(float)); //3
				file.write((char*)&f_He1[i][j][k],         sizeof(float)); //4
				file.write((char*)&f_H2 [i][j][k],         sizeof(float)); //5
				file.write((char*)&f_He2[i][j][k],         sizeof(float)); //6
				file.write((char*)&f_He3[i][j][k],         sizeof(float)); //7
				file.write((char*)&gamma_H1_tot[i][j][k],  sizeof(float)); //8
				file.write((char*)&gamma_He1_tot[i][j][k], sizeof(float)); //9
				file.write((char*)&mfp[i][j][k],           sizeof(float)); //10
				file.write((char*)&mfp_912[i][j][k],       sizeof(float)); //11
				file.write((char*)&heat_rate[i][j][k],     sizeof(float)); //12
				file.write((char*)&cool_rate[i][j][k],     sizeof(float)); //13
				file.write((char*)&avg_nh1[i][j][k],       sizeof(float)); //14
				file.write((char*)&fion[i][j][k],          sizeof(float)); //15
				file.write((char*)&tcross[i][j][k],        sizeof(float)); //16
				file.write((char*)&treion_out[i][j][k],    sizeof(float)); //17
				//Added by Chris 01/03/22: write gamma_nu to gas file
				for (int nu = 0; nu < Nfreq; nu++)  {
					file.write((char*)&ngam_nu[i][j][k][nu], sizeof(float));
				}
			}
		}
	}
	file.close();
}

//Added 08/31/22 - write vIF file
void write_vIF_binary()  {
	string s1,vIF_out;
	
	float a_rt = scale_factor(t + cosmic_time(1./(1.+zinit))); //ADDED 08/31/22
        float z_rt = 1./a_rt - 1.;

	string zz_string = to_string(z_rt);
	
	if (zz < 10.)  {
		s1      = std::string(output_dir)+"vIF_z=0";
		vIF_out = s1 + zz_string;
	}
	else  {
		s1      = std::string(output_dir)+"vIF_z=";
		vIF_out = s1 + zz_string;
	}

	vIF_out.resize(((int) vIF_out.length()) - 2);
	const char* vIF_file = vIF_out.c_str();
	
	ofstream file;
	printf("Writing gas file\n");
	file.open(vIF_file, ios::out | ios::binary);
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				file.write((char*)&IF_speed[i][j][k], sizeof(float)); 
				file.write((char*)&Inc_Flux[i][j][k], sizeof(float));
			        file.write((char*)&fion[i][j][k],sizeof(float));	
				file.write((char*)&Delta_fion[i][j][k],sizeof(float));
			}
		}
	}
	file.close();
}
	

//MODIFIED 05/27/22 - update the MP-Gadget variables
void update_gadget_var()  {
	
	gadget_step_time += dt;

	//get average for the arrays
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			for (int k = 0; k < Nz; k++) {
				//count the total photoionzations in each step
				gadget_gamma[i][j][k]  += (gamma_H1_tot[i][j][k]*dt);
				gadget_heat[i][j][k] += (gadget_ini_heat[i][j][k] + heat_rate[i][j][k]*dt/rho[i][j][k]);
			}
		}
	}
}

//MODIFIED 05/27/22 - Write the smaller gas grid on every RT step
void write_gadget_binary(float current_time)  {

	//get average for heating and photoionization
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			for (int k = 0; k < Nz; k++){
				// in 1/sec
				gadget_gamma[i][j][k] /= gadget_step_time;
				// in erg/sec/grams
				gadget_heat[i][j][k] /= gadget_step_time;
			}
		}
	}

	ofstream file;
	// the file time is the current time minus half-way in each gadget time step plus
	string small_gas_file = "output_files/small_gas_rt_z=" +
	to_string(toRedshift(current_time - gadget_step_time/2.0));

	printf("Writing gas file %s\n", small_gas_file);
	file.open(small_gas_file.c_str(), ios::out | ios::binary);
	file.write((char*)gadget_gamma, Nx*Ny*Nz*sizeof(float));
	// this heating rate plus the initial heat injection due to cell ionized to Treion.
	file.write((char*)gadget_heat, Nx*Ny*Nz*sizeof(float));
	file.close();
}


//update the gadget variables
//ADDED 05/27/22 - reset MP-Gadget heating variables
void reset_gadget_var()  {

	gadget_step_time = 0.0;

	//reset the gadget arrays to zero for the next cycle
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			for (int k = 0; k < Nz; k++) {
				gadget_gamma[i][j][k] = 0.0;
				gadget_heat[i][j][k] = 0.0;
				gadget_ini_heat[i][j][k] = 0.0;
			}
		}
	}
	
}

//Write the on the fly data
//MODIFIED 03/06 to include averages over binary ionization field (fully ionized cells only)
void write_otf(void)  {
	FILE *file = NULL;
	
	float gamma_min_intg = 3e-14;
	float inH1[Nx][Ny][Nz];
	float inH[Nx][Ny][Nz];
	float binary[Nx][Ny][Nz];
	float rho_mean[Nx][Ny][Nz];
	float temp0[Nx][Ny][Nz];
	

        printf("please work write_otf()\n");
        fflush(stdout);
	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				if (fion[i][j][k] == 1.)  {
					inH1[i][j][k] = nH1[i][j][k];
					inH[i][j][k] = nH[i][j][k];
				}
				else  {
					inH1[i][j][k] = 0.;
					inH[i][j][k] = 0.;
				}
				
				if ( (fion[i][j][k] == 1.) && (avg_gamma[i][j][k] >= gamma_min_intg) && (avg_gamma_prev[i][j][k] >= gamma_min_intg)
						&& (abs(log10(avg_gamma[i][j][k]) - log10(avg_gamma_prev[i][j][k])) < 1.) )  {
					binary[i][j][k] = 1.;
				}
				else  {
					binary[i][j][k] = 0.;
				}
				
				//Fixed temp @ mean density averaging in no subgrid case on 09/27/21
				if ( (rho[i][j][k]/avg_rho > 0.9) && (rho[i][j][k]/avg_rho < 1.1) )  {
					if (subgrid == TRUE)  {
						rho_mean[i][j][k] = 1.;
						temp0[i][j][k] = temp[i][j][k]*fion[i][j][k] + temp_0*(1. - fion[i][j][k]);
					}
					else  {
						rho_mean[i][j][k] = 1.;
						temp0[i][j][k] = temp[i][j][k]*f_H2[i][j][k] + temp_0*(1. - f_H2[i][j][k]);
					}
 				}
				else  {
					rho_mean[i][j][k] = 0.;
					temp0[i][j][k] = temp_0;
				}
				
			}
		}
	}
	}

	file = fopen(otf_output, "a");

	fprintf(file, "Step Number    : %d\n", step);
	fprintf(file, "Redshift       : %le\n", zz);
	fprintf(file, "Time [Myr]     : %le\n", t/yr_to_s/1e6);
	fprintf(file, "Time Step [Myr]: %le\n", dt/yr_to_s/1e6);
	fprintf(file, "\n");
    
	float fH1_avg_vol   = calc_vol_avg(f_H1);
	float fH2_avg_vol   = calc_vol_avg(f_H2);
	float fHe1_avg_vol  = calc_vol_avg(f_He1);
	float fHe2_avg_vol  = calc_vol_avg(f_He2);
	float fHe3_avg_vol  = calc_vol_avg(f_He3);
	
	float fH1_avg_mass  = calc_mass_avg(f_H1, nH);
	float fH2_avg_mass  = calc_mass_avg(f_H2, nH);
	float fHe1_avg_mass = calc_mass_avg(f_He1, nHe);
	float fHe2_avg_mass = calc_mass_avg(f_He2, nHe);
	float fHe3_avg_mass = calc_mass_avg(f_He3, nHe);
	
	float nh1_avg  = calc_vol_avg(avg_nh1);
	float nhe1_avg = calc_vol_avg(avg_nhe1);
	
	float nH_avg  = calc_vol_avg(nH);
    float nHe_avg = calc_vol_avg(nHe);
    float ne_avg  = calc_vol_avg(ne);
	
	float nH_avg_fion  = calc_mass_avg(nH, fion);
    float nHe_avg_fion = calc_mass_avg(nHe, fion);
	float nH_avg_binary  = calc_mass_avg(nH, binary);
    float nHe_avg_binary = calc_mass_avg(nHe, binary);

    float gam_h1_vol   	  = calc_vol_avg(gamma_H1_tot);
    float gam_h1_vol_fion = calc_mass_avg(gamma_H1_tot, fion);
	float gam_h1_vol_binary = calc_mass_avg(gamma_H1_tot, binary);
	float gam_h1_mass1 	  = calc_mass_avg(gamma_H1_tot, inH);
    float gam_h1_mass2 	  = calc_mass_avg(gamma_H1_tot, nH1);
	
	float rec_h2_vol      = calc_vol_avg(recrate_H2);
	float rec_h2_vol_fion = calc_mass_avg(recrate_H2, fion);
	float rec_h2_vol_binary = calc_mass_avg(recrate_H2, binary);
	float rec_h2_mass1    = calc_mass_avg(recrate_H2, inH);
    float rec_h2_mass2    = calc_mass_avg(recrate_H2, nH2);
	
	float temp_vol      = calc_vol_avg(temp);
	float temp_vol_fion = calc_mass_avg(temp, fion);
	float temp_vol_binary = calc_mass_avg(temp, binary);
	float temp_mass     = calc_mass_avg(temp, nH);
	float temp_mean_den = calc_mass_avg(temp0, rho_mean);
	
	float mfp_vol      = calc_vol_avg(mfp);
	float mfp_vol_fion = calc_mass_avg(mfp, fion);
	float mfp_vol_binary = calc_mass_avg(mfp, binary);
	float mfp_mass1    = calc_mass_avg(mfp, inH);
	
	float mfp_912_vol        = calc_vol_avg(mfp_912);
	float mfp_912_vol_fion   = calc_mass_avg(mfp_912, fion);
	float mfp_912_vol_binary = calc_mass_avg(mfp_912, binary);
	float mfp_912_mass1      = calc_mass_avg(mfp_912, inH);
	
	float heat_vol           = calc_vol_avg(heat_rate);
	float heat_vol_fion      = calc_mass_avg(heat_rate, fion);
	float heat_vol_binary    = calc_mass_avg(heat_rate, binary);
	float heat_mass1         = calc_mass_avg(heat_rate, nH);
	float heat_mass2         = calc_mass_avg(heat_rate, inH);
	
	float cool_vol           = calc_vol_avg(cool_rate);
	float cool_vol_fion      = calc_mass_avg(cool_rate, fion);
	float cool_vol_binary    = calc_mass_avg(cool_rate, binary);
	float cool_mass1         = calc_mass_avg(cool_rate, nH);
	float cool_mass2         = calc_mass_avg(cool_rate, inH);
	
	double photon_count = count_photons_io();
	float fion_avg_vol  = calc_vol_avg(fion);
	float fion_avg_mass = calc_mass_avg(fion, nH);

	fprintf(file, "Global Averages\n");
	fprintf(file, "fH1     : %le    %le\n", fH1_avg_vol, fH1_avg_mass);
	fprintf(file, "fH2     : %le    %le\n", fH2_avg_vol, fH2_avg_mass);
	fprintf(file, "fHe1    : %le    %le\n", fHe1_avg_vol, fHe1_avg_mass);
	fprintf(file, "fHe2    : %le    %le\n", fHe2_avg_vol, fHe2_avg_mass);
	fprintf(file, "fHe3    : %le    %le\n", fHe3_avg_vol, fHe3_avg_mass);
	fprintf(file, "nH1     : %le       \n", nh1_avg); //these are the time-averaged quantities
	fprintf(file, "nHe1    : %le       \n", nhe1_avg); //
	fprintf(file, "nH      : %le    %le    %le\n", nH_avg, nH_avg_fion, nH_avg_binary);
    fprintf(file, "nHe     : %le    %le    %le\n", nHe_avg, nHe_avg_fion, nHe_avg_binary);
    fprintf(file, "ne      : %le       \n", ne_avg);
    fprintf(file, "GH1     : %le    %le    %le    %le    %le    %le    %le\n", gam_h1_vol, gam_h1_vol_fion, gam_h1_vol_binary, gam_h1_mass1, gam_h1_mass2, minarr(gamma_H1_tot), maxarr(gamma_H1_tot));
    fprintf(file, "RH2     : %le    %le    %le    %le    %le\n", rec_h2_vol, rec_h2_vol_fion, rec_h2_vol_binary, rec_h2_mass1, rec_h2_mass2);
	fprintf(file, "T       : %le    %le    %le    %le    %le\n", temp_vol , temp_vol_fion, temp_vol_binary, temp_mass, temp_mean_den);
	fprintf(file, "MFP     : %le    %le    %le    %le    %le    %le\n", mfp_vol,     mfp_vol_fion,     mfp_vol_binary,     mfp_mass1,  minarr(mfp), maxarr(mfp));
	fprintf(file, "MFP_912 : %le    %le    %le    %le\n", mfp_912_vol, mfp_912_vol_fion, mfp_912_vol_binary, mfp_912_mass1);
	fprintf(file, "Heat    : %le    %le    %le     %le    %le\n", heat_vol, heat_vol_fion, heat_vol_binary, heat_mass1, heat_mass2);
	fprintf(file, "Cool    : %le    %le    %le     %le    %le\n", cool_vol, cool_vol_fion, cool_vol_binary, cool_mass1, cool_mass2);
    fprintf(file, "npht    : %le       \n", photon_count);
	fprintf(file, "nsrc    : %d        \n", num_src);
	fprintf(file, "fion    : %le    %le\n", fion_avg_vol, fion_avg_mass);
	fprintf(file, "\n");

	fclose(file);
}
