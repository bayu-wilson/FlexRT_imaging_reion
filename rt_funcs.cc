#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "ray_tracing.cc"

//Update time step.  
void update_dt(void)  {
	printf("updating dt\n");
	dt = tstep_factor*dx[0]/cl_factor/cl;
}

//update the ionization state fractions used to calculate the i-front position and 
//speed (this is only helpful for single source tests).  
void update_step(void)  {
	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				f_H1_step[i][j][k]  = f_H1[i][j][k];
				f_H2_step[i][j][k]  = f_H2[i][j][k];
				f_He1_step[i][j][k] = f_He1[i][j][k];
				f_He2_step[i][j][k] = f_He2[i][j][k];
				f_He3_step[i][j][k] = f_He3[i][j][k];
			}
		}
	}
	}
	t_step = t;
}

//Initialize RT variables.  
void init_rt(void)  {
	
	//Initialze domain and buffer variables
	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				
				ray_merge[i][j][k].resize(num_directions, -1); //merging domain
				ray_merge_flag[i][j][k].resize(num_directions, FALSE);
				
				if (input_grid == FALSE)  {
					if (ifront == TRUE)  {
						fion[i][j][k]     = 0.;
						fion_max[i][j][k] = 0.; //Added 08/23/22
					}
					else  {
						fion[i][j][k]     = 1.;
						fion_max[i][j][k] = 1.;
					}
					tcross[i][j][k] = 0.;
					
					//Added 08/19/22
					for (int tt = 0; tt < Ntcross; tt++)  {
						tcross_new[i][j][k][tt]   = 0.;
						tcross_check[i][j][k][tt] = FALSE; //Aedded 08/23/22
					}
					
					//set RT flags
					treion_flag[i][j][k] = FALSE;
					tcross_flag[i][j][k] = FALSE;
					neutral_flag[i][j][k] = TRUE;
				}
				
				for (long int pix = 0; pix < size_cell; pix++)  {
					long int x = i*Ny*Nz*size_cell + j*Nz*size_cell + k*size_cell + pix;
					ray_end[x] = 0;
					ray_tags[x] = FALSE; //buffer tags
				}
			}
		}
	}
	}
	
	#pragma omp parallel
	{
	#pragma omp for 
	for (int src = 0; src < num_src; src++)  {
		if (source_lum[src].lum != 0)  {
			int i = source_lum[src].i;
			int j = source_lum[src].j;
			int k = source_lum[src].k;
			
			//MODIFIED 07/26/22 to work with new restart setup
			if (input_grid == FALSE)  {
				ray_dt[i][j][k] = dt; //time for source
			}
			else  {
				ray_dt[i][j][k] = 0.;
			}

		}
	}
	}
	
	//If reading from an input grid, read in the existing rays.  
	if (input_grid == TRUE)  {
		printf("Reading ray data\n");
		read_rays_binary(ray_file_init);
	}
	
	//if no input grid, set the first open spot to n = nproc
	if (input_grid == FALSE)  {
		//get the first free spot
		for (int n = 0; n < Ncpu; n++)   {
			ray_free[n] = (long int) n*size_cpu;
		}
	}
	//otherwise, locate a free index for each processor
	else  {
		for (int n = 0; n < Ncpu; n++)   {
			ray_free[n] = (long int) n*size_cpu;
			long int x = select_index(n); //Not sure if this is correct for restarts
			if (ray[x].nphot[0] == 0.)   {
				remove_ray(x);
			}
		}
	}
}

//set the speed of light and the number of cell crossings per ray update.  
void set_clight(void)  {
	
	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				clight[i][j][k] = cl_factor*cl;
				nc    [i][j][k] = nl_const;
			}
		}
	}
	}
}

//set active flag to true for all rays that exist
void set_to_true(void)  {
	#pragma omp parallel
	{
	#pragma omp for
	for (long int x = 0; x < size_buffer; x++)  {
		if ( (ray_tags[x] == TRUE) && (ray[x].active == FALSE) )  {
			ray[x].active = TRUE;
		}
	}
	}
}

//Loop over all rays and update their positions.  
void update_rays(void)  {
	int ray_counter = 0;
	double photon_counter = 0.;
	
	#pragma omp parallel
	{
	#pragma omp for 
	for (int src = 0; src < num_src; src++)  {
		int i = source_lum[src].i;
		int j = source_lum[src].j;
		int k = source_lum[src].k;
		ray_dt[i][j][k] += dt; //time for source
	}
	}
	
	if (c2ray_iter == FALSE)  {
		#pragma omp parallel 
		{
		#pragma omp for
		for (int i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  {
					//if not using the c2ray scheme, initialize the energy here
					for (int nu = 0; nu < Nfreq; nu++)  {
						u_nu[i][j][k][nu] = 0.;
					}
				}
			}
		}
		}
	}
	//advance all the rays by one time step
	#pragma omp parallel reduction (+ : ray_counter, photon_counter)
	{
	#pragma omp for
	for (long int x = 0; x < size_buffer; x++)  {
		if (ray_tags[x] == TRUE)  {
			advance_ray(x);
			ray_counter += 1; //update ray and photon counts
			photon_counter += ((double) sum_nphot(x))*phot_norm;
		}
	}
	}
	
	printf("total rays: \n");
	printf("%d\n", ray_counter);
	printf("%d\n", total_rays);
	printf("total photons propagated: ");
	printf("%le\n", photon_counter);
}

//Check each ray and split if needed
void split_rays(void)  {
	
	#pragma omp parallel
	{
	#pragma omp for 
	for (long int x = 0; x < size_buffer; x++)  {
		if ( (ray_tags[x] == TRUE) && (ray[x].active == TRUE) )  { //split only if a ray exists and is active (2nd condition avoids repeat splits)
			short int l = ray[x].order;
			float r 	= ray[x].r;
			//splitting condition: if r > r(l) and l is less than the maximum and if the resulting rays would not immediately be merged or deleted
			if ( (r > rmax_from_l(l)) && (l < l_max) && (sum_nphot(x) > maxd(4.*nphot_merge, 4.*nphot_min)) )  {
			//if ( (r > rmax_from_l(l)) && (l < l_max) )  { //TESTING!!
				int nproc = omp_get_thread_num();
				ray_move_counter[nproc] += 1;
				
				int ind = ray[x].path[0].ind; //grab index
				int ic = (ind - mod(ind, Ny*Nz))/(Ny*Nz);
				int jc = ind - ic*Ny*Nz;
				jc = (jc - mod(jc, Nz))/Nz;
				int kc  = ind - ic*Ny*Nz - jc*Nz;
				split_ray(ic, jc, kc, x);
			}
		}
	}
	}
	set_to_true(); //set all rays to active
}

//Create new rays at source cells
void new_rays(void)  {
	long int rays_created = 0;
	double lum_tot = 0.;
	
	#pragma omp parallel reduction (+ : rays_created, lum_tot)
	{
	#pragma omp for 
	for (int src = 0; src < num_src; src++)  {
		int i = source_lum[src].i;
		int j = source_lum[src].j;
		int k = source_lum[src].k;
		double lum = source_lum[src].lum;
		
		init_rays(i, j, k, lum);
		rays_created += 12*pow(2, 2*l_source[i][j][k]);
		lum_tot += lum;
		ray_dt[i][j][k] = 0.;
	}
	}
	total_rays += rays_created;
	printf("photon production rate: %le\n", lum_tot);
}
