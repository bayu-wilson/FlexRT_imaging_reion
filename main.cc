#include <ctime>
#include <cstdlib>
#include <omp.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "rt_funcs.cc"

int main(void)
{
	double start, end;
	double rt_start, rt_end;
	double move_start, move_end;
	double iter_start, iter_end;
	double split_start, split_end;
	double merge_start, merge_end;
	double rt_time_used = 0.; 
	double move_time_used = 0.;
	double iter_time_used = 0.;
	double split_time_used = 0.;
	double merge_time_used = 0.;
	double cpu_time_used;

printf("please work 1\n");
fflush(stdout);    

	//initialize random number generator for random angles
	init_rand();
	start = omp_get_wtime();

	//initialize the locks
	init_lock();
	
	//MODFIED 04/08/21 to initialize cosmic time function before initialization step on restart
	import_cosmic_time();
	make_scale_factor_table(); //ADDED 05/26/22 to enable effecient interpolation in both directions
printf("please work 2\n");
fflush(stdout);      	
	//initialize the grid IF restarting the simulation
	if (input_grid == TRUE)  {
		init();
	}
printf("please work 3\n");
fflush(stdout);      
	//update gas and source field
	read_hydro_steps();
	update_hydro_step();
	get_random_angles(); //generate random rotation angles
printf("please work 4\n");
fflush(stdout);	
	//initialize the grid IF NOT restarting the simulation
	if (input_grid == FALSE)  {
		init();
	}
	
	printf("Grid initialized.\n");
	printf("Gas initialized.\n");

	//grab the time step
	update_dt();
	
	printf("dt: %le\n", dt);
	printf("t_tot: %le\n", t_tot);
	
	if (multi_freq == TRUE)  {
		read_spectrum_table(); //read in the spectrum table, if doing multi_freq
	}
printf("please work 5\n");
fflush(stdout);
	set_spectrum(); //get the spectrum if doing multifrequency RT
	init_rt(); //initialize the RT grid
	
	set_ray_per_cell(); //initialize the maximum number of rays per cell
	set_lsource(); //set the healpix level for casting rays from the sources
printf("please work 6\n");
fflush(stdout);
	//MODIFIED 07/26/22 to cast an initial set of rays only if not reading from a file.  
	if (input_grid == FALSE)  {
		new_rays(); //cast the first set of rays
	}
	
	printf("RT data initialized.\n");

	//initialize sub-grid model data
	if (subgrid == TRUE)  {
		calc_mfp_boost(); //calculate the boost factor for the MFP (ADDED 03/14)
		import_delta_over_sigma_table(); //get the table for going from delta_lin/sigma <-> deltaNL
		init_radhydro_data(); //get the data for RadHydro
	}
printf("please work 7\n");
fflush(stdout);
	make_output(); //make output files
printf("please work 8\n");
fflush(stdout);
	if (input_grid == TRUE)  {
		write_otf();
	}

	end   = omp_get_wtime();
	cpu_time_used = ((double) (end - start));
	printf("%le\n",cpu_time_used);
	
	//Get the minimum and maximum photon count from the source distribution
	get_nphot_min_max();
	
	//while (t < 2.*dt)  {
	while (t < t_tot)  {
		step += 1;

		//source spectrum may evolve with time
		if (evolv_spec == TRUE)  {
			set_spectrum();
		}
		
		//set the speed of light
		set_clump(); //set the clumping factor
		set_clight(); //set the speed of light
		rt_start = omp_get_wtime();
	

		//update time
		t += dt;
		
		//update the conversion between delta/sigma <-> deltaNL
		if (subgrid == TRUE)  {
			delta_dom_to_deltaNL_dom(maxd(cosmic_time(1./(1.+12)), t + tinit));
		}
		
		//advance rays
		update_rays();
		
		rt_end = omp_get_wtime();
		rt_time_used += ((double) (rt_end - rt_start));
		
		iter_start = omp_get_wtime();

		//update energy
		if (c2ray_iter == TRUE)  {
			find_avg_unu(); //if the C2Ray is on, update energy here
		}
		else  {
			set_to_true(); //if not, set all rays to active
		}
		
		iter_end = omp_get_wtime();
		iter_time_used += ((double) (iter_end - iter_start));
		
		merge_start = omp_get_wtime();

		//merge rays
		
		set_lmerge(); //set the merging level
		set_merge_domain(); //initialize the domain for ray merging
		
		//get the normalized photon count for the merging threshold
		get_num_thresh(); //number of rays for threshold
		long int count = count_rays();
		if (count > num_thresh)  {
			calc_nphot_merge(); //calulate photon PDF and get merging thershold
		}
		

		//TEMPORARILY TURNED OFF MERGING!
		merge(); //merge rays
		
		merge_end = omp_get_wtime();
		merge_time_used += ((double) (merge_end - merge_start));
		
		split_start = omp_get_wtime();
		
		//split rays if the split rays would not be under the merging and/or deletion thresholds
		split_rays();
		
		set_ray_per_cell(); //set the maximum number of rays per cell
		set_lsource(); //source healpix level
		srcstep = mod(srcstep + num_src/Ncpu, num_src); //source step offset
		
		split_end = omp_get_wtime();
		split_time_used += ((double) (split_end - split_start));
		
		//update photo-ionization rate here if the c2ray scheme is off
		if (c2ray_iter == FALSE)  {
			update_gamma(); 
		}
		
		//update chemistry equations if C2Ray is off
		if (c2ray_iter == FALSE)  {
			update_chem();
		}
		
		//update heating/cooling rates and then temperature
		update_thermal();
		
		//update_dt(); //update time step
		
		//MODIFIED 07/26/22 to cast new rays prior to comsmo expansion/hydro updates
		new_rays(); //make new rays for the next time step
		
		//correct for cosmo expansion (ADDED 05/26/22)
		update_cosmo_expansion();
		
		//now update the hydro step if needed
		if (t >= t_tot)  {
			update_hydro_step(); //update the source and gas files
			get_nphot_min_max(); //set the minimum/maximum photon count
			//write_vIF_binary(); //try it here // commented by BAYU 11/10/22
		}

		//if ( (zz < 5.9) && (zz > 5.3) )  { //changed from 5.6 to 5.9 BAYU 11/10/22
		if ( (zz < 6.1) && (zz > 5.7) )  { 
			write_vIF_binary();
		}

		//MODIFIED 01/27/22 to fix fencepost error in timestep update.  
		update_dt(); //update time step
		
		move_start = omp_get_wtime();
		
		move_end = omp_get_wtime();
		move_time_used += ((double) (move_end - move_start));
		
		//if ( (mod(step, NumStep) == 0) || (t == t_tot) )  {
		//	printf("Writing on-the-fly output\n");
		//	write_otf();
		//	update_step();
		//}
		
		//MODIFIED 05/27/22 - get the total heating and Gamma for MP-Gadget
		if (MP_GADGET_HEAT == TRUE)  {
			update_gadget_var();
		}
		
		if ( (mod(step, NumStep) == 0) || (t == t_tot) )  {
			printf("Writing on-the-fly output\n");
			//FAHAD
			if ( (MP_GADGET_HEAT==TRUE) && (mod(step, gadget_step)==0) ){
				//Write file on every RT step
				write_gadget_binary(toTime(zinit)+t);
				reset_gadget_var();
			}
			write_otf();
			update_step();
		}
		
		//Added 08/31/22 - write vIF stuff. Will add condition to write statment later
		//write_vIF_binary(); //un-commented by BAYU, 2022/11/10 
		
		//update_ray_counts(); //update the cell-by-cell ray counts
		//get_cpu_ranks(); //sort the cpus by ray count
		balance_cpus(); //ADDED 05/30/22 - balance rays across cpus
		
		for (int nproc = 0; nproc < Ncpu; nproc++)  {
			ray_move_counter[nproc] = 0;
			long int proc_counter = 0;
			for (long int x = nproc*size_cpu; x < (nproc + 1)*size_cpu; x++)  {
				if (ray_tags[x] == TRUE)  {
					proc_counter += 1;
				}
			}
		}
		printf("step: %d\n", step);
		printf("time: %le\n", t/yr_to_s/1e6);
	}
	
	end   = omp_get_wtime();
	cpu_time_used = ((double) (end - start));
	printf("%le\n",t/yr_to_s/1e6);
	printf("Total time run: ");
	printf("%le\n",cpu_time_used);
	printf("RT time: ");
	printf("%le\n", rt_time_used);
	printf("ray creation time: ");
	printf("%le\n", move_time_used);
	printf("Iteration time: ");
	printf("%le\n", iter_time_used);
	printf("Merging time: ");
	printf("%le\n", merge_time_used);
	printf("Splitting time: ");
	printf("%le\n", split_time_used);
	printf("Total RT time: ");
	printf("%le\n", rt_time_used + move_time_used + iter_time_used + merge_time_used + split_time_used);
}
