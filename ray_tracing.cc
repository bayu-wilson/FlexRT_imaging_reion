#include <omp.h>
#include <stdio.h>
#include <iostream> 
#include <vector>
#include <array>
#include <algorithm>
#include <thread>
#include <pthread.h>
//#include </home/ccain002/Healpix_3.40/src/cxx/Healpix_cxx/healpix_base.h>
//#include </home/bwilson/scratch16-ansond/chris/Bayu_sims/Healpix_3.40/src/cxx/optimized_gcc/include/healpix_base.h>
#include </expanse/lustre/projects/uot171/bwils033/FlexRT/OTHER/Healpix_3.40/src/cxx/optimized_gcc/include/healpix_base.h>
#include "mem_funcs.cc" // /home/bwilson/scratch16-ansond/chris/Bayu_sims/Healpix_3.40/src/cxx/optimized_gcc/include
using namespace std;

//ray parameters
//order = healpix level
//pixel = healpix pixel number
//xhat = direction vector x component
//yhat = direction vector y component
//zhat = direction vector z component
//r    = radius (comoving units)
//dr   = distance to next split (comoving units)
//nphot  = normalized photon count
//x  = x position within a cell
//y  = y position
//z = z position
//time = time
//path = path through the box on a time step
//path.ind = cell index along path
//path.dist = distance traveled in cell ind

//MODIFIED 03/07 all Healpix_Base calls to explicitly initialize an order between l_min and l_max inclusive

//count the total number of rays in the box. 
long int count_rays(void)  {
	long int ray_count = 0;
	#pragma omp parallel reduction (+ : ray_count)
	{
	#pragma omp for
	for (long int x = 0; x < size_buffer; x++)  {
		if (ray_tags[x] == TRUE)  {
			ray_count += 1;
		}
	}
	}
	return ray_count;
}

//Set the maximum number of rays per cell.  This quantity can be spatially inhomogeneous, but we'll 
//keep it constant for now.  
void set_ray_per_cell(void)  {
	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < Nx; i ++)   {
		for (int j = 0; j < Ny; j++)   {
			for (int k = 0; k < Nz; k++)   {
				ray_per_cell[i][j][k] = ray_per_cell_max;
			}
		}
	}
	}
}

//Compute the minimum value of radius for order l
float rmin_from_l(short int l)  {
	//MODIFIED 05/31/22 to use comoving units for r
	float r = dx0*pow(2, l-1)*pow(3./pi/safety_factor, 0.5);
	if (l == 0)  { //honestly not sure if this should be here
		r = 0.; //TESTING leaving this part out, 02/24/21
	}
	return r;
}

//Maximum radius for order l
float rmax_from_l(short int l)  {
	//MODIFIED 05/31/22 to use comoving units for r
	return dx0*pow(2, l)*pow(3./pi/safety_factor, 0.5);
}

//Get the distance til the next split
float get_dr(short int l, float r)  {
	return rmax_from_l(l) - r;
}

//Get the healpix order from radius
short int l_from_r(float r)  {
	//MODIFIED 05/31/22 to use comoving units for r
	short int l = (short int) min((int) ceil(log(pi*safety_factor/3.*pow(r/dx0, 2.))/log(4.)), (int) l_max);
	l = mins(l_max, maxs(l_min, l)); //Added condition for l_max 02/25/21
	if (r/dx0 < 1.)  {
		l = l_min;
	}
	return l;
}

//Set the level for merging rays.  For now keep this as the casting level
//but eventually this can be made spatially/temporally variable
void set_lmerge(void)  {
	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < Nx; i++)   {
		for (int j = 0; j < Ny; j++)   {
			for (int k = 0; k < Nz; k++)   {
				l_merge[i][j][k] = hpx_lvl;
			}
		}
	}
	}
}

//Set the level for casting rays.  
void set_lsource(void)  {
	#pragma omp parallel
	{
	#pragma omp for
	for (int src = 0; src < num_src; src++)  {
		int i = source_lum[src].i;
		int j = source_lum[src].j;
		int k = source_lum[src].k;
		l_source[i][j][k] = hpx_lvl;
	}
	}
}

//Get the minimum and maximum photon count from the source distribution.  This is 1e-6, 1e+3 times
//the logarithmic average of the source base ray photon count.  
void get_nphot_min_max(void)  {
	int counter = 0;
	nphot_base_avg = 0.;
	#pragma omp parallel reduction (+ : nphot_base_avg, counter)
	{
	#pragma omp for
	for (int src = 0; src < num_src; src++)  {
		int i = source_lum[src].i;
		int j = source_lum[src].j;
		int k = source_lum[src].k;
		if (source_lum[src].lum > 0.)  {
			double nphot_temp = 0.;
			for (int nu = 0; nu < Nfreq; nu++)  {
				nphot_temp += spectrum[nu]*source_lum[src].lum*dt/12./pow(2.,2.*l_source[i][j][k])/phot_norm;
			}
			nphot_base_avg += (double) log10(nphot_temp);
			counter += 1;
		}
	}
	}
	nphot_base_avg /= counter;
	nphot_base_avg = pow(10., nphot_base_avg);
	
	nphot_max = 1e3*nphot_base_avg;
	nphot_min = 1e-10*nphot_base_avg;
}

//Get the number of rays above which merging starts
void get_num_thresh(void)  {
	num_thresh = 0;
	#pragma omp parallel reduction (+ : num_thresh)
	{
	#pragma omp for
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				num_thresh += (long int) (ray_per_cell[i][j][k] - 12*pow(2, 2*l_merge[i][j][k]));
			}
		}
	}
	}
}

//Compute the photon threshold for merging.  The photon count pdf of the rays is computed.  The first 
//(ray_per_cell_max - num_directions)*N^3 rays are excluded from merging, and a photon threshold is selected to 
//consider the remaining rays for merging.  This ensures that the total ray count will converge to ray_per_cell_max*N^3
//when the box is fully ionized.  
void calc_nphot_merge(void)  {
	
	double lg_max = log10(nphot_max);
	double lg_min = log10(nphot_min);
	double del = (lg_max - lg_min)/(10000 - 1);
	int nphot_pdf[10000];
	int ray_cum = 0;
	
	//initialize the pdf on a logarithmic grid
    #pragma omp parallel
	{
	#pragma omp for
	for (int n = 0; n < 10000; n++)  {
		nphot_pdf[n] = 0;
	}
    }
	
	//compute the pdf for all rays
	nphot_merge = 0.;
	#pragma omp parallel
	{
	#pragma omp for
	//loop over all rays in the buffer
	for (long int x = 0; x < size_buffer; x++)  {
		if (ray_tags[x] == TRUE)  { //check if ray exists
			int pdf_ind = (int) maxd((log10(sum_nphot(x)) - lg_min)/del, 0.);
			
			pdf_ind = min(pdf_ind, 9999);
			omp_set_lock(&(pdf_ind_lock[pdf_ind])); //use lock to avoid race conditions
			nphot_pdf[pdf_ind] += 1;
			omp_unset_lock(&(pdf_ind_lock[pdf_ind]));
		}
	}
	}
	
	//Start at the top of the pdf and iterate down until the cumulative photon count exceeds
	//the number of rays that should not be considered for merging.  
	int pdf_ind = 9999;
	while (pdf_ind >= 0) {
		ray_cum += nphot_pdf[pdf_ind];
		if (ray_cum >= num_thresh)  {
			nphot_merge = pow(10., pdf_ind*del + lg_min);
			printf("nphot_merge: %le\n", nphot_merge);
			pdf_ind = -1;
		}
		else  {
			pdf_ind -= 1;
		}
	}
	printf("number of rays counted: %d\n", ray_cum);
	printf("nphot_merge end: %le\n", nphot_merge);
}

//Initialize rays on a source cell at the end of a time step
void init_rays(int i, int j, int k, double lum)  {
	
	//get index for random angles
	unsigned short int itheta = siRand(Ntheta);
	//select random rotation angles
	float alpha0 = random_angles[itheta][0];
	float beta0  = random_angles[itheta][1];
	float gamma0 = random_angles[itheta][2];
	
	//will help when there is a very large number of source cells that are still opaque.  
	int Nside;
	short int lpix;
	if ( (subgrid == TRUE) && (lum*ray_dt[i][j][k] < (1. - fion[i][j][k])*nH[i][j][k]*pow(dx[0], 3.)) )  {
		Nside = 1;
		lpix = 0;
	}
	else if ( (subgrid == FALSE) && (nH1[i][j][k]/nH[i][j][k] > 0.1) )  {
		Nside = 1;
		lpix = 0;
	}
	else  {
		Nside = (int) pow(2, l_source[i][j][k]);
		lpix = l_source[i][j][k];
	}
	
	T_Healpix_Base<int> b((int) lpix, NEST);
	
	//loop over pixels
	for (int f = 0; f < 12; f++)  {
		for (int pnp = 0; pnp < pow(Nside, 2); pnp++)  {
			
			//loop over all pixel numbers on the Healpix sphere
			int pix = f*pow(Nside, 2) + pnp;
			//grab an index from the buffer allocated to nproc
			int nproc = omp_get_thread_num();
			omp_set_lock(&(proc_lock[nproc]));
			long int x = select_index(nproc);
			omp_unset_lock(&(proc_lock[nproc]));
			
			//if no index is available, move rays to another processor that is not full.  Pick the one with the fewest rays on it.  
			
			//03/08 undid parallelization change due to Healpix error
			//check other processors if the current one is full
			//04/08/21 put critical back inside the if statement after adding proc locks to the free_cpu subroutine
			
			if (x == -1)  {
				#pragma omp critical 
				{
				printf("processor full\n");
				update_ray_counts();
				get_cpu_ranks();
				free_cpu(nproc);
				
				omp_set_lock(&(proc_lock[nproc]));
				x = select_index(nproc);
				omp_unset_lock(&(proc_lock[nproc]));
				}
			}
			
			//grab the unit vector and rotate through random angles
			vec3 v = b.pix2vec(pix);
			v = get_unit_vector(v, alpha0, beta0, gamma0, 0);
			
			//assign ray to the buffer
			ray[x].active = TRUE;
			ray[x].order  = lpix;
			ray[x].pixel  = pix;
			ray[x].itheta = itheta;
			ray[x].xhat   = v.x;
			ray[x].yhat   = v.y;
			ray[x].zhat   = v.z;
			
			for (int nu = 0; nu < Nfreq; nu++)  {
				ray[x].nphot[nu] = (float) (spectrum[nu]*lum*ray_dt[i][j][k]/12./pow(Nside,2.)/phot_norm);
			}
			
			ray[x].r    = 0.; //rmin_from_l(lpix); FIXED 08/22/22
			ray[x].dr   = get_dr(lpix, ray[x].r);
			ray[x].x    = dx[i]/2.; //assign to the center of the source cell
			ray[x].y    = dy[j]/2.;
			ray[x].z    = dz[k]/2.;
			ray[x].time = 0.;
			
			ray[x].path[0].ind = i*Ny*Nz + j*Nz + k; 
			ray[x].path[0].dist = 0.;
			for (int a = 1; a < path_length; a++)  { //initialize path
				ray[x].path[a].ind = 0;
				ray[x].path[a].dist = 0.;
			}
		}
	}
}

//SUBROUTINE VERIFIED: 02/28 at 9:36 PM -> direct checks using external HealPix/rotation functions and direct math checks.  
void merge(void)  {
	
	//get index for random angles
	unsigned short int itheta = siRand(Ntheta);
	//select random rotation angles
	float alpha0 = random_angles[itheta][0];
	float beta0  = random_angles[itheta][1];
	float gamma0 = random_angles[itheta][2];
	
	#pragma omp parallel
	{
	#pragma omp for
	for (long int x = 0; x < size_buffer; x++)  {
	//loop over all rays in the buffer
		if ( (ray_tags[x] == TRUE) && (ray[x].active == TRUE) )  { //check if ray exists
			if (sum_nphot(x) < nphot_merge)  {
				int ind = ray[x].path[0].ind; //current location should always be stored at beginning of path
				//get cell indices
				int ic = (ind - mod(ind, Ny*Nz))/(Ny*Nz);
				int jc = ind - ic*Ny*Nz;
				jc = (jc - mod(jc, Nz))/Nz;
				int kc  = ind - ic*Ny*Nz - jc*Nz;
				
				//Healpix object.  
				T_Healpix_Base<int> b((int) l_merge[ic][jc][kc], NEST);
				
				vec3 v;
				v.x = ray[x].xhat; //get unit vector
				v.y = ray[x].yhat;
				v.z = ray[x].zhat;
				
				//rotate backwards through random angle and grab the pixel number in the rotated frame
				v = get_unit_vector(v, alpha0, beta0, gamma0, 6); 
				int pixel = (int) b.vec2pix(v); 
						
				//if there is no ray in the pixel bin, put the current one there.  
				//assign level/pixel number and photon-weighted quantities.  
				omp_set_lock(&(ray_lock[ic][jc][kc])); //lock to avoid race conditions on a cell
				if (ray_merge[ic][jc][kc][pixel] == -1)  {
					
					ray_merge[ic][jc][kc][pixel] 	  = x; //assign buffer index
					ray_merge_flag[ic][jc][kc][pixel] = FALSE; //should be false if only one ray in bin
		
					float nphot_total = sum_nphot(x);
		
					//MODIFIED 05/31/22 to use comoving units for r and dr
					ray[x].r    /= dx0;///= dx[ic]; //normalize to ensure no floating point error
					ray[x].r     = pow(ray[x].r, 2.)*nphot_total;
					ray[x].dr   *= nphot_total;
					ray[x].x    *= nphot_total;
					ray[x].y    *= nphot_total;
					ray[x].z    *= nphot_total;
					ray[x].time *= nphot_total;
				}
				//If there's already a base ray in the pixel bin, append the current ray to the base ray and delete it.  
				else  {
							
					long int x_merge = ray_merge[ic][jc][kc][pixel]; //grab index
					
					float nphot_total = sum_nphot(x);
					//append photon weighted quantities
					for (int nu = 0; nu < Nfreq; nu++)  {
						ray[x_merge].nphot[nu] += ray[x].nphot[nu];
					}
					
					//MODIFIED 05/31/22 to use comoving units for r
					ray[x].r          /= dx0;
					ray[x_merge].r    += pow(ray[x].r, 2.)*nphot_total;
					ray[x_merge].dr   += ray[x].dr*nphot_total;
					ray[x_merge].x    += ray[x].x*nphot_total;
					ray[x_merge].y    += ray[x].y*nphot_total;
					ray[x_merge].z    += ray[x].z*nphot_total;
					ray[x_merge].time += ray[x].time*nphot_total;
							
					ray_merge_flag[ic][jc][kc][pixel] = TRUE; //set to true to indicate merging should happen
					
					//remove merged ray
					remove_ray(x);
				}
				omp_unset_lock(&(ray_lock[ic][jc][kc]));
			}
		}
	}
	}
	
	//Loop over the merging domain and compute the properties of all outgoing rays
	#pragma omp parallel
	{
	#pragma omp for
	for (int ic = 0; ic < Nx; ic++)  {
		for (int jc = 0; jc < Ny; jc++)  {
			for (int kc = 0; kc < Nz; kc++)  {
				for (int pix = 0; pix < ray_merge[ic][jc][kc].size(); pix++)  {
					
					if (ray_merge[ic][jc][kc][pix] != -1)  { //check if at least one ray in buffer
						
						//Healpix object.  
						T_Healpix_Base<int> b((int) l_merge[ic][jc][kc], NEST);
						
						omp_set_lock(&(ray_lock[ic][jc][kc])); //avoid race conditions
						
						//if so, re-normalize everything and reset the merging bin
						long int x = ray_merge[ic][jc][kc][pix];
						
						//MODIFIED 05/31/22 to use comoving units for r and dr
						float nphot_total = sum_nphot(x);
						ray[x].r    /= nphot_total;
						ray[x].r     = pow(ray[x].r, 0.5);
						ray[x].r    *= dx0; //*= dx[ic];
						ray[x].dr   /= nphot_total;
						ray[x].x    /= nphot_total;
						ray[x].y    /= nphot_total;
						ray[x].z    /= nphot_total;
						ray[x].time /= nphot_total;
						
						//If multiple rays were combined, compute the new direction, order, radius, and pixel number of the outgoing ray
						if (ray_merge_flag[ic][jc][kc][pix] == TRUE)  {
							//normalize direction vector
							
							//grab the direction vector from the pixel number on l_merge sphere
							vec3 v      = b.pix2vec(pix);
							
							//update random angle index
							ray[x].itheta   = itheta; //update new angle index
							
							//MODIFIED on 03/07 to explicitly require l <= l_max
							//Get order from the merged radius, then recalculate r and dr
							ray[x].order = mins(l_max, l_from_r(ray[x].r));
							float rmin   = rmin_from_l(ray[x].order);
							float rmax   = rmax_from_l(ray[x].order);
							ray[x].r     = maxd(rmax - ray[x].dr, rmin);
							ray[x].dr    = rmax - ray[x].r;
							
							//Grab the pixel number from the outgoing vector
							T_Healpix_Base<int> b_merge((int) ray[x].order, NEST);
							ray[x].pixel = (int) b_merge.vec2pix(v);
							
							v = b_merge.pix2vec(ray[x].pixel);
							v = get_unit_vector(v, alpha0, beta0, gamma0, 0);
							ray[x].xhat 	= v.x;
							ray[x].yhat 	= v.y;
							ray[x].zhat 	= v.z;
							
							ray_merge_flag[ic][jc][kc][pix] = FALSE; //reset flag
						}
						ray_merge[ic][jc][kc][pix] = -1; //reset domain
						omp_unset_lock(&(ray_lock[ic][jc][kc]));
					}
				}
			}
		}
	}
	}
}

//Split ray at position x in the buffer in cell (ic, jc, kc).  The split ray is unpacked and deleted, 
//then the four child rays are created by increasing the order by 1. 
//VERIFIED SUBROUTINE: 02/28 at 5:48 PM 
void split_ray(int ic, int jc, int kc, long int x)  { 
	
	//unpack ray
	short int l  = ray[x].order;
	int pn 		 = ray[x].pixel;
	float r0 	 = ray[x].r;
	float x0_old = ray[x].x;
	float y0_old = ray[x].y;
	float z0_old = ray[x].z;
	float time 	 = ray[x].time;
	
	//unpack random angles
	unsigned short int itheta = ray[x].itheta;
	float alpha0 = random_angles[itheta][0];
	float beta0  = random_angles[itheta][1];
	float gamma0 = random_angles[itheta][2];
	
	//photons
	float nphot[Nfreq];
	for (int nu = 0; nu < Nfreq; nu++)  {
		nphot[nu] = ray[x].nphot[nu] / 4.;
	}
	
	//save the old direction vector
	vec3 v_old;
	v_old.x = ray[x].xhat;
	v_old.y = ray[x].yhat;
	v_old.z = ray[x].zhat;
	
	//delete the current ray from the buffer
	remove_ray(x);
	
	//initialize healpix objects for orders l, l+1
	T_Healpix_Base<int> b_old((int) l, NEST);
	T_Healpix_Base<int> b_new((int) l + 1, NEST);
	
	//Get the direction vector corresponding to the pixel of the parent ray
	//in the un-rotated frame.  
	vec3 v_old_pix = b_old.pix2vec(pn);
	
	//loop over child rays.  
	for (int n = 0; n < 4; n++)   {
		
		//new indices
		int ic_new = ic;
		int jc_new = jc;
		int kc_new = kc;
		
		int pix = 4*pn + n; //new pixel
		
		//direction vector of the new pixel in the un-rotated frame
		vec3 v_new_pix = b_new.pix2vec(pix);  
		
		//Rotate through the same rotation used on the parent ray
		vec3 v_new = get_unit_vector(v_new_pix, alpha0, beta0, gamma0, 0);
		
		//get displacement vector
		vec3 v_disp = v_new - v_old;
		
		//get new positions
		//MODIFIED 05/31/22 to account for the fact that r is now in comoving units
		float scale_factor = dx[0]/dx0;
		float x0 = x0_old + v_disp.x*(r0*scale_factor);
		float y0 = y0_old + v_disp.y*(r0*scale_factor);
		float z0 = z0_old + v_disp.z*(r0*scale_factor);
		
		//update positions and indices to find the new position and index of the child ray. 
		while (x0 > dx[ic_new])  {
			x0 -= dx[ic_new];
			ic_new = mod(ic_new + 1, Nx);
		}
		while (x0 < 0.)    {
			x0 += dx[ic_new];
			ic_new = mod(ic_new - 1, Nx);
		}
		while (y0 > dy[jc_new])  {
			y0 -= dx[jc_new];
			jc_new = mod(jc_new + 1, Ny);
		}
		while (y0 < 0.)    {
			y0 += dy[jc_new];
			jc_new = mod(jc_new - 1, Ny);
		}
		while (z0 > dz[kc_new])  {
			z0 -= dz[kc_new];
			kc_new = mod(kc_new + 1, Nz);
		}
		while (z0 < 0.)    {
			z0 += dz[kc_new];
			kc_new = mod(kc_new - 1, Nz);
		}
		
		//grab a free index from the buffer
		int nproc = omp_get_thread_num();
		omp_set_lock(&(proc_lock[nproc]));
		long int x_new = select_index(nproc);
		omp_unset_lock(&(proc_lock[nproc]));
		
		//MODIFIED 04/10/21 to put critical statement back inside if statement
		//check other processors if the current one is full
		if (x_new == -1)  {
			#pragma omp critical 
			{
			printf("processor full\n");
			update_ray_counts();
			get_cpu_ranks();
			free_cpu(nproc);
			omp_set_lock(&(proc_lock[nproc]));
			x_new = select_index(nproc);
			omp_unset_lock(&(proc_lock[nproc]));
			}
		}
		ray[x_new].active = FALSE;
		
		//Re-compute the radius
		float rmax = rmax_from_l(l + 1);
		
		//Assign child ray to buffer
		ray[x_new].order  = l + 1;
		ray[x_new].pixel  = pix;
		ray[x_new].itheta = itheta;
		ray[x_new].xhat   = v_new.x;
		ray[x_new].yhat   = v_new.y;
		ray[x_new].zhat   = v_new.z;
		
		for (int nu = 0; nu < Nfreq; nu++)  {
			ray[x_new].nphot[nu] = nphot[nu];
		}
		
		ray[x_new].r  = r0;
		ray[x_new].dr = rmax - r0;
		ray[x_new].x  = x0;
		ray[x_new].y  = y0;
		ray[x_new].z  = z0;
		ray[x_new].time = time;
		ray[x_new].path[0].ind = ic_new*Ny*Nz + jc_new*Nz + kc_new; //new cell location
	}
}

//Loop over all rays that moved and compute the radiation energy deposited in each cell.  This section uses
//the C2Ray procedure (if specified) to iteratively compute the photo-ionization rate and neutral fraction.  
void find_avg_unu(void)  {
	
	for (int iter = 0; iter < itercount; iter++)  { //iterate
	
		//Loop over cells.  Set gamma and energy to 0 and if iter = 0, set an initial guess for the neutral fraction.  
		//initialize the recombination rates on the first iteration, and if sub-grid = TRUE interpolate to get the initial clumping
		//factor and mfp
		#pragma omp parallel
		{
		#pragma omp for
		for (int i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  {
					avg_gamma[i][j][k] = 0.;
					avg_gamma_he1[i][j][k]  = 0.;
					ngam[i][j][k][0]              = 0.;
					ngam[i][j][k][1]              = 0.;
					ngam[i][j][k][2]              = 0.;
                                        
					//Initialize vIF variables
                                        IF_speed[i][j][k] = 0.;
                                        Inc_Flux[i][j][k] = 0.;
                                        Delta_fion[i][j][k] = 0.;
					
					for (int nu = 0; nu < Nfreq; nu++)  {
						ngam_nu[i][j][k][nu] = 0.;
						avg_gamma_nu[i][j][k][nu] 		= 0.;
						avg_gamma_nu_he1[i][j][k][nu] 	= 0.;
					}
					for (int nu = 0; nu < Nfreq; nu++)  {
						u_nu[i][j][k][nu] = 0.;
					}
					if (iter == 0)  {
						avg_nh1[i][j][k] = nH1[i][j][k]; //get initial guess for HI for neutral fraction
						avg_nhe1[i][j][k] = nHe1[i][j][k];
						
						//MODFIED 03/04/21 to correct for really low (or maybe negative?) temperatures
						//recombination coefficients
						if (strcmp(rec_case, "A") == 0)  {
							recomb_H2[i][j][k]  = (float) alphaA_H2(maxd(temp_min, temp[i][j][k]));
							recomb_He2[i][j][k] = (float) alphaA_He2(maxd(temp_min, temp[i][j][k]));
							recomb_He3[i][j][k] = (float) alphaA_He3(maxd(temp_min, temp[i][j][k]));
						}
						else  {
							recomb_H2[i][j][k]  = (float) alphaB_H2(maxd(temp_min, temp[i][j][k]));
							recomb_He2[i][j][k] = (float) alphaB_He2(maxd(temp_min, temp[i][j][k]));
							recomb_He3[i][j][k] = (float) alphaB_He3(maxd(temp_min, temp[i][j][k]));
						}
						//diaelectric recombination of He II.  
						if ( (temp[i][j][k] >= 3e4) && (temp[i][j][k] <= 1e6) )  {
							recomb_He2[i][j][k] += (float) Dalpha_He2(maxd(temp_min, temp[i][j][k]));
						}
						
						//time derivative of the electron density
						if (subgrid == FALSE)  {
							dne_dt[i][j][k]   = -ne[i][j][k];
						}
						else  { //assume the gas is fully ionized, so the only contribution to the derivative will be changes in density
							dne_dt[i][j][k] = -(nH[i][j][k] + nHe[i][j][k]);
						}
						
						//Initialize the sub-grid physics.  If sub-grid is true and the cell has already been ionized, initialize
						//interpolation functions and grab the MFP as a function of the gamma from the previous iteration.  
						//MODIFIED 03/05 to remove unnecessary MFP update step
						if ( (subgrid == TRUE) && (treion_flag[i][j][k] == TRUE) )  { //get mfp and xHI at as function of t, treion, and delta. 
							get_gamma_interp(i, j, k, treion[i][j][k] + tinit, rho[i][j][k]/avg_rho - 1., t - treion[i][j][k]);
						}
						else  { //otherwise assume uniform density and guess the MFP at the lowest frequency
							mfp[i][j][k] = 1./nH1[i][j][k]/sigmapi_H1(nu_phot[0]); 
						}
					}
				}
			}
		}
		}
		
		printf("iteration: %d\n", iter);
		
		#pragma omp parallel
		{
		#pragma omp for
		for (long int x = 0; x < size_buffer; x++)  {
			if ( (ray_tags[x] == TRUE) && (ray[x].active == FALSE) )   {
				//initialize opacity, photon count, and path index
				double nphot[Nfreq];
				double tau_h1[Nfreq], dtau_h1[Nfreq], tau_he1[Nfreq], dtau_he1[Nfreq];
				for (int nu = 0; nu < Nfreq; nu++)  {
					nphot[nu]    = ray[x].nphot[nu]*phot_norm;
					tau_h1[nu]   = 0.; //initialize opacity and path index
					dtau_h1[nu]  = 0.;
					tau_he1[nu]  = 0.;
					dtau_he1[nu] = 0.;
				}
				
				short int p = 0;
				short int p_prev = 0;
				float delta_s_eff;
				
				int nproc = omp_get_thread_num();
				
				while (p < ray_end[x])  { //loop over path
					int ind  = ray[x].path[p].ind;
					int ip = (ind - mod(ind, Ny*Nz))/(Ny*Nz);
					int jp = ind - ip*Ny*Nz;
					jp     = (jp - mod(jp, Nz))/Nz;
					int kp = ind - ip*Ny*Nz - jp*Nz;
					
					//If sub-grid is off, use the average h1 fraction to compute the change in 
					//opacity across the segment.  
					//Fixed helium opacity (for fractional absorption calculation) on 09/24/21 (in he2_tr_h2 == TRUE case)
					if (subgrid == FALSE)  {
						for (int nu = 0; nu < Nfreq; nu++)  {
							dtau_h1[nu]  = avg_nh1[ip][jp][kp]*sigmapi_H1(nu_phot[nu])*ray[x].path[p].dist;
							if (he2_tr_h2 == TRUE)  {
								dtau_he1[nu] = nhe_to_nh*dtau_h1[nu];
							}
					else  {
								dtau_he1[nu] = avg_nhe1[ip][jp][kp]*sigmapi_He1(nu_phot[nu])*ray[x].path[p].dist;
							}
						}
					}
					
					//initialize evolving photon count and cell volume
					double nphot_in[Nfreq];
					double vol = ((double) dx[ip])*((double) dy[jp])*((double) dz[kp]);
					
					//using the total tau from previous iterations, tally the number of photons left 
					//in the ray at it's current point along the path.  Also tally, in normalized units, 
					//the number of photons over all frequency bins.  
					double nphot_in_total = 0.;
					for (int nu = 0; nu < Nfreq; nu++)  {
						nphot_in[nu]     = nphot[nu]*exp(-tau_h1[nu] - tau_he1[nu]);
						nphot_in_total += nphot_in[nu]/phot_norm;
					}
					
					//initialize fraction of radiation absorbed by HI and HeI.  
					double frac_h1[Nfreq], frac_he1[Nfreq];
					int nproc = omp_get_thread_num();
					
					//get absorption fractions for each species.  
					for (int nu = 0; nu < Nfreq; nu++)  {
						if ( (dtau_h1[nu] > 0.) || (dtau_he1[nu] > 0.) || (subgrid == TRUE) )  {
							if ( (subgrid == TRUE) || (he2_tr_h2 == TRUE) )  {
								frac_h1[nu]  = 1. - (Y/4.)/(1. - Y + Y/4.);
								frac_he1[nu] = (Y/4.)/(1. - Y + Y/4.);
							}
							else  {
								frac_h1[nu]  = (1. - exp(-dtau_h1[nu]))*exp(-dtau_he1[nu])/((1.- exp(-dtau_h1[nu]))*exp(-dtau_he1[nu]) + (1.- exp(-dtau_he1[nu]))*exp(-dtau_h1[nu]));
								frac_he1[nu] = (1. - exp(-dtau_he1[nu]))*exp(-dtau_h1[nu])/((1.- exp(-dtau_h1[nu]))*exp(-dtau_he1[nu]) + (1.- exp(-dtau_he1[nu]))*exp(-dtau_h1[nu]));
							}
							
						}
						else  {
							frac_h1[nu]  = 0.;
							frac_he1[nu] = 0.;
						}
					}
					
					//If the number of photons in the ray is currently below the minimum, dump the entire ray in the cell
					if (nphot_in_total < nphot_min)  { 
						omp_set_lock(&(gam_lock[ip][jp][kp]));
						for (int nu = 0; nu < Nfreq; nu++)  {
							ngam_nu[ip][jp][kp][nu]          += (float) (nphot_in[nu]/phot_norm); //keep track of spectrum
							if (subgrid == TRUE)  {
								if (fion[ip][jp][kp] > 0.)  {
									if (iter == 0)  {
										delta_s_eff = fion[ip][jp][kp]*ray[x].path[p].dist/frac_h1[nu]; //set initial guess on the high end
									}
									else  {
										delta_s_eff = mfp[ip][jp][kp]*(1. - exp(-fion[ip][jp][kp]*ray[x].path[p].dist/mfp[ip][jp][kp]/frac_h1[nu])); //account for low mfp absorbing rays before they can escape
									}
									//omp_set_lock(&(gam_lock[ip][jp][kp]));
									avg_gamma_nu[ip][jp][kp][nu]	 += (float) (frac_h1[nu]*nphot_in[nu]*sigma_eff_rh*delta_s_eff/dt/vol/fion[ip][jp][kp]);
									avg_gamma[ip][jp][kp] 			 += (float) (frac_h1[nu]*nphot_in[nu]*sigma_eff_rh*delta_s_eff/dt/vol/fion[ip][jp][kp]);
									avg_gamma_nu_he1[ip][jp][kp][nu] += (float) (frac_he1[nu]*nphot_in[nu]*sigma_eff_rh*delta_s_eff/dt/vol/fion[ip][jp][kp]);
									avg_gamma_he1[ip][jp][kp] 		 += (float) (frac_he1[nu]*nphot_in[nu]*sigma_eff_rh*delta_s_eff/dt/vol/fion[ip][jp][kp]);
									//omp_unset_lock(&(gam_lock[ip][jp][kp]));
								}
							}
							else  {
								//omp_set_lock(&(gam_lock[ip][jp][kp]));
								avg_gamma_nu[ip][jp][kp][nu]	 += (float) (frac_h1[nu]*nphot_in[nu]/dt/vol/avg_nh1[ip][jp][kp]);
								avg_gamma[ip][jp][kp] 			 += (float) (frac_h1[nu]*nphot_in[nu]/dt/vol/avg_nh1[ip][jp][kp]);
								avg_gamma_nu_he1[ip][jp][kp][nu] += (float) (frac_he1[nu]*nphot_in[nu]/dt/vol/avg_nhe1[ip][jp][kp]);
								avg_gamma_he1[ip][jp][kp] 		 += (float) (frac_he1[nu]*nphot_in[nu]/dt/vol/avg_nhe1[ip][jp][kp]);
								//omp_unset_lock(&(gam_lock[ip][jp][kp]));
							}
						}
						omp_unset_lock(&(gam_lock[ip][jp][kp]));
						p = 99;
					}
					else  { //otherwise compute the photo-ionization rate from the number of photons deposited
						omp_set_lock(&(gam_lock[ip][jp][kp]));
						for (int nu = 0; nu < Nfreq; nu++)  {
							ngam_nu[ip][jp][kp][nu]          += (float) (nphot_in[nu]/phot_norm); //keep track of spectrum
							//if the sub-grid model is on, get gamma in the optically thin limit
							if (subgrid == TRUE)  {
								if (fion[ip][jp][kp] > 0.)  {
									if (iter == 0)  {
										delta_s_eff = fion[ip][jp][kp]*ray[x].path[p].dist/frac_h1[nu]; //set initial guess on the high end
									}
									else  {
										delta_s_eff = mfp[ip][jp][kp]*(1. - exp(-fion[ip][jp][kp]*ray[x].path[p].dist/mfp[ip][jp][kp]/frac_h1[nu])); //account for low mfp absorbing rays before they can escape
									}
									//omp_set_lock(&(gam_lock[ip][jp][kp]));
									avg_gamma_nu[ip][jp][kp][nu]	 += (float) (frac_h1[nu]*nphot_in[nu]*sigma_eff_rh*delta_s_eff/dt/vol/fion[ip][jp][kp]);
									avg_gamma[ip][jp][kp] 			 += (float) (frac_h1[nu]*nphot_in[nu]*sigma_eff_rh*delta_s_eff/dt/vol/fion[ip][jp][kp]);
									avg_gamma_nu_he1[ip][jp][kp][nu] += (float) (frac_he1[nu]*nphot_in[nu]*sigma_eff_rh*delta_s_eff/dt/vol/fion[ip][jp][kp]);
									avg_gamma_he1[ip][jp][kp] 		 += (float) (frac_he1[nu]*nphot_in[nu]*sigma_eff_rh*delta_s_eff/dt/vol/fion[ip][jp][kp]);
									//omp_unset_lock(&(gam_lock[ip][jp][kp]));
								}
							}
							else  {
								//omp_set_lock(&(gam_lock[ip][jp][kp]));
								avg_gamma_nu[ip][jp][kp][nu] 	 += (float) (frac_h1[nu]*nphot_in[nu]*(1. - exp(-dtau_h1[nu] - dtau_he1[nu]))/dt/vol/avg_nh1[ip][jp][kp]);
								avg_gamma[ip][jp][kp] 		  	 += (float) (frac_h1[nu]*nphot_in[nu]*(1. - exp(-dtau_h1[nu] - dtau_he1[nu]))/dt/vol/avg_nh1[ip][jp][kp]);
								avg_gamma_nu_he1[ip][jp][kp][nu] += (float) (frac_he1[nu]*nphot_in[nu]*(1. - exp(-dtau_h1[nu] - dtau_he1[nu]))/dt/vol/avg_nhe1[ip][jp][kp]);
								avg_gamma_he1[ip][jp][kp] 	  	 += (float) (frac_he1[nu]*nphot_in[nu]*(1. - exp(-dtau_h1[nu] - dtau_he1[nu]))/dt/vol/avg_nhe1[ip][jp][kp]);
								//omp_unset_lock(&(gam_lock[ip][jp][kp]));
							}
						}
						omp_unset_lock(&(gam_lock[ip][jp][kp]));
						p += 1;
					}
					
					//use the mfp from the sub-grid model if the cell is at least partially ionized
					if ( (subgrid == TRUE) && (fion[ip][jp][kp] > 0.) )  {
						for (int nu = 0; nu < Nfreq; nu++)  {
							dtau_h1[nu]  = fion[ip][jp][kp]*ray[x].path[p_prev].dist/mfp[ip][jp][kp];
							dtau_he1[nu] = nhe_to_nh*dtau_h1[nu];
							//if (source_tau[ip][jp][kp] > 0.)  { //add opacity contributing from the sources
							//	dtau_h1[nu] += source_tau[ip][jp][kp];
							//}
						}
					}
					
					//update the total opacity
					for (int nu = 0; nu < Nfreq; nu++)  {
						tau_h1[nu]  += dtau_h1[nu];
						tau_he1[nu] += dtau_he1[nu];
					}
					
					//If the ionized fraction is below 1, count how many photons will contribute to advancing the i-front
					//Seems as if there is a BUG here because removing the iter condition caused NANs in the output file
					if (fion[ip][jp][kp] < 1.)  {
						omp_set_lock(&(ray_lock[ip][jp][kp])); //lock to avoid race conditions on a cell
						double sum_ngam = ngam[ip][jp][kp][0]; //sum_over_procsd(ngam_proc[ip][jp][kp][0]); //number of photons absorbed so far by neutral part
						
						double chi = nhe_to_nh;
						//neglect nHI_res here.  Either front will pass quickly and nHI_res will be small, or it will pass slowly and the error
						//will contribute to only a small fraction of the total number of absorptions.  
						double free_ngam = vol*nH[ip][jp][kp]*(1. + chi)*(1. - fion[ip][jp][kp]); //HI atoms available
						double nphot_left = 0.; //photons remaining in the ray after absorption by the neutral part
					
						for (int nu = 0; nu < Nfreq; nu++)  {	
							//ngam[ip][jp][kp][2] += nphot_in[nu];
							ngam[ip][jp][kp][2] += nphot_in[nu]*exp(-(dtau_h1[nu] + dtau_he1[nu])); //Chris Feb. 29, 2024
						}

						//MODIFIED 04/10/21 to fix optically thin photo-ionization rate in newly ionized cells
						if (sum_ngam < free_ngam)  { //check if any more photons can be absorbed by the neutral part of the cell on this time step
							for (int nu = 0; nu < Nfreq; nu++)  {
								//double ngam_abs = mindd(maxd(0., nphot_in[nu] - recombs), free_ngam - sum_ngam);
								double ngam_abs = mindd(nphot_in[nu]*exp(-(dtau_h1[nu] + dtau_he1[nu])), free_ngam - sum_ngam); //number of photons absorbed into the NEUTRAL part
								//double ngam_abs = mindd(nphot_in[nu], free_ngam - sum_ngam); //ignore attenuation behind the front (because we can't estimate it accurately)
								ngam[ip][jp][kp][0] += ngam_abs; //count photons
								if (subgrid == TRUE)  {
									ngam[ip][jp][kp][1] += ngam_abs*sigma_eff_rh*ray[x].path[p_prev].dist;
								}
								else  {
									ngam[ip][jp][kp][1] += ngam_abs*sigmapi_H1(nu_phot[nu])*ray[x].path[p_prev].dist;
								}
								//count how many photons are left
								//nphot_in[nu] = maxd(0., nphot_in[nu] - recombs) - ngam_abs;
								nphot_in[nu] = nphot_in[nu]*exp(-(dtau_h1[nu] + dtau_he1[nu])) - ngam_abs;
								//nphot_in[nu] = nphot_in[nu] - ngam_abs;
								nphot_left += nphot_in[nu];
							}
							
							//If none are left, absorb the photons and delete the ray
							if (nphot_left == 0.)  {
								p = 99;
								for (int nu = 0; nu < Nfreq; nu++)  {
									tau_h1[nu] = 999.;
								}
							}
							else  { //otherwise, spread the effective opacity over the spectrum
								//renormalize spectrum
								//MODIFIED 05/24/22 to correct under-counting of absorptions for rays that finish ionizing a cell for the first time
								for (int nu = 0; nu < Nfreq; nu++)  { //need to subtract off dtau here because we are neglecting absorptions behind the i-front
									//tau_h1[nu] += -log(nphot_in[nu]/nphot[nu]/exp(-(tau_h1[nu] - dtau_h1[nu] + tau_he1[nu] - dtau_he1[nu])));
									tau_h1[nu] += -log(nphot_in[nu]/nphot[nu]/exp(-(tau_h1[nu] + tau_he1[nu])));
								}
							}
						}
						omp_unset_lock(&(ray_lock[ip][jp][kp])); //lock to avoid race conditions on a cell
					}
					p_prev = p;
				}
				//At the end of the last iteration, update the photon count
				if (iter == itercount - 1)  {
					for (int nu = 0; nu < Nfreq; nu++)  {
						ray[x].nphot[nu] = (float) (nphot[nu]*exp(-tau_h1[nu] - tau_he1[nu])/phot_norm);
					}
					if (sum_nphot(x) > nphot_min)  { //if the ray was not dumped, recompute the photon count and set to active
						ray[x].active = TRUE;
						ray[x].path[0].ind = ray[x].path[ray_end[x] - 1].ind; //re-label 1st path index
					}
					else  { //if below photon minimum, remove the ray
						remove_ray(x);
					}
					
				}
			}
		}
		}
		
		//Update the grid quantities.  This includes gamma, the ionization states of the gas, etc. 
		#pragma omp parallel
		{
		#pragma omp for
		for (int i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  { //loop over grid cells
					int nproc = omp_get_thread_num();
					
					//Fixed helium opacity (for fractional absorption calculation) on 09/(24-25)/21 (in he2_tr_h2 == TRUE case)
					//if (he2_tr_h2 == TRUE)  {
					//	avg_gamma[i][j][k] /= 1. + nhe_to_nh;
					//}
					
					//electron number density and cell volume
					double nelec = (double) (nH[i][j][k] - avg_nh1[i][j][k] + nHe[i][j][k] - avg_nhe1[i][j][k]);
					double vol = ((double) dx[i])*((double) dy[j])*((double) dz[k]);
					
					//Sub-grid model case.   
					if (subgrid == TRUE)  {
						
						double gamma_fion;
						bool first = FALSE;
						
						//if this is the first time photons have been absorbed by a cell, we will need to initialize the sub-grid
						//data and tag the cell as "ionized".  
						if ( (iter == itercount - 1) && (treion_flag[i][j][k] == FALSE) && ((ngam[i][j][k][1] > 0.) || (avg_gamma[i][j][k] > 0.)) ) {
							treion[i][j][k] = t - dt;
							treion_flag[i][j][k] = TRUE;
							get_gamma_interp(i, j, k, treion[i][j][k] + tinit, rho[i][j][k]/avg_rho - 1. , t - treion[i][j][k]);
							first = TRUE;
						}
						
						//Otherise, interpolate the MFP for the next iteration
						if ( (treion_flag[i][j][k] == TRUE) && (fion[i][j][k] > 0.) )  {
							
							if (equilibrium == FALSE)  {
								interp_mfp(i, j, k, avg_gamma_prev[i][j][k], avg_gamma[i][j][k]); //get final mfp/number densities
							}

							//MODIFIED 02/25/22 to adjust the MFP for density-dependent clumping factor and temperature variations
							if (equilibrium == TRUE)  {
									mfp[i][j][k] = maxd(gamma_min,avg_gamma[i][j][k])/clump[i][j][k]/sigma_eff_rh/alphaB_H2(maxd(temp_min,temp[i][j][k]))/pow(nH[i][j][k],2.)/1.082;
									mfp_912[i][j][k] = maxd(gamma_min,avg_gamma[i][j][k])/clump[i][j][k]/sigmapi_H1(13.6/h_eV)/alphaB_H2(maxd(temp_min,temp[i][j][k]))/pow(nH[i][j][k],2.)/1.082;
									mfp[i][j][k] = maxd(mfp_min, mfp[i][j][k]);
									mfp_912[i][j][k] = maxd(mfp_min, mfp_912[i][j][k]);
							}

							//interp_mfp(i, j, k, avg_gamma_prev[i][j][k], avg_gamma[i][j][k]); //get final mfp/number densities
							
							//MODIFIED 02/25/22 to adjust the MFP for density-dependent clumping factor and temperature variations
							//if (equilibrium == TRUE)  {
							//	mfp[i][j][k]     /= clump[i][j][k]; // const_clump;
							//	mfp_912[i][j][k] /= clump[i][j][k]; // const_clump;
							//	if (temp_ev == TRUE)  {
							//		mfp[i][j][k]     /= alphaB_H2(maxd(temp_min,temp[i][j][k]))/alphaB_H2(1e4);
							//		mfp_912[i][j][k] /= alphaB_H2(maxd(temp_min,temp[i][j][k]))/alphaB_H2(1e4);
							//	}
							//}
						}
						
						//MODIFIED on 03/02 to fix the t_cross calculation - removed ngam > 0 flag and replaced with treion == TRUE
						//I-front correction
						if ( (ifront == TRUE) && (treion_flag[i][j][k] == TRUE) && (iter == itercount - 1) && (fion[i][j][k] < 1.) )  {
							
							gamma_fion = ngam[i][j][k][1]/dt/vol;
							
							//update the ionization front position within the cell
							update_fion(i, j, k);
							if (fion[i][j][k] >= 1.)  {
								neutral_flag[i][j][k] = FALSE;
							}
							
							//MODIFIED 02/22/22 to avoid the "empty cell" recombination correction on the very first time step.
							//MODIFIED 07/28/22 to use cell-wise clumping factor for equilibrium runs
							float clump_fion = 0.;
							if (equilibrium == FALSE)  {
									clump_fion = 3.; //this is CR
									if ( (fion[i][j][k] < 1.) && (avg_gamma[i][j][k] == 0.) && (first == FALSE) )  {
											fion[i][j][k] = maxd(0., fion[i][j][k] - clump_fion*alphaB_H2(1e4)*(1.+nhe_to_nh)*nH[i][j][k]*fion[i][j][k]*dt);
									}
							}
							else  {
									clump_fion = clump[i][j][k]; // this is C_HII
									if ( (fion[i][j][k] < 1.) && (avg_gamma[i][j][k] == 0.) && (first == FALSE) )  {
											fion[i][j][k] = maxd(0., fion[i][j][k] - clump_fion*alphaB_H2(maxd(temp_min,temp[i][j][k]))*(1.+nhe_to_nh)*nH[i][j][k]*fion[i][j][k]*dt);
									}
							}

							
							//if a cell is being ionized for the first time, the MFP in the ionized part should be initialized tusing the incident flux on the cells
							//MODIFIED 03/06 to ensure that cells that are ionized in one time step do not end up with pathological initial conditions.   
							if (first == TRUE)  {
								mfp_prev[i][j][k] = 0.; //initial condition for mfp
								avg_gamma_prev[i][j][k] = 0.; //initial condition for gamma
								avg_gamma[i][j][k] = gamma_fion;
								//interp_mfp(i, j, k, avg_gamma_prev[i][j][k], avg_gamma[i][j][k]);
								if (equilibrium == FALSE)  {
										interp_mfp(i, j, k, avg_gamma_prev[i][j][k], avg_gamma[i][j][k]);
								}
								if (equilibrium == TRUE)  {
										mfp[i][j][k] = maxd(gamma_min,avg_gamma[i][j][k])/clump[i][j][k]/sigma_eff_rh/alphaB_H2(maxd(temp_min,temp[i][j][k]))/pow(nH[i][j][k],2.)/1.082;
										mfp_912[i][j][k] = maxd(gamma_min,avg_gamma[i][j][k])/clump[i][j][k]/sigmapi_H1(13.6/h_eV)/alphaB_H2(maxd(temp_min,temp[i][j][k]))/pow(nH[i][j][k],2.)/1.082;
										mfp[i][j][k] = maxd(mfp_min, mfp[i][j][k]);
										mfp_912[i][j][k] = maxd(mfp_min, mfp_912[i][j][k]);
								}

							}
						}
						
						//ADDED 07/28/22 to account for fully ionized empty cells
						if ( (equilibrium == TRUE) && (fion[i][j][k] == 1.) && (avg_gamma[i][j][k] == 0.) && (iter == itercount - 1) && (ifront == TRUE) )  {
							float clump_fion_2 = clump[i][j][k];
							fion[i][j][k] = maxd(0., fion[i][j][k] - clump_fion_2*alphaB_H2(maxd(temp_min,temp[i][j][k]))*(1.+nhe_to_nh)*nH[i][j][k]*fion[i][j][k]*dt);
						}

						if ( (equilibrium == FALSE) && (fion[i][j][k] == 1.) && (avg_gamma[i][j][k] == 0.) && (iter == itercount - 1) && (ifront == TRUE) )  {
							float clump_fion_2 = 3.;
							fion[i][j][k] = maxd(0., fion[i][j][k] - clump_fion_2*alphaB_H2(maxd(temp_min,temp[i][j][k]))*(1.+nhe_to_nh)*nH[i][j][k]*fion[i][j][k]*dt);
						}

						
						//update relevant quantities on the last iteration.  
						if ( (iter == itercount - 1) && (treion_flag[i][j][k] == TRUE) )  {
							interp_nh1(i, j, k, avg_gamma[i][j][k]);
							interp_nhe1(i, j, k, avg_gamma[i][j][k]);
							nH2[i][j][k] = nH[i][j][k]*fion[i][j][k]; // - nH1[i][j][k];
							nHe2[i][j][k] = nHe[i][j][k]*fion[i][j][k]; // - nHe1[i][j][k];
							nHe3[i][j][k] = 0.;
							//update previous quantities
							avg_gamma_prev[i][j][k] = avg_gamma[i][j][k];
							mfp_prev[i][j][k] = mfp[i][j][k];
							mfp_912_prev[i][j][k] = mfp_912[i][j][k];
							
						}
					}
					//If not using the sub-grid model, proceed with the C2Ray scheme.  
					else  {
						mfp[i][j][k]        = 1./sigmapi_H1(nu_phot[0])/avg_nh1[i][j][k];
						mfp_prev[i][j][k]   = 1./sigmapi_H1(nu_phot[0])/avg_nh1[i][j][k];
						mfp_912[i][j][k]        = 1./sigmapi_H1(13.6/h_eV)/avg_nh1[i][j][k];
						mfp_912_prev[i][j][k]   = 1./sigmapi_H1(13.6/h_eV)/avg_nh1[i][j][k];
						
						//Compute the average neutral fraction using the C2Ray procedure.  
						double col = 0.;
						if ( (coll_ion == TRUE) && (temp[i][j][k] >= 1e4) && (temp[i][j][k] <= 1e9) )  {
							col = cic_H1(temp[i][j][k])*nelec;
						}
						//solve for average HI fraction assuming constant mean photoionization rate
						double x0    = (double) (1. - f_H1[i][j][k]); //initial ionized fraction
						
						double pr;
						if (he2_tr_h2 == TRUE)  {
							pr = (double) avg_gamma[i][j][k] + col + clump[i][j][k]*(recomb_H2[i][j][k] + nhe_to_nh*recomb_He2[i][j][k])/(1. + nhe_to_nh)*nelec;
						}
						else  {
							pr = (double) avg_gamma[i][j][k] + col + clump[i][j][k]*recomb_H2[i][j][k]*nelec;
						}
						
						double xeq   = (avg_gamma[i][j][k] + col)/pr; //equilibrium ionized fraction
						double ti    = 1./pr; //timescale
						
						//average nHI over the time step
						avg_nh1[i][j][k] = (float) (nH[i][j][k]*(1. - (xeq + ti/dt*(x0 - xeq)*(1. - exp(-dt/ti))))); //average neutral fraction
						
						//Do the same thing for HeI if we are tracking it seperately.  
						double x0_he, pr_he, xeq_he, ti_he;
						if (he2_tr_h2 == FALSE)   {
							x0_he    = (double) (1. - f_He1[i][j][k]); //initial ionized fraction
							pr_he    = (double) avg_gamma_he1[i][j][k] + clump[i][j][k]*recomb_He2[i][j][k]*nelec;
							xeq_he   = avg_gamma_he1[i][j][k]/pr_he;
							ti_he    = 1./pr_he; //timescale
							avg_nhe1[i][j][k] = (float) (nHe[i][j][k]*(1. - (xeq_he + ti_he/dt*(x0_he - xeq_he)*(1. - exp(-dt/ti_he))))); //average neutral fraction
						}
						
						//On the last iteration, update number densities and treion
						if (iter == itercount - 1)  {
							//Update ionization states of H and He.  
							nH1[i][j][k]  = (float) nH[i][j][k]*(1. - (xeq + (x0 - xeq)*exp(-dt/ti)));
							nH2[i][j][k]  = (float) nH[i][j][k] - nH1[i][j][k];
							if (he2_tr_h2 == TRUE)  {
								nHe1[i][j][k] = nH1[i][j][k]/nH[i][j][k]*nHe[i][j][k];
								nHe2[i][j][k] = nH2[i][j][k]/nH[i][j][k]*nHe[i][j][k];
							}
							else  {
								nHe1[i][j][k] = (float) nHe[i][j][k]*(1. - (xeq_he + (x0_he - xeq_he)*exp(-dt/ti_he)));
								nHe2[i][j][k] = (float) nHe[i][j][k] - nHe1[i][j][k];
							}
							nHe3[i][j][k] = 0.; //no HeIII
							
							if ( (avg_nh1[i][j][k]/nH[i][j][k] < 0.5) && (treion_flag[i][j][k] == FALSE) )  {
								treion[i][j][k] = t;
								treion_out[i][j][k] = t;
								treion_flag[i][j][k] = TRUE;
							}
						}
					}
					
					//On the last iteration, update other stuff
					if (iter == itercount - 1)  {
						
						//update radiation energy density, remove later!
						for (int nu = 0; nu < Nfreq; nu++)  {
							u_nu[i][j][k][nu] = (float) (spectrum[nu]*avg_gamma[i][j][k]/sigmapi_H1(nu_phot[nu]))*h*nu_phot[nu]/clight[i][j][k];
						}
						
						for (int nu = 0; nu < Nfreq; nu++)  {
							gamma_H1_nu[i][j][k][nu]  = avg_gamma_nu[i][j][k][nu];
							gamma_He1_nu[i][j][k][nu] = avg_gamma_nu_he1[i][j][k][nu];
							gamma_He2_nu[i][j][k][nu] = 0.;
						}
						gamma_H1_tot[i][j][k]  = avg_gamma[i][j][k];
						gamma_He1_tot[i][j][k] = avg_gamma_he1[i][j][k];
						gamma_He2_tot[i][j][k] = 0.;
						//collisional ionization counts as an effective photoionization rate with coeff
						//cic x ne
						if ( (coll_ion == TRUE) && (temp[i][j][k] >= 1e4) && (temp[i][j][k] <= 1e9) )  {
							gamma_H1_tot[i][j][k]  += cic_H1(temp[i][j][k])*ne[i][j][k];
							gamma_He1_tot[i][j][k] += cic_He1(temp[i][j][k])*ne[i][j][k];
							gamma_He2_tot[i][j][k] += cic_He2(temp[i][j][k])*ne[i][j][k];
						}
						
						recrate_H2[i][j][k]  = clump[i][j][k]*recomb_H2[i][j][k]*nelec*nH2[i][j][k];
						recrate_He2[i][j][k] = clump[i][j][k]*recomb_He2[i][j][k]*nelec*nHe2[i][j][k];
						recrate_He3[i][j][k] = 0.;
						
						//Electron density and total number density
						ne[i][j][k]    = nH2[i][j][k] + nHe2[i][j][k] + 2.*nHe3[i][j][k];
						
						//MODIFIED 03/04 to give a more accurate dne/dt and ntot for the sub-grid model
						//Derivative of the electron number density (for temp calculation if needed).  
						if (subgrid == FALSE)  {
							dne_dt [i][j][k] += ne[i][j][k];
							n_tot[i][j][k] = nH[i][j][k] + nHe[i][j][k] + ne[i][j][k];
						}
						else  {
							dne_dt [i][j][k] += nH[i][j][k] + nHe[i][j][k];
							n_tot[i][j][k] = 2.*(nH[i][j][k] + nHe[i][j][k]);
						}
						dne_dt [i][j][k] /= dt; //dne_dt should be accurate in both the sub-grid model and the default mode now
						
						//Update previous values of ionization states
						//Updated 08/31/22 to fix vIF calculation in non-subgrid mode
						if (subgrid == TRUE)  {
							nH1_prev[i][j][k]  = nH1[i][j][k];
							nH2_prev[i][j][k]  = nH2[i][j][k];
							nHe1_prev[i][j][k] = nHe1[i][j][k];
							nHe2_prev[i][j][k] = nHe2[i][j][k];
							nHe3_prev[i][j][k] = nHe3[i][j][k];
						}
						
						//update ionization state fractions
						f_H2[i][j][k]  = nH2[i][j][k]/nH[i][j][k];
						f_He2[i][j][k] = nHe2[i][j][k]/nHe[i][j][k];
						f_He3[i][j][k] = nHe3[i][j][k]/nHe[i][j][k];
						
						if (ifront == FALSE)  {
							f_H1[i][j][k]  = nH1[i][j][k]/nH[i][j][k];
							f_He1[i][j][k] = nHe1[i][j][k]/nHe[i][j][k];
						}
						else {
							if (fion[i][j][k] > 0.)  {
								f_H2[i][j][k] *= fion[i][j][k];
								f_He2[i][j][k] *= fion[i][j][k];
								f_He3[i][j][k] *= fion[i][j][k];
							}
							
							f_H1[i][j][k] = (float) ((double) 1. - f_H2[i][j][k]);
							f_He1[i][j][k] = (float) ((double) 1. - f_He2[i][j][k] - f_He3[i][j][k]);
						}
					}
				}
			}
		}
		}
	}
}

//Advance a ray across the RT grid from it's initial to final position.  Break the ray into 
//segments intersecting each cell along the path, and save the distance along each segment.  
//If the C2Ray procedure is off, sum the photons deposited by each ray along the path to get
//the photo-ionization rate.  
void advance_ray(long int x)  {
	float tx, ty, tz, tmin;
	int tag;
	int i_new, j_new, k_new;
	float clgt;
	float scale_factor;
	
	//get cell indices
	int ind = ray[x].path[0].ind; 
	int i = (ind - mod(ind, Ny*Nz))/(Ny*Nz);
	int j = ind - i*Ny*Nz;
	j = (j - mod(j, Nz))/Nz;
	int k  = ind - i*Ny*Nz - j*Nz;
	
	//direction vector
	float xhat = ray[x].xhat;
	float yhat = ray[x].yhat;
	float zhat = ray[x].zhat;
	
	//positions
	float x0 = ray[x].x;
	float y0 = ray[x].y;
	float z0 = ray[x].z;
	
	//ray time
	ray[x].time += dt;
	float deltat = ray[x].time;
	float dist = clight[i][j][k]*deltat; //this is the total distance
	
	//If the ray moved during this time step, trace it's path.  
	if ( (x0 + dist*xhat > nc[i][j][k]*dx[i]) || (y0 + dist*yhat > nc[i][j][k]*dy[j])
	  || (z0 + dist*zhat > nc[i][j][k]*dz[k]) || (x0 + dist*xhat < (1. - nc[i][j][k])*dx[i]) 
	  || (y0 + dist*yhat < (1. - nc[i][j][k])*dy[j]) || (z0 + dist*zhat < (1. - nc[i][j][k])*dz[k]) )  {
		
		float r0 = ray[x].r;
		float dr = ray[x].dr;
		
		int p = 0;
		ray[x].time = 0.; 
		//MODIFIED 05/31/22 to use comoving units for r
		scale_factor = dx[0]/dx0;
		r0 += dist/scale_factor; //update radius (in comoving)
		dr -= dist/scale_factor; //update distance to next split (in comoving)
		i_new = i; //new indices
		j_new = j;
		k_new = k;
		
		float tlast = 0.;
		while (tlast < deltat)  { //loop over the path
			//Compute the time to the boundary of the next cell.  This is different for positive/negative
			//direction vector components
			clgt = clight[i_new][j_new][k_new];
			if (xhat > 0.)  { //x direction
				tx = (dx[i_new] - x0)/xhat/clgt;
			}
			else if (xhat < 0.) {
				tx = -1.*x0/xhat/clgt;
			}
			else {
				tx = 9e99;
			}
			
			if (yhat > 0.)  { //y direction
				ty = (dy[j_new] - y0)/yhat/clgt;
			}
			else if (yhat < 0.) {
				ty = -1.*y0/yhat/clgt;
			}
			else  {
				ty = 9e99;
			}
			
			if (zhat > 0.)  { //z direction
				tz = (dz[k_new] - z0)/zhat/clgt;
			}
			else if (zhat < 0.) {
				tz = -1.*z0/zhat/clgt;
			}
			else  {
				tz = 9e99;
			}
			
			//Identify which of the crossing times was smallest.  
			//tag = 0, 1, 2 --> x, y, z
			if ( (tx <= ty) && (tx <= tz) )  {
				tmin = tx;
				tag = 0;
			}
			else if ( (ty <= tx) && (ty <= tz) )  {
				tmin = ty;
				tag = 1;
			}
			else  {
				tmin = tz;
				tag = 2;
			}
			
			//If the ray still has more time on it after tmin, update the path and 
			//indices.  
			if ((tlast + tmin) < deltat)  {
				dist = tmin*clgt;
				ray[x].path[p].ind  = i_new*Ny*Nz + j_new*Nz + k_new;
				//If using the C2Ray scheme, save the distance along the segment
				if (c2ray_iter == TRUE)  {
					ray[x].path[p].dist = dist;
				}
				//Otherwise, compute the photo-ionization rate
				else  {
					double nphot[Nfreq], nphot_add[Nfreq], dtau[Nfreq];
					
					for (int nu = 0; nu < Nfreq; nu++)  {
						nphot[nu] = (double) ray[x].nphot[nu]*phot_norm;
					}
					
					//This mode is not set up to run with the sub-grid model yet, so just get
					//opacity from the grid and compute gamma from that ray
					double ddx = (double) dx[i_new];
					double ddy = (double) dy[j_new];
					double ddz = (double) dz[k_new];
					double vol = ddx*ddy*ddz;
					
					for (int nu = 0; nu < Nfreq; nu++)  {
						dtau[nu] = (double) nH1[i_new][j_new][k_new]*sigmapi_H1(nu_phot[nu])*dist;
						nphot_add[nu] = nphot[nu] * (1. - exp(-dtau[nu]));
						double gamma = nphot_add[nu]/dt/vol/nH1[i_new][j_new][k_new];
						
						//append to radiation energy 
						omp_set_lock(&(gam_lock[i_new][j_new][k_new]));
						u_nu[i_new][j_new][k_new][nu] += (float) (spectrum[nu]*gamma*h*nu_phot[nu]/clight[i_new][j_new][k_new]/sigmapi_H1(nu_phot[nu]));
						omp_unset_lock(&(gam_lock[i_new][j_new][k_new]));
						
						ray[x].nphot[nu] *= (float) exp(-dtau[nu]); //update photon count
					}
				}
				p += 1; //update ray path index and positions
				x0 += dist*xhat;
				y0 += dist*yhat;
				z0 += dist*zhat;
				
				//Check which direction the cell crossing occured in.  
				if (tag == 0)  {
					//If the crossing was in the positive direction, set x = 0 and add to index
					if (xhat > 0.)  {
						i_new = mod(i_new + 1, Nx);
						x0 = min_dist*dx[i_new]; //0.;
					}
					else  { //if negative, set x = dx and subtract from index
						i_new = mod(i_new - 1, Nx);
						x0 = dx[i_new]*(1. - min_dist); //dx[i_new];
					}
				}
				//if y direction
				else if (tag == 1)  {
					if (yhat > 0.)  {
						j_new = mod(j_new + 1, Ny);
						y0 = min_dist*dy[j_new]; //0.;
					}
					else  {
						j_new = mod(j_new - 1, Ny);
						y0 = dy[j_new]*(1. - min_dist); //dy[j_new];
					}
				}
				//if z direction
				else if (tag == 2)  {
					if (zhat > 0.)  {
						k_new = mod(k_new + 1, Nz);
						z0 = min_dist*dz[k_new]; //0.;
					}
					else  {
						k_new = mod(k_new - 1, Nz);
						z0 = dz[k_new]*(1. - min_dist); //dz[k_new];
					}
				}
				else  {
						printf("Invalid tag\n");
				}
				
				//CHRIS 07/05/22: corner correction
				if (x0 <= 0.)  {
						x0 = min_dist*dx[i_new];
				}
				if (x0 >= dx[i_new]*(1. - min_dist/2.)) {
						x0 = dx[i_new]*(1. - min_dist);
				}
				if (y0 <= 0.)  {
						y0 = min_dist*dy[j_new];
				}
				if (y0 >= dy[j_new]*(1. - min_dist/2.)) {
						y0 = dy[j_new]*(1. - min_dist);
				}
				if (z0 <= 0.)  {
						z0 = min_dist*dz[k_new];
				}
				if (z0 >= dz[k_new]*(1. - min_dist/2.)) {
						z0 = dz[k_new]*(1. - min_dist);
				}
			}
			tlast += tmin; //update time
		}
		
		//find the left-over distance traveled in the final cell, update the path and 
		//either save distance or update photo-ionization rate as above.  
		
		//MODFIED 04/05/22: if cell is partially ionized, we want the total path length through the cell, not just the path
		//during the time step, to get the number of absorptions behind the front right on average.  Otherwise we still
		//want the path during the timestep only.  This change will only affect rays that arrive at a partially ionized cell 
		//on their last path position.  
		//float t_remain = (deltat - tlast + tmin);
		if (fion[i_new][j_new][k_new] < 1.)  {
			dist = tmin*clight[i_new][j_new][k_new];
			dist = maxd(dist, min_dist*dx[i_new]); //CHRIS 07/05/22
		}
		else  {
			dist = (deltat - tlast + tmin)*clight[i_new][j_new][k_new];
			dist = maxd(dist, min_dist*dx[i_new]); //CHRIS 07/05/22
		}
		
		ray[x].path[p].ind  = i_new*Ny*Nz + j_new*Nz + k_new;
		if (c2ray_iter == TRUE)  {
			ray[x].path[p].dist = dist;
		}
		else  {
			double dtau[Nfreq], nphot[Nfreq], nphot_add[Nfreq];
			
			double ddx = (double) dx[i_new];
			double ddy = (double) dy[j_new];
			double ddz = (double) dz[k_new];
			double vol = ddx*ddy*ddz;
			
			for (int nu = 0; nu < Nfreq; nu++)  {
				nphot[nu] = (double) ray[x].nphot[nu]*phot_norm;
				dtau[nu]  = 0.;
				if ( (treion[i_new][j_new][k_new] == 0.) || (t == treion[i_new][j_new][k_new]) || (subgrid == FALSE) || (mfp[i_new][j_new][k_new] == 0.))  {
					dtau[nu] = (double) nH1[i_new][j_new][k_new]*sigmapi_H1(nu_phot[nu])*dist;
				}
				else  {
					dtau[nu] = (double) dist/mfp[i_new][j_new][k_new]; //use the mfp from the sub-grid model if the cell has already ionized
				}
			
				nphot_add[nu] = nphot[nu] * (1. - exp(-dtau[nu]));
				double gamma = nphot_add[nu]/dt/vol/nH1[i_new][j_new][k_new];
			
				omp_set_lock(&(gam_lock[i_new][j_new][k_new]));
				u_nu[i_new][j_new][k_new][nu] += (float) (spectrum[nu]*gamma*h*nu_phot[nu]/clight[i_new][j_new][k_new]/sigmapi_H1(nu_phot[nu]));
				omp_unset_lock(&(gam_lock[i_new][j_new][k_new]));
			
				ray[x].nphot[nu] *= (float) exp(-dtau[nu]);
			}
		}
		ray_end[x] = p + 1; //end of the path
		
		//update position for the last time
		x0 += dist*xhat;
		y0 += dist*yhat;
		z0 += dist*zhat;
		
		//CHRIS 06/24/22: corner correction
		if (x0 <= 0.)  {
				x0 = min_dist*dx[i_new];
		}
		if (x0 >= dx[i_new]*(1. - min_dist/2.)) {
				x0 = dx[i_new]*(1. - min_dist);
		}
		if (y0 <= 0.)  {
				y0 = min_dist*dy[j_new];
		}
		if (y0 >= dy[j_new]*(1. - min_dist/2.)) {
				y0 = dy[j_new]*(1. - min_dist);
		}
		if (z0 <= 0.)  {
				z0 = min_dist*dz[k_new];
		}
		if (z0 >= dz[k_new]*(1. - min_dist/2.)) {
				z0 = dz[k_new]*(1. - min_dist);
		}
		
		//update ray parameters
		ray[x].active = FALSE;
		ray[x].r  = r0;
		ray[x].dr = dr;
		ray[x].x  = x0;
		ray[x].y  = y0;
		ray[x].z  = z0;
		
		//If C2Ray is off, update the first path index to the new location of the ray 
		//(if C2Ray is on, this happens in the energy function).  
		if (c2ray_iter == FALSE)  {
			ray[x].path[0].ind = ray[x].path[ray_end[x] - 1].ind; //re-label 1st path index
		}
	}
}
