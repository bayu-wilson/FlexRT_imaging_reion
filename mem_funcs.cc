#include <omp.h>
#include <stdio.h>
#include <iostream> 
#include <vector>
#include <array>
#include <algorithm>
#include <thread>
#include <pthread.h>
//#include </home/ccain002/Healpix_3.40/src/cxx/Healpix_cxx/healpix_base.h>
//#include </home/ccain002/scr16_ansond/chris/Bayu_sims/Healpix_3.40/src/cxx/optimized_gcc/include/healpix_base.h>
//#include </home/bwilson/scratch16-ansond/chris/Bayu_sims/Healpix_3.40/src/cxx/optimized_gcc/include/healpix_base.h>
#include </expanse/lustre/projects/uot171/bwils033/FlexRT/OTHER/Healpix_3.40/src/cxx/optimized_gcc/include/healpix_base.h>
#include "gas_funcs.cc"

//Count the number of rays on a single procesor. 
void update_ray_count(int nproc)  {
	ray_count[nproc] = 0;
	for (long int x = nproc*size_cpu; x < (nproc + 1)*size_cpu; x++)  {
		if (ray_tags[x] == TRUE)  {
			ray_count[nproc] += 1;
		}
	}
	ray_cap[nproc] = ray_count[nproc] / ( (float) size_cpu );
}

//update counts on all the processors
void update_ray_counts(void)  {
	#pragma omp parallel
	{
	#pragma omp for
	for (int n = 0; n < Ncpu; n++)  {
		update_ray_count(n);
	}
	}
}

long int ray_count_avg(void)  {
	long int avg = 0;
	for (int n = 0; n < Ncpu; n++)  {
		avg += ray_count[n];
	}
	avg /= Ncpu;
	return avg;
}

//initialize the parallel locks
void init_lock(void)  {
	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < Nx; i ++)   {
		for (int j = 0; j < Nx; j++)   {
			for (int k = 0; k < Nz; k++)   {
				omp_init_lock(&(ray_lock[i][j][k]));
				omp_init_lock(&(gam_lock[i][j][k]));
			}
		}
	}
	}
	
	for (int n = 0; n < Ncpu; n++)  {
		omp_init_lock(&(proc_lock[n]));
	}
	
	#pragma omp parallel
	{
	#pragma omp for
	for (int n = 0; n < 10000; n++)  {
		omp_init_lock(&(pdf_ind_lock[n]));
	}
	}
	
	omp_init_lock(&(select_lock));
}

//Select a free index from the buffer for processor nproc.  Slots are free if ray_tags[nproc] = FALSE.  
//Grab the current free slot from ray_free[nproc] and then find another slot.  Return -1 if there are
//no free slots.  
long int select_index(int nproc)  {
	long int x = ray_free[nproc];
	long int xx = x;
	bool tag = TRUE;
	while (tag == TRUE)  {
		xx += 1;
		xx = modlong(xx - nproc*size_cpu, size_cpu) + nproc*size_cpu;
		tag = ray_tags[xx];
	}
	
	if (xx == x)  {
		return -1;
	}
	else  {
		ray_tags[x] = TRUE; //flag the free spot as used
		ray_free[nproc] = xx; //assign the next free spot
		return x; //return the old free spot and assign it to a ray
	}
}

//ADDED 05/30/22: locates an occupied ray slot on processor n
long int find_occupied_ray(long int start, int nproc)  {
	long int x = start + nproc*size_cpu; //start is relative to the processor
	int tag = ray_tags[x];
	
	while (tag == FALSE)  {
		x += 1;
		x = modlong(x - nproc*size_cpu, size_cpu) + nproc*size_cpu;
		tag = ray_tags[x];
	}
	return x;
}

//remove a ray from the buffer at position x by setting it's tag to false
//reduce the ray count for cell (i, j, k)
void remove_ray(long int x)  {
	ray_tags[x]   = FALSE;
	ray[x].active = FALSE;
}

//ADDED 05/30/22 to isolate the operation of transferring a ray from one location to another
//note: this function INCLUDES the remove_ray operation, so it does not need to be done again afterwards
void move_ray(long int x_old, long int x_new)  {
	ray[x_new].active = ray[x_old].active;
	ray[x_new].order  = ray[x_old].order;
	ray[x_new].pixel  = ray[x_old].pixel;
	ray[x_new].itheta = ray[x_old].itheta;
	ray[x_new].xhat   = ray[x_old].xhat;
	ray[x_new].yhat   = ray[x_old].yhat;
	ray[x_new].zhat   = ray[x_old].zhat;
	ray[x_new].r      = ray[x_old].r;
	ray[x_new].dr     = ray[x_old].dr;
	
	for (int nu = 0; nu < Nfreq; nu++)  {
		ray[x_new].nphot[nu] = ray[x_old].nphot[nu];
	}
	
	ray[x_new].x   		    = ray[x_old].x;
	ray[x_new].y   		    = ray[x_old].y;
	ray[x_new].z   		    = ray[x_old].z;
	ray[x_new].time		    = ray[x_old].time;
	ray[x_new].path[0].ind  = ray[x_old].path[0].ind;
	ray[x_new].path[0].dist = ray[x_old].path[0].dist;
	
	remove_ray(x_old);
}

//Set the domain for merging.  At each cell, the size of the merging domain
//should be equal to the number of directions being tracked.  Not really useful now but 
//will be when I implement adaptive angular resolution.  
void set_merge_domain(void)  {
	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < Nx; i ++)   {
		for (int j = 0; j < Nx; j++)   {
			for (int k = 0; k < Nz; k++)   {
				int size = (int) 12*pow(2, 2*l_merge[i][j][k]);
				if ((int) ray_merge[i][j][k].size() != size) {
					ray_merge[i][j][k].resize(size, -1);
					ray_merge_flag[i][j][k].resize(size, FALSE);
				}
			}
		}
	}
	}
}

//rank the cpus in order of the number of rays they have on them.  
void get_cpu_ranks(void)  {
	int total = 1;
	cpu_rank[0] = 0;
	for (int n = 1; n < Ncpu; n++)  {
		cpu_rank[n] = 0;
		for (int nn = 0; nn < n; nn++)  {
			if (ray_count[nn] < ray_count[n])  {
				cpu_rank[nn] += 1;
			}
			else  {
				cpu_rank[n] += 1;
			}
		}
	}
}

//NEW SUBROUTINE
//Evenly distributes all rays to processors at the end of each RT time step
void balance_cpus(void)  {
	int top_rank    = 0;
	int bottom_rank = Ncpu - 1;
	long int x, x_new;
	long int start = 0;
	int sorted_cpu_inds[Ncpu];
	
	update_ray_counts();
	get_cpu_ranks();
	long int ray_avg = ray_count_avg();
	
	for (int n = 0; n < Ncpu; n++)  { //target rank
		printf("cpu rank %d: %d\n", n, cpu_rank[n]);
		printf("ray count %d: %d\n", n, ray_count[n]);
		for (int nn = 0; nn < Ncpu; nn++)  { //cpu index
			if (cpu_rank[nn] == n)  {
				sorted_cpu_inds[n] = nn;
				break;
			}
		}
		printf("sorted cpu ind %d: %d\n", n, sorted_cpu_inds[n]);
	}
	
	while (top_rank < bottom_rank)  {
		x_new = select_index(sorted_cpu_inds[bottom_rank]);
		x = find_occupied_ray(start, sorted_cpu_inds[top_rank]);
		start = x - sorted_cpu_inds[top_rank]*size_cpu;
		
		//printf("ray tag (should be true): %d\n", ray_tags[x]);
		
		move_ray(x, x_new);
		
		ray_count[sorted_cpu_inds[top_rank]]    -= 1;
		ray_count[sorted_cpu_inds[bottom_rank]] += 1;
		
		//printf("ray count on processor %d (top): %d\n", sorted_cpu_inds[top_rank], ray_count[sorted_cpu_inds[top_rank]]);
		//printf("ray count on processor %d (bottom): %d\n", sorted_cpu_inds[bottom_rank], ray_count[sorted_cpu_inds[bottom_rank]]);
		
		if (ray_count[sorted_cpu_inds[top_rank]] <= ray_avg)  {
			top_rank += 1;
			//printf("top rank: %d\n", top_rank);
		}
		
		if (ray_count[sorted_cpu_inds[bottom_rank]] >= ray_avg)  {
			bottom_rank -= 1;
			start = 0;
			//printf("bottom rank: %d\n", bottom_rank);
		}
	}
	
}

//if a cpu is more than 95% full, free up some space on it by locating the processor with the fewest rays until the
//fraction of used space drops below 95%.  
void free_cpu(int nproc)  {
	printf("freeing cpu %d\n", nproc);
	long int x = nproc*size_cpu;
	int free_proc;
	for (int n = 0; n < Ncpu; n++)  {
		if (cpu_rank[n] == Ncpu - 1)  {
			free_proc = n;
		}
	}
	
	//MODIFIED 04/10/21 to include processor lock on select_index routine
	//as long as the full cpu is over max_cpu_cap % full, move rays to the cpu that initially had the 
	//fewest rays.  
	long int x_new = 0;
	while ( (ray_cap[nproc] > max_cpu_cap) && (x_new != -1) )  {
		
		//assume the lowest-ranked processor is not full...obviously
		omp_set_lock(&(proc_lock[free_proc]));
		x_new = select_index(free_proc);
		omp_unset_lock(&(proc_lock[free_proc]));
		
		bool tag = FALSE;
		omp_set_lock(&(proc_lock[nproc]));
		while (tag == FALSE)  {
			x += 1;
			x = modlong(x - nproc*size_cpu, size_cpu) + nproc*size_cpu;
			tag = ray_tags[x];
		}
		omp_unset_lock(&(proc_lock[nproc]));
		
		//copy ray to the new location
		ray[x_new].active = ray[x].active;
		ray[x_new].order  = ray[x].order;
		ray[x_new].pixel  = ray[x].pixel;
		ray[x_new].itheta = ray[x].itheta;
		ray[x_new].xhat   = ray[x].xhat;
		ray[x_new].yhat   = ray[x].yhat;
		ray[x_new].zhat   = ray[x].zhat;
		ray[x_new].r      = ray[x].r;
		ray[x_new].dr     = ray[x].dr;
		
		for (int nu = 0; nu < Nfreq; nu++)  {
			ray[x_new].nphot[nu] = ray[x].nphot[nu];
		}
		
		ray[x_new].x   		    = ray[x].x;
		ray[x_new].y   		    = ray[x].y;
		ray[x_new].z   		    = ray[x].z;
		ray[x_new].time		    = ray[x].time;
		ray[x_new].path[0].ind  = ray[x].path[0].ind;
		ray[x_new].path[0].dist = ray[x].path[0].dist;
		
		remove_ray(x);
		
		ray_count[nproc] -= 1;
		ray_cap[nproc] = ray_count[nproc] / ( (float) size_cpu );
		ray_count[free_proc] += 1;
	}
	ray_cap[free_proc] = ray_count[free_proc] / ( (float) size_cpu );
}

//functions for summing over all processors with floats 
float sum_over_procs(float f[Ncpu])  {
	float sum = 0.;
	for (int n = 0; n < Ncpu; n++)  {
		sum += f[n];
	}
	return sum;
}

//and doubles
double sum_over_procsd(double f[Ncpu])  {
	double sum = 0.;
	for (int n = 0; n < Ncpu; n++)  {
		sum += f[n];
	}
	return sum;
}
