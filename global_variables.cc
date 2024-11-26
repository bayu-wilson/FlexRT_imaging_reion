#include <math.h>
#include <vector>
#include <omp.h>
#include <pthread.h>
#include "user_inputs.h"

struct path_type  { //ray path structure
	int ind;
	float dist;
};

struct source_type  {
	int i;
	int j;
	int k;
	double lum;
};

//data structure for storing rays
struct ray_type  { //1 + 2 + 2 + 4*19 = 81 bytes
	bool active;
	short int order;
	int pixel;
	unsigned short int itheta;
	float xhat;
	float yhat;
	float zhat; //maybe can convert to two angles
	float nphot[Nfreq];
	float r;
	float dr; //might be able to drop
	float x;
	float y;
	float z;
	float time;
	path_type path[path_length];
};

float dt;

short int const ray_per_cell_max = 22; //maximum number of rays per cell
short int const ray_buff	 	= 8;
short int const size_cell		= ray_per_cell_max + ray_buff; //maximum size each cell can hold in memory
long int const size_buffer		= ((long int) Nx*Ny*Nz) * size_cell; //buffer size
long int const size_cpu 		= size_buffer / Ncpu; //memory per cpu
int num_directions 				= 12*pow(2, 2*hpx_lvl); //number of directions tracked
float max_cpu_cap = 0.95;

long int size_cache = 750;

//initialize order for merging, and min/max order (max = 13 always)
//short int l_merge 	= hpx_lvl;
short int l_merge[Nx][Ny][Nz];
short int l_source[Nx][Ny][Nz];
short int ray_per_cell[Nx][Ny][Nz];
short int l_max 	= 12; //to try to avoid Planck error, set below 13
short int l_min 	= hpx_lvl;

//smoothing length for the gamma smoothing algorithm (not using currently)
int smooth_length = 4;

//time & stepping variables
float t 		= 0.;
float t_step	= 0.;
int    step		= 0;
int   srcstep	= 0; //testing
float tinit;

//ADDED 05/31/22 - cell size in comoving units
float dx0 = Lx0/Nx*kpc_to_cm;
float dy0 = Ly0/Ny*kpc_to_cm;
float dz0 = Lz0/Nz*kpc_to_cm;

//speed of light, cell crossings per update, and time since last ray cast from sources
float clight	[Nx][Ny][Nz];
float nc		[Nx][Ny][Nz];
float ray_dt    [Nx][Ny][Nz];

//grid variables
float dx	[Nx]; //x differential
float dy	[Ny]; //y differential
float dz	[Nz]; //z differential

//number densities
float avg_rho = 0.;
float nhe_to_nh = (Y/4.)/(1. - Y);
float  rho		[Nx][Ny][Nz];
float  nH		[Nx][Ny][Nz]; //number density of hydrogen
float  nHe		[Nx][Ny][Nz]; //number density of helium
float  nH1		[Nx][Ny][Nz]; //number density of neutral hydrogen
float  nHe1		[Nx][Ny][Nz]; //number denstiy of neutral helium
float  nH2		[Nx][Ny][Nz]; //number density of ionized hydrogen //can remove
float  nHe2		[Nx][Ny][Nz]; //number density of singly ionized helium
float  nHe3		[Nx][Ny][Nz]; //number density of doubly ionized helium //can remove
float  ne		[Nx][Ny][Nz]; //number density of free electrons
float  n_tot	[Nx][Ny][Nz]; //total number density

/////////////////////////////////////////////
//ADDED 05/27/22, Fahad's changes FlexRT to work with MP-GADGET
//Added by Fahad Nasir 09/16/2021
//ergs/s/grams
float  gadget_ini_heat  [Nx][Ny][Nz]; //the photoheating rates + the heat injection rate
//[erg/s] 3/2 Delta_fion kB Treion(vIF)/(mu mP dt)  required by MP-Gadget output files
float  gadget_heat      [Nx][Ny][Nz];
float  gadget_gamma     [Nx][Ny][Nz];
float gadget_step_time;
/////////////////////////////////////////////

//previous time step
float rho_prev  [Nx][Ny][Nz];
float nH1_prev	[Nx][Ny][Nz]; //number density of neutral hydrogen
float nHe1_prev	[Nx][Ny][Nz]; //number denstiy of neutral helium
float nH2_prev	[Nx][Ny][Nz]; //number density of ionized hydrogen //can remove
float nHe2_prev	[Nx][Ny][Nz]; //number density of singly ionized helium
float nHe3_prev	[Nx][Ny][Nz]; //number density of doubly ionized helium //can remove

//time derivatives
float dne_dt	[Nx][Ny][Nz];
float ddelta_dt [Nx][Ny][Nz];

//abundances (can remove all "steps")
float f_H1		[Nx][Ny][Nz]; //HI fraction
float f_H1_step [Nx][Ny][Nz];
float f_He1	 	[Nx][Ny][Nz]; //HeI fraction
float f_He1_step[Nx][Ny][Nz];
float f_H2		[Nx][Ny][Nz]; //HII fraction
float f_H2_step [Nx][Ny][Nz];
float f_He2	 	[Nx][Ny][Nz]; //HeII fraction
float f_He2_step[Nx][Ny][Nz];
float f_He3	 	[Nx][Ny][Nz]; //HeIII fraction
float f_He3_step[Nx][Ny][Nz];

//thermal parameters
float temp	[Nx][Ny][Nz]; //temperature

//photoionization rates
float gamma_H1_nu	[Nx][Ny][Nz][Nfreq]; //photoionization rate of HI
float gamma_He1_nu 	[Nx][Ny][Nz][Nfreq]; //photoionization rate of HeI
float gamma_He2_nu	[Nx][Ny][Nz][Nfreq]; //photionization rate of He2
float gamma_H1_tot	[Nx][Ny][Nz]; //photoionization rate of HI
float gamma_He1_tot [Nx][Ny][Nz]; //photoionization rate of HeI
float gamma_He2_tot [Nx][Ny][Nz]; //photionization rate of He2

//recombination coefficients
float recomb_H2 	[Nx][Ny][Nz];
float recomb_He2	[Nx][Ny][Nz];
float recomb_He3	[Nx][Ny][Nz];

//recombination rate
float recrate_H2        [Nx][Ny][Nz];
float recrate_He2       [Nx][Ny][Nz];
float recrate_He3       [Nx][Ny][Nz];

//heating/cooling rates
float heat_rate	[Nx][Ny][Nz]; //total heating rate
float cool_rate	[Nx][Ny][Nz]; //total cooling rate

//energy and spectrum
float spectrum 		[Nfreq];
float u_nu			[Nx][Ny][Nz][Nfreq];  //specific energy density of radiation

/////////////////////////////////////////////////////////////////////////////////
//RT variables for ray tracing program

long int total_rays = 0;
long int num_thresh = 0;
double phot_norm = 1e60;
double nphot_base_avg = 0.;
double nphot_min;
double nphot_max;
float nphot_merge;

int num_src, num_inc;
vector<struct source_type> source_lum; //luminosity of source cells.  0 for cells that are not a source  

struct				ray_type ray[size_buffer]; //ray buffer
bool				ray_tags [size_buffer]; //tags for positions within cell
vector<long int> 	ray_merge[Nx][Ny][Nz]; //domain vectors for ray merging (size is adjustable)
vector<bool> 		ray_merge_flag[Nx][Ny][Nz]; //keep track of whether merging happened in a cell
short int			ray_end[size_buffer]; //keep track of the index where each ray path terminates
long int 			ray_free[Ncpu];
long int 			ray_count[Ncpu];
float 				ray_cap[Ncpu];
short int 			cpu_rank[Ncpu];

int ray_move_counter[Ncpu]; //count how many rays each processor moved

//locks to prevent race conditions
omp_lock_t select_lock;
omp_lock_t proc_lock[Ncpu];
omp_lock_t pdf_ind_lock[10000];
omp_lock_t ray_lock[Nx][Ny][Nz];
omp_lock_t gam_lock[Nx][Ny][Nz];

//average gamma and xh1 for C2Ray scheme and treion
float avg_nh1[Nx][Ny][Nz];
float avg_nhe1[Nx][Ny][Nz];
float treion[Nx][Ny][Nz];
float treion_out[Nx][Ny][Nz];
float avg_gamma[Nx][Ny][Nz];
float avg_gamma_prev[Nx][Ny][Nz];
float avg_gamma_he1[Nx][Ny][Nz];
double ngam[Nx][Ny][Nz][3];

//some flags
bool treion_flag[Nx][Ny][Nz];
bool tcross_flag[Nx][Ny][Nz];
bool neutral_flag[Nx][Ny][Nz];

//testing
float avg_gamma_nu	[Nx][Ny][Nz][Nfreq];
float avg_gamma_nu_he1	[Nx][Ny][Nz][Nfreq];
float ngam_nu	  	[Nx][Ny][Nz][Nfreq];

bool rh_init[Nx][Ny][Nz];

//testing
float fion[Nx][Ny][Nz];
float fion_prev[Nx][Ny][Nz];
float tcross[Nx][Ny][Nz];
float tcross_prev[Nx][Ny][Nz];
float gamma_rh_ss[Nx][Ny][Nz];

//Added 08/19/22 - new zreion averaging scheme
int const Ntcross = 10; //cannot be larger than 100
float tcross_new[Nx][Ny][Nz][Ntcross];
bool  tcross_check[Nx][Ny][Nz][Ntcross];
float fion_max[Nx][Ny][Nz];

//Stuff for Bayu's project - vIF and incident flux.  Added 08/31/22
float IF_speed[Nx][Ny][Nz];
float Inc_Flux[Nx][Ny][Nz];
float Delta_fion[Nx][Ny][Nz];

float xh1_thresh = 0.01;

//////////////////////////////////
//Sub-grid model stuff

float gamma_thresh = 3e-16;
float t_res = 15.; //time after zreion to calculate the residual HI. 
float sigma_eff_rh = 2.517139346618065e-18;
float nu_rh[5] = {14.48/h_eV, 16.70/h_eV, 20.05/h_eV, 25.84/h_eV, 39.50/h_eV};
float nu_eff_rh = (14.48*sigmapi_H1(14.48/h_eV) + 16.70*sigmapi_H1(16.70/h_eV) + 20.05*sigmapi_H1(20.05/h_eV) + 25.84*sigmapi_H1(25.84/h_eV) + 39.50*sigmapi_H1(39.50/h_eV))/5./h_eV/sigma_eff_rh;
float beta_col = 1.3; //assumed slope of the column density distribution
float max_subgrid_frac = 0.49;

float mfp_boost = 0.; //factor by which the MFP should be boosted over the 912 value in the subgrid model

//clumping factor and mean free path
float clump[Nx][Ny][Nz];
float mfp  [Nx][Ny][Nz];
float mfp_prev[Nx][Ny][Nz];
float mfp_next[Nx][Ny][Nz];

//Source opacity subgrid model
double source_tau[Nx][Ny][Nz];

//MFP at 912 angstroms
float mfp_912     [Nx][Ny][Nz];
float mfp_912_prev[Nx][Ny][Nz];

//interpolation arrays for Radhydro data
const int Nzre = 3;
const int Ndelta = 3;
const int Ngamma = 3;
const int Nt = 5000;

float zre_dom[Nzre], delta_dom[Ndelta], deltaNL_dom[Ndelta], gamma_dom[Ngamma+2];
float tre_dom[Nzre];

//radhydro variables (renamed 08/19/22)
float time_RH[Nt][Nzre][Ndelta][Ngamma];
float mfp_RH [2][Nt][Nzre][Ndelta][Ngamma]; //0 = frequency-averaged, 1 = at 912 angstroms
float fh1_RH [Nt][Nzre][Ndelta][Ngamma];
float temp_RH[Nt][Nzre][Ndelta][Ngamma];

//gamma interpolation arrays
float mfp_gamma_interp_b  [2][Nx][Ny][Nz][Ngamma+2];
float mfp_gamma_interp_f  [2][Nx][Ny][Nz][Ngamma+2];
float fh1_gamma_interp    [Nx][Ny][Nz][Ngamma+2];
float temp_gamma_interp   [Nx][Ny][Nz][Ngamma+2];

//time slope for i-front correction
float mfp_t_tcross_gamma_interp_b    [2][Nx][Ny][Nz][Ngamma+2];
float mfp_t_tcross_gamma_interp_f    [2][Nx][Ny][Nz][Ngamma+2];
float fh1_t_tcross_gamma_interp      [Nx][Ny][Nz][Ngamma+2];
float temp_t_tcross_gamma_interp     [Nx][Ny][Nz][Ngamma+2];

//Added 08/19/22: use tcross bins
//float mfp_t_tcross_new_gamma_interp_b    [5][2][Nx][Ny][Nz][Ngamma+2];
//float mfp_t_tcross_new_gamma_interp_f    [5][2][Nx][Ny][Nz][Ngamma+2];
//float fh1_t_tcross_new_gamma_interp      [5][Nx][Ny][Nz][Ngamma+2];
//float temp_t_tcross_new_gamma_interp     [5][Nx][Ny][Nz][Ngamma+2];

//for delta/sigma <-> deltaNL conversion
const int n_dlin = 500;
float t_deltaNL[n_dlin];
float delta_over_sigma_interp[n_dlin];
float deltaNL_interp[n_dlin][n_dlin];

//minimum temperature and MFP
float mfp_min  = 1e-2*dx[0];
float gamma_min = 1e-16;
float temp_min = 1e2;
float min_dist = 1e-4; //CHRIS 07/05/22

//hydro variables
int ihydro = 0;
int num_hydro_steps;
float Lx, Ly, Lz;
float t_tot = 0.;
float zz = 99.;
float zz_end;
vector<float> hydro_steps;

//Mao+19 clumping arrays
int const n_mao = 76;
float mao_z[n_mao];
float mao_a0[n_mao];
float mao_a1[n_mao];
float mao_a2[n_mao];

//spectrum stuff
int const Nalph = 15;
string spectrum_table = "input_files/spectrum_tables_5bins.txt";
double alpha[Nalph];
double spectrum_input[Nalph][Nfreq];
