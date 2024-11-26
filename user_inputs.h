#include "rates.cc"

#define FALSE 0
#define TRUE 1

char const rec_case[]   = "B"; //Use case "A" or "B" recombination coefficients

int const input_grid    = TRUE; //User-defined input grid flag

bool const c2ray_iter    = TRUE; //implement C2Ray iteration scheme flag
bool const ifront        = TRUE; //track i-fronts
bool const he2_tr_h2     = FALSE; //Set xHeII = xHII flag
bool const subgrid       = TRUE; //sub-grid model flag 
bool const relaxed       = FALSE; //use relaxed sub-grid limit flag
bool const equilibrium   = TRUE; //use an equilibrium subgrid model
bool const recomb_cool   = TRUE; //Recombination cooling flag
bool const coll_ion      = FALSE; //Collisional ionization flag
bool const coll_exc_cool = FALSE; //Collisional excitation cooling flag
bool const free_free     = TRUE; //free-free cooling flag
bool const compton       = TRUE; //compton heating/cooling flag
bool const temp_ev       = TRUE; //temperature evolution flag

//Chris' multi frequency additions. Added by Bayu Jan. 4 2023
bool const evolv_spec = FALSE;
bool const multi_freq = FALSE;

//MP-GADGET (ADDED 05/27/22)
bool const MP_GADGET_HEAT = FALSE;
int const gadget_step     = 1;

//initial source and gas files (always set these)
const int snapinit             = 56;
///home/bwilson/scratch16-ansond/chris/Bayu_sims/
//const char ugas_file_init[] ="/home/bwilson/scratch4-ansond/chris/P3M_Nbody/subres_sources/input_files_400/density.000.400";
const char ugas_file_init[] ="/expanse/lustre/projects/uot171/bwils033/FlexRT/OTHER/subres_sources/input_files_400/density.056.400";
//BAYU: 200 is cell per side, change everwhere
//const char source_field_init[] = "/home/bwilson/scratch4-ansond/chris/P3M_Nbody/subres_sources/ndot_catalogs_40_3e8_maxrapid/ndot_catalog_000";  //BAYU: emissivity catalog, where source model/ emissivity history is
const char source_field_init[] = "/expanse/lustre/projects/uot171/bwils033/FlexRT/OTHER/subres_sources/ndot_catalogs_40_3e8_maxrapid/ndot_catalog_056";  //BAYU: emissivity catalog, where source model/ emissivity history is
//BAYU: emissivity catalog, where source model/ emissivity history is
const char hydro_steps_file[]  = "./hydro_Snapshots_200_restart.txt"; 
//const char num_src_file_init[] = "/home/bwilson/scratch4-ansond/chris/P3M_Nbody/subres_sources/ndot_catalogs_40_3e8_maxrapid/numhalos.000";//BAYU folder might change
const char num_src_file_init[] = "/expanse/lustre/projects/uot171/bwils033/FlexRT/OTHER/subres_sources/ndot_catalogs_40_3e8_maxrapid/numhalos.056";//BAYU folder might change
//BAYU folder might change // /home/bwilson/scratch16-ansond/chris/P3M_Nbody/subres_sources/ndot_catalogs_40_3e8_rapid25/numihalos.000

//input files for a restart (if needed)
//alpha1.5_monoFreq_maxrapid_N400_f0.75_cl0.3 or start_alpha1.5_monoFreq_rapid25_N400_f0.65_cl0.3
//start_alpha1.5_monoFreq_rapid25_N400_f0.65_cl0.25/ 06.1328
//start_alpha1.5_monoFreq_rapid25_N400_f0.5_cl0.3/ray_z=06.2394
const char start_file_init[]   = "./output_files/start_alpha1.5_monoFreq_maxrapid_N400_f0.75_cl0.3/gas_z=06.0809";
const char ray_file_init[]     = "./output_files/start_alpha1.5_monoFreq_maxrapid_N400_f0.75_cl0.3/ray_z=06.0809";
const char treion_file_init[]  = "./output_files/start_alpha1.5_monoFreq_maxrapid_N400_f0.75_cl0.3/tre_z=06.0809";

//output directory
const char output_dir[]      = "./output_files/restart_alpha1.5_monoFreq_maxrapid_N400_f1.5_cl0.3/";

//more directories - ADDED 05/28/22
//const char hydro_base[]   = "/home/bwilson/scratch4-ansond/chris/P3M_Nbody/subres_sources/";
const char hydro_base[]    = "/expanse/lustre/projects/uot171/bwils033/FlexRT/OTHER/subres_sources/";
//"/home/ccain002/scr16_ansond/chris/P3M_Nbody/subres_sources/"; //BAYU folder might change
const char density_dir[] = "input_files_400/"; //BAYU folder might change
const char source_dir[] = "ndot_catalogs_40_3e8_maxrapid/"; //BAYU folder might change
const char subgrid_dir[]  = "subgrid_files/";

//file bases
const char density_base[]   = "density.";
const char source_base[]    = "ndot_catalog_";
const char sourcenum_base[] = "numhalos.";

//on-the-fly output file
const char otf_output[]      = "output_files/restart_alpha1.5_monoFreq_maxrapid_N400_f1.5_cl0.3/otf.txt";
//"output_files_f0.7/otf_test_f0.7.txt"; //BAYU change for every run

//grid sizes
int const Ncpu    = 48; //number of threads
int const Nx      = 400; //Number of grid points in the x direction
int const Ny      = 400; //Number of grid points in the x direction
int const Nz      = 400; //Number of grid points in the x direction
int const Nfreq = 1; //Number of frequency bins 

//time stepping  //BAYU CHANGED FROM 0 to 2 and back to 0
short int const hpx_lvl = 0; //healpix level for ray casting/merging /BAYU angular resolution parameter. controls merging
float const tstep_factor = 1.0; //fraction of light crossing time per time step
float const nl_const     = 0.5; //update and casting frequency parameter
float const cl_factor    = 0.3; //simulation speed of light / true speed of light 
int const itercount      = 5; //number of iterations for the photo-ionization rate solver
float const safety_factor = 2.; //minimum number of rays per cell //increase safety factor will increase angular res. splitting.
float const const_clump = 1.0; //constant clumping factor (in subgrid mode this scales down the MFP even more) 

float const zinit    = 12.0038; //initial redshift

float const fesc                = 1.5;//escape fraction for sources //BAYU: THIS IS f in the filename.
float const Lx0      		= 40*1e3/hh; //Box length in the x direction (in kpc)
float const Ly0      		= 40*1e3/hh; //Box length in the y direction
float const Lz0      		= 40*1e3/hh; //Box length in the z direction
double const nu_phot[Nfreq] = {19/h_eV} ; //Photon frequency (in Hz) //BAYU: CHANGE THIS FOR DIFFERENT ALPHA 
float const temp_0  		= 1e2; //Initial gas temperature/temperature of neutral gas

float const fH1_0   = 1.- 1.2e-3; //Initial HI fraction
float const fHe1_0  = 1.- 1.2e-3; //Initial HeI fraction
float const fHe2_0  = 1.2e-3; //Initial HeII fraction

float const alpha_mfp = 0.67; //Power law index of MFP vs. Gamma (if in sub-grid mode)
float const t_relax = 100.; //Relaxtion timescale in Myr

int const NumStep = 1; //Number of steps between on the fly output writes (set to 1 for cosmo sims)

int const path_length = 4; //change for c2ray iteration mode

