import numpy as np

#user input
Ngas = 400
z = 5.7
L_cMpc_over_h = 40

nb_depth = 26 #mpc/h
sim_depth = 40 #mpc/h
#index = int(nb_depth/sim_depth*Ngas+0.5)
#index_half = int(index/2)
omega_b = 0.048
rho_crit_z = 8.688e-30*(1+z)**3 #g/cm3
rhoBaryons_0 = rho_crit_z*omega_b
Y=0.24
m_H = 1.672622e-24 #g


#constants
h = 6.626068e-27 #cgs
c = 2.998e10 #cm/s
lambda_lya_cm = 1.21567e-5 #cm
E_lya = h*c/lambda_lya_cm #erg
sr_to_arcsec2 = 4.25e10
kpc2cm = 3.08560e21
Mpc2km = 3.08560e19
yr2sec = 3.154e+7
Myr2sec = 1e+6*yr2sec

#cosmology
hh = 0.68

cell_size_proper = L_cMpc_over_h/Ngas*1000/hh/(1+z)*kpc2cm

def h_of_a(a, h=0.68, omega_m=0.305147, t_cmb0=2.7255):
    # Assume flat universe.
    omega_r = 4.48e-7 * (1 + 0.69) * t_cmb0 ** 4 / h ** 2
    omega_l = 1 - omega_m - omega_r
    # Calculate hubble parameter squared at scale factor a.
    hsq = (100 * h) ** 2 * (omega_m / a ** 3 + omega_r / a ** 4 + omega_l)
    return hsq ** 0.5
# Get the cosmic time since the big bang.  This function lets us work in the time domain when interpolating the
# clumping factor and its derivative.
def time_global_of_z(h=0.68):
    # Define a redshift domain the runs from the big bang to z = 5
    z = np.logspace(5, 0, int(1e6))
    a = 1 / (1 + z)
    # Integrate to get cosmic time
    time = cumint(1 / (a * h_of_a(a, h=h) / Mpc2km), a)
    # Return an interpolation function.
    return interp1d(z, time)

### FOR RECOMBINATIONS ###
k_B = 1.38065e-16
eV_cgs = 1.60218e-12
Eth_H1 =1.360e1

def llambda(T, Eth):
    K2eV = k_B/eV_cgs
    eV2K = 1.0/K2eV
    xHI  = Eth*eV2K/T
    x    = 2*xHI
    return x

def alphaA_H2(T):
    x = llambda(T, Eth_H1)
    return 1.269e-13*pow(x,1.503)/pow((1.0+pow((x/0.522),0.470)),1.923)

#Case B recombination coefficient of HII from Hui & Gnedin 1996
def alphaB_H2(T):
     x = llambda(T, Eth_H1)
     return 2.753e-14*pow(x,1.500)/pow((1.0+pow((x/2.740),0.407)),2.2420)

#fraction of recombinations producing lya phootons
def epsilonB_lya(T):
    T4 = T/1e+4
    return 0.686 - 0.106*np.log10(T4) - 0.009*T4**(-0.44)


