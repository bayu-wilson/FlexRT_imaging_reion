import numpy as np
import pp_vars_and_funcs as ppvf 
from scipy.interpolate import griddata

###### CUBES ######
### User input 
alpha = 1.5
Ngas = ppvf.Ngas #400
### loading vIF and incident LyC radiation
base_path = "/expanse/lustre/projects/uot171/bwils033/FlexRT/output_files/"
#sim_path = "start_alpha1.5_monoFreq_maxrapid_N400_f0.75_cl0.3/" 
#IF_info = "vIF_z=05.8029"
#sim_path = "restart_alpha1.5_monoFreq_maxrapid_N400_f0.4_cl0.3/" 
#IF_info = "vIF_z=05.7767"
sim_path = "restart_alpha1.5_monoFreq_maxrapid_N400_f1.5_cl0.3/"
IF_info = "vIF_z=05.8913" #"vIF_z=05.8624"
#IF_info = "vIF_z=05.7767" #05.8029" #"05.7767" #"vIF_z=05.7701"#05.7832" #"vIF_z=05.7996" ####MAKE
#IF_info = "vIF_z=05.8029"
#sim_path = "restart_alpha1.5_monoFreq_maxrapid_N400_f0.4_cl0.3/"

path_to_IF_info = base_path + sim_path + IF_info

### loading interpolation table
interp_table = np.loadtxt("interp_tables/fd_parameter_space_interp_table.txt")
bb = np.column_stack((interp_table[:,0],interp_table[:,1]))
Lya_efficiency = interp_table[:,2]



vIF_moving_screen = np.zeros((Ngas,Ngas,Ngas))
F_lyC_inc = np.zeros((Ngas,Ngas,Ngas))
xHII = np.zeros((Ngas,Ngas,Ngas))
zlist = ["vIF_z=05.8891","vIF_z=05.8902","vIF_z=05.8924","vIF_z=05.8935","vIF_z=05.8913"]
for i in range(len(zlist)):
    path_i = base_path + sim_path + zlist[i]
    print(f"loading: {path_i}")
    IF = np.fromfile(path_i,dtype=np.float32).reshape((4,Ngas,Ngas,Ngas),order='F')
    vIF_moving_screen += IF[0]/1e5  #km/s
    F_lyC_inc += IF[1] #ph/s/cm2
    xHII += IF[2]
#xHII = IF[2]
xHII /= len(zlist)
mask_IF = xHII < 0.99
vIF_moving_screen *= mask_IF / len(zlist) #divide by three to get the average 
F_lyC_inc *= mask_IF / len(zlist)

np.save("results/May_xHII.npy",xHII)

"""
print(f"loading: {path_to_IF_info}")
IF = np.fromfile(path_to_IF_info,dtype=np.float32).reshape((4,Ngas,Ngas,Ngas),order='F')
vIF_moving_screen = IF[0]/1e5  #km/s
F_lyC_inc = IF[1] #ph/s/cm2
xHII = IF[2]
mask_IF = xHII < 0.99
print(len(xHII.flatten()))
print(np.sum(mask_IF.flatten()),flush=True)
vIF_moving_screen *= mask_IF
F_lyC_inc *= mask_IF
#xHII *= mask_IF
#vIF_array = vIF_array[(fion_array[mask_vIF]<0.99)]
"""
Flya=griddata(bb,Lya_efficiency,(np.log10(vIF_moving_screen),alpha),method='linear',fill_value=0)*F_lyC_inc
print("FLya:{}".format(np.nanmedian(Flya)))
pi = np.pi
E_lya = ppvf.E_lya
SB = Flya*E_lya/4/np.pi/(1.+float(ppvf.z))**4./ppvf.sr_to_arcsec2 #erg/s/cm2/arcsec 

### saving data cubes
ff = sim_path.split("_")[-2][1:]
np.save(f"results/dataCubes/vIF_a{alpha}_f{ff}.npy", vIF_moving_screen)
np.save(f"results/dataCubes/FlyC_a{alpha}_f{ff}.npy", F_lyC_inc)
np.save(f"results/dataCubes/Flya_a{alpha}_f{ff}.npy", Flya) #no redshift dimming
np.save(f"results/dataCubes/SB_a{alpha}_f{ff}.npy", SB) #yes redshift dimming

###### SLICES ######
#get indices corresponding to the center of the filter
index = int(ppvf.nb_depth/ppvf.sim_depth*Ngas+0.5)
index_half = int(index/2)

### get gas file
#gas_info = "gas_z=05.7882" #05.7882" #"gas_z=05.8824"
gas_info = "gas_z=05.8349" #"gas_z=05.8349"
path_to_gas_info = base_path + sim_path + gas_info

print(f"loading: {path_to_gas_info}")
gas = np.fromfile(path_to_gas_info,dtype=np.float32)[1:].reshape((18+1,Ngas,Ngas,Ngas),order='F')#23/01/05 BAYU: include +5 for multifrequency
rho = gas[0]
Delta = rho/ppvf.rhoBaryons_0

#choose slice
Delta_slice = (Delta)[index_half,:,:]#[:,:,index_half]
xHII_slice = xHII[index_half,:,:] #[:,:,index_half] 
F_lyC_inc_slice = F_lyC_inc[index_half,:,:] #[:,:,index_half] 
vIF_slice = vIF_moving_screen[index_half,:,:]      #[:,:,index_half] 
F_lya_slice = Flya[index_half,:,:]

stack = np.stack((Delta_slice, xHII_slice, F_lyC_inc_slice, vIF_slice,F_lya_slice))
np.save("results/slices/FEB_slice_stack.npy",stack)


###### RECOMBINATIONS AND AVERAGES ######
#gas suffix because slightly different redshift than vIF file
xHII_gas = gas[5]
xHeII_gas = gas[6]
xHeIII_gas = gas[7]
m_H = ppvf.m_H
Y = ppvf.Y
nH = (1-Y)*rho/m_H
nHe = Y*rho/m_H
nHII = nH*xHII_gas
#ne = nH*xHII + nHe*(xHeII_gas+xHeIII_gas)
ne = nHII + nHe*(xHeII_gas+xHeIII_gas)
T_OOM = 1.0e+4
alpha_eff_lya = ppvf.epsilonB_lya(T_OOM)*ppvf.alphaB_H2(T_OOM)
j_lya_recombB = nHII * ne * alpha_eff_lya * E_lya/4/pi / (1.+ppvf.z)**4
RR = j_lya_recombB*ppvf.cell_size_proper #erg/s/cm2/sr
np.save(f"results/dataCubes/RR_a{alpha}_f{ff}.npy", RR)

avg_xHI = np.mean(1-xHII[:,0:index,:],axis=1)
np.save("results/averages/FEB_avg_neutral_fraction.npy", avg_xHI)

avg_Delta = np.mean(Delta[:,10:10+index,:],axis=1)
np.save("results/averages/FEB_avg_Delta.npy", avg_Delta)



