"""
@author: Sleeba Varghese

@description: Python script to compute the Laplace transform of a data set
and the value for Navier friction coefficient for the selected wall. 
"""

import sys
import numpy as np
import time
import concurrent.futures

# ########### User Packages ########### #
sys.path.append(" ../../custom_pkgs")
from num_integration import *
from data_analysis import Analysis
from slip import *

# ########## User Input ########## #
nrun = int(sys.argv[1])     # No. of correlation windows; 
                            # nrun = Total run time/correlation interval time 

lagTime = int(sys.argv[2])  # Correlation interval time [in ps]
skip = int(sys.argv[3])     # Frames to be skipped
irun = int(sys.argv[4])     # Number of independent simulations
nproc = int(sys.argv[5])     # Number of CPU processors

s_end = 1.0     # End value for the Laplacian variable 's'
s_num = 50002   # No. of datapoints for the Laplace transform

dt = 1.0*skip*0.001		    # value of 1 time frame(in pico-second [ps]); 
                            # Simulation dt = 1.0 femto-second	
tau = lagTime*0.001*skip
setlength = lagTime 
t_0 = 0  
t_end = setlength-1 

start_time = time.time()

# To convert LAMMPS real units [(kcal*ps)/(mol*A^2)] to SI units [(J*s)/(m^2)]
conversion_factor = 4184*1.66046*pow(10,-24)*pow(10,8)    
Area = 36.868*35.613*pow(10,-20)    # Wall surface area [m^2]

xi_M1_runs = np.zeros(irun)
xi_M2_runs = np.zeros(irun)
xi_M3_runs = np.zeros(irun)

def main(k):
	cuu_in_fname = ' ../correlation/r' + str(int(k+1)) + '_cuu.dat'
	cuf_in_fname = ' ../correlation/r' + str(int(k+1)) + '_cuf.dat'
	
	file1 = np.loadtxt(fname=cuu_in_fname)
	file2 = np.loadtxt(fname=cuf_in_fname)
	
	Time = file1[t_0:t_end,0]
	cuu = file1[t_0:t_end,1]    # <u(0).u(t)> data
	cuf = file2[t_0:t_end,1]    # <u(0).F(t)> data 
	
	frequ, Lcuu = Analysis.laplace(Time, cuu[t_0:t_end], 0.0,\
            s_end, s_num)   
	freqf, Lcuf = Analysis.laplace(Time, cuf[t_0:t_end],\
            0.0, s_end, s_num)
	
	zeta_M1 = friction_coefficient.method1(Lcuu, Lcuf)              # Method-1
	xi_M1 = (zeta_M1*conversion_factor)/Area                   	
	
	zeta_M2 = friction_coefficient.method2(frequ, Lcuu, Lcuf)       # Method-2
	xi_M2 = (zeta_M2*conversion_factor)/Area                       
	
	zeta_M3 = friction_coefficient.method3(frequ, Lcuu, Lcuf)       # Method-3
	xi_M3 = (zeta_M3*conversion_factor)/Area                       
	
	
	Lcuu_out_fname = 'r' + str(int(k+1)) + '_Lcuu.dat'          
	Lcuf_out_fname = 'r' + str(int(k+1)) + '_Lcuf.dat'          
	
	np.savetxt(Lcuu_out_fname, np.c_[frequ, Lcuu], fmt='%.10f \t %.20f')
	np.savetxt(Lcuf_out_fname, np.c_[freqf, Lcuf], fmt='%.10f \t %.20f')
	
	return k, xi_M1, xi_M2, xi_M3
	
# ########### Parallelizing ############# #
with concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as executor:
	run_ids = range(irun)
	results = executor.map(main, run_ids)
	for result in results:
		xi_M1_runs[result[0]] = result[1]
		xi_M2_runs[result[0]] = result[2]
		xi_M3_runs[result[0]] = result[3]

xi_M1_mean = np.mean(xi_M1_runs)
xi_M1_std = np.std(xi_M1_runs, ddof=1)
xi_M1_sterr = xi_M1_std/np.sqrt(irun)

xi_M2_mean = np.mean(xi_M2_runs)
xi_M2_std = np.std(xi_M2_runs, ddof=1)
xi_M2_sterr = xi_M2_std/np.sqrt(irun)

xi_M3_mean = np.mean(xi_M3_runs)
xi_M3_std = np.std(xi_M3_runs, ddof=1)
xi_M3_sterr = xi_M3_std/np.sqrt(irun)

np.savetxt('xi_M1.dat', np.c_[xi_M1_mean, xi_M1_std, xi_M1_sterr],\
        header='Mean, Std, St.Err')
np.savetxt('xi_M2.dat', np.c_[xi_M2_mean, xi_M2_std, xi_M2_sterr],\
        header='Mean, Std, St.Err')
np.savetxt('xi_M3.dat', np.c_[xi_M3_mean, xi_M3_std, xi_M3_sterr],\
        header='Mean, Std, St.Err')

end_time = time.time()
print('Finished in ' + str(round(end_time-start_time)) + ' seconds')
