"""
@author: Sleeba Varghese

@description: Python script to compute the 
time correlation between two data sets.
"""

import sys						
import numpy as np
import time
import concurrent.futures

# ########### User Packages ######### #
sys.path.append(" ../../custom_pkgs")
from data_analysis import Analysis

# ########### User Input ########## #
nrun = int(sys.argv[1])         # No. of correlation windows; 
                                # nrun = Total run time/correlation interval time 
startTime = int(sys.argv[2])    # Start frame of post-process stage
endTime = int(sys.argv[3])      # End frame of post-process stage
skip = int(sys.argv[4])         # Frames to be skipped
irun = int(sys.argv[5])         # Number of independent simulations
dt = 1.0*skip*0.001				# value of 1 time frame(in pico-second [ps]); 
                                # Simulation dt = 1.0 femto-second
start_time = time.time()

def main(k):
    in_fname = ' ../../../r' + str(int(k+1)) + '_top.dat'   # Top wall LAMMPS data
	cuu_out_fname = 'r' + str(int(k+1)) + '_cuu.dat'
	cuf_out_fname = 'r' + str(int(k+1)) + '_cuf.dat'
	file1 = np.loadtxt(fname=in_fname, skiprows=2)
	
	num = file1[startTime:endTime:skip,4]
	uslab_x = file1[startTime:endTime:skip,5]*1000  # Convert A/fs to A/ps
	force_x = file1[startTime:endTime:skip,1]		# Select the x-comp of force [Kcal/(mol-A)]
	
	cuu, cuu_norm = Analysis.correlate(nrun, uslab_x, uslab_x)  # <u(0).u(t)>
	cuf, cuf_norm = Analysis.correlate(nrun, uslab_x, force_x)	# <u(0).F(t)>

	np.savetxt(cuu_out_fname, np.c_[np.arange(len(cuu))*dt, cuu, cuu_norm], fmt='%.10f \t %.20f  \t %.20f')  
	np.savetxt(cuf_out_fname, np.c_[np.arange(len(cuf))*dt, cuf, cuf_norm], fmt='%.10f \t %.20f \t %.20f') 
	return cuu, cuf

# ########### Parallelizing ########## #
with concurrent.futures.ProcessPoolExecutor() as executor:
	run_ids = range(irun)
	results = executor.map(main, run_ids)

end_time = time.time()
print('Finished in ' + str(round(end_time-start_time)) + ' seconds')
