import sys
import numpy as np

sys.path.append(".")
from num_integration import *   # For performing numerical intgration

class Analysis:
    s_0 = 0
    s_end = 0
    s_num = 0
    
    @staticmethod
    def correlate(no_of_sets, data1, data2):
        """
        Function routine to compute the time correlation between 
        two datasets.
        """
        no_of_cols = int(len(data1)/no_of_sets)
        data1 = np.array_split(data1, no_of_sets)
        data2 = np.array_split(data2, no_of_sets)
        corr_array = np.zeros((no_of_sets, no_of_cols))
        for j in range(no_of_sets):
            a = data1[j]
            b = data2[j]
            N1 = len(data1[j])
            corr = []
            for k in range(N1):
                Cab = np.sum(a[0:(N1-k)]*b[k:N1])/(N1-k)
                corr = np.append(corr,Cab)
            corr_array[j,:] = corr
        corravg = np.mean(corr_array, axis=0)
        corravg_norm = corravg/np.max(np.absolute(corravg))
        return corravg, corravg_norm

# ################# Laplace Transform ############## #
    @staticmethod
    def laplace(time, f_t, s_0, s_end, s_num):
        """
        Function routine to compute the Laplace transform of a dataset
        """
        dt = time[1] - time[0]
        s = np.linspace(s_0, s_end, num=s_num)
        Analysis.s_0, Analysis.s_end, Analysis.s_num = s_0, s_end, s_num
        Lfs = np.zeros(np.size(s))
        func_t = []
        for j in range(np.size(s)):
            e_st = np.exp(-s[j]*time)
            func_t = np.multiply(f_t, e_st)
            Lfs[j] = Integrate.simpson(func_t, dt)  # Simpson's rule       
        return s, Lfs
