import numpy as np
import sys
from scipy.optimize import curve_fit

sys.path.append(".")
from num_integration import *

class friction_coefficient:
    s_0 = 0
    s_end = 0
    s_num = 0
    
    @staticmethod
    def method1(Lcuu, Lcuf):
        """
        Function routine to compute the Navier friction coefficient
        according to METHOD-1
        """
        memory_kernel = -Lcuf/Lcuu
        return memory_kernel[0]
    
    @staticmethod
    def method2(frequ, Lcuu, Lcuf):
        """
        Function routine to compute the Navier friction coefficient
        according to METHOD-2
        """
        def func(X,B,L):
            s, jj = X
            return -(B*jj)/(s+L)
        
        popt,pcov = curve_fit(func, (frequ,Lcuu), Lcuf, maxfev=100000)
        zeta0 = popt[0]/popt[1]
        return zeta0
    
    @staticmethod
    def method3(s, Lcuu, Lcuf):
        """
        Function routine to compute the Navier friction coefficient
        according to METHOD-3
        """
        s_mid = int(len(s)/2)
        s_last = int(len(s))
        ds = s[1] - s[0]
        Sigma1 = s[0]
        Sigma2 = s[s_mid]
        Sigma3 = s[s_mid]
        Sigma4 = s[s_last-1]
        A = -1*np.divide(Lcuu, Lcuf)
        A1 = Integrate.simpson(A[0:s_mid], ds)
        A2 = Integrate.simpson(A[s_mid:s_last],ds)
        C = A1/A2
        lamda1 = ((Sigma2**2) - (Sigma1**2) + (C*((Sigma3**2)-(Sigma4**2))))/(2*((C*(Sigma4-Sigma3)) - Sigma2 + Sigma1))
        B1 = ((Sigma2**2) + (2*lamda1*Sigma2) - (Sigma1**2) - (2*lamda1*Sigma1))/(2*A1)
        zeta0 = B1/lamda1
        return zeta0
