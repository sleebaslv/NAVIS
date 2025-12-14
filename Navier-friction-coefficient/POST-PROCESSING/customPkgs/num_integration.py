import numpy as np

class Integrate:
	@staticmethod
	def trapz(func_x, delta_x):
        """
        Function routine for performing numerical intgration 
        according to Trapezoidal rule
        """
	        n = np.size(func_x)
	        f_x0 = func_x[0]
	        f_xn = func_x[n-1]
	        f_xmid = np.sum(func_x[1:(n-1)])
	        def_integral = (delta_x/2)*(f_x0 + 2*f_xmid + f_xn)
	        return def_integral
	
	@staticmethod
	def simpson(func_x, delta_x):
        """
        Function routine for performing numerical intgration 
        according to Simpson's rule
        """
	        n = np.size(func_x)
	        if n%2 != 0:
	                f_x0 = func_x[0]
	                f_xn = func_x[n-1]
	                f_xmidOdd = np.sum(func_x[1:(n-1):2])
	                f_xmidEven = np.sum(func_x[2:(n-2):2])
	                def_integral = (delta_x/3)*(f_x0 + 2*f_xmidEven + 4*f_xmidOdd + f_xn)
	                return def_integral
	        else:
	                raise Exception ('The number of data points for Simpson integration is even')
