
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def mainfitting(timedata, peakdata, a_fit, c_fit, tau_fit, t0_fit, plotlabel):
    times = np.array(timedata)   # convert whatever into array
    peaks = np.array(peakdata)

    def piecewise_exp_function(t, a, c, tau, t0):
        if t < t0:
            return a + c
        else:
            return a*np.exp(-1.0*(t - t0)/tau) + c

    fitfunc_vec = np.vectorize(piecewise_exp_function)
    def fitfunc_vec_self(t,a,c,tau,t0):
        y = np.zeros(t.shape)
        for i in range(len(y)):
            y[i]=piecewise_exp_function(t[i],a,c,tau,t0)
        return y

    popt, pcov = curve_fit(fitfunc_vec_self, times, peaks, p0 =(a_fit,c_fit,tau_fit,t0_fit))
    # print("-----" + plotlabel)
    # print("fitted a= " + str(round(popt[0],3)))
    # print("fitted c= " + str(round(popt[1],3)))
    # print("fitted tau= " + str(round(popt[2],3)))
    # print("fitted t0= " + str(round(popt[3],3)))

    x_out = np.arange(min(times),max(times), 0.01)
    y_out = fitfunc_vec_self(x_out, popt[0],popt[1],popt[2],popt[3])

    return x_out,y_out,popt,pcov
