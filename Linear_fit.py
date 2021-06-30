
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def mainfit(x, y, slope_guess, intercept_guess, include_intercept):

    def intercept_fit(x,m,c):
        return m*x + c
    def no_intercept_fit(x,m):
        return m*x

    if include_intercept:
        popt, pcov = curve_fit(intercept_fit, x, y, p0 =(slope_guess, intercept_guess))
        x_out = np.arange(0, max(x), 0.01)
        y_out = intercept_fit(x_out, popt[0], popt[1])
    else:
        popt, pcov = curve_fit(no_intercept_fit, x, y, p0 =(slope_guess))
        x_out = np.arange(0, max(x), 0.01)
        y_out = no_intercept_fit(x_out, popt)

    return x_out, y_out, popt, pcov
