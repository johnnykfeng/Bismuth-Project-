import os
import glob
from skimage import io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.signal import find_peaks
from scipy.signal import argrelmin
from scipy import interpolate
from matplotlib.font_manager import FontProperties


def removeduplicates(sample_array):
    temp_list = []
    for i in sample_array:
        if i not in temp_list:
            temp_list.append(i)
    return temp_list

def BaseSubtraction(rad_avg1, do_second_run, plot_find_peaks):

    peakposn1, properties = find_peaks(rad_avg1, threshold=None, distance=20, prominence=1, width=2)  # finds peaks and properties
    bases1 = properties['left_bases']  # finds the bases for each peak
    bases1 = np.append(bases1, [375, 420]) # adds some extra points to define the long tail
    f = interpolate.interp1d(bases1, rad_avg1[bases1], kind='cubic')  # creates the spline function that goes through the base points
    x_bases1 = np.arange(min(bases1), max(bases1), 1) #
    spline_bg1 = f(x_bases1)  # spline of bases, to be subtracted as background
    rad_avg2 = rad_avg1[min(bases1): max(bases1)] - spline_bg1   # create new curve that subtracts the generated background

    if do_second_run:
        peakposn2, properties = find_peaks(rad_avg2, threshold=None, distance=20, prominence=10, width=5, height=20)  # repeat the process
        # bases2 = properties['left_bases']  # finds the bases for each peak
        bases2 =  argrelmin(rad_avg2, order=10)[0]  # way better way of finding bases
        f2 = interpolate.interp1d(bases2, rad_avg2[bases2], kind= "cubic")
        x_bases2 = np.arange(min(bases2), max(bases2), 1)
        spline_bg2 = f2(x_bases2)
        rad_avg3 = rad_avg2[min(bases2): max(bases2)] - spline_bg2

    # final peakposition and bases using find_peaks
    peakposition3, properties3 = find_peaks(rad_avg3, threshold=None, distance=20, prominence=2, width=2)  # repeat the process
    bases3 = argrelmin(rad_avg3, order=20)[0]

    print_peaksandbases = 0
    if print_peaksandbases:
        print('~peakposn1:')
        print(peakposn1)
        print('~bases1:')
        print(bases1)
        print('~peakposn2:')
        print(peakposn2)
        print('~bases2:')
        print(bases2)
        print('~peakposition3:')
        print(peakposition3 + min(x_bases2))
        print('~bases3:')
        print(bases3 + min(x_bases2))

    plt.figure('BaseSubtraction code')

    # --- plot the curves --- #
    plt.plot(rad_avg1, label = 'original' , color = 'k')
    plt.plot(x_bases1 , rad_avg2, label = 'rad_avg2', color = 'r')
    plt.plot(x_bases2 + min(bases1), rad_avg3, label = 'rad_avg3', color = 'b')
    plt.plot(x_bases1, spline_bg1, '--')
    plt.plot(x_bases2 + min(bases1), spline_bg2, '--')

    # --- plot PEAKS with 'x' ---#
    plt.plot(peakposn2 + min(x_bases1), rad_avg2[peakposn2], 'x', color = 'r')
    plt.plot(peakposition3 + min(x_bases1) + min(bases2), rad_avg3[peakposition3], 'x', color = 'b')

    # --- plot BASES with 'o' ---#
    plt.plot(bases1, rad_avg1[bases1], 'o', color = 'k')
    plt.plot(bases2 + min(bases1), rad_avg2[bases2], '^', color = 'r')
    plt.plot(bases3 + min(bases1) + min(bases2), rad_avg3[bases3], '*', color = 'b')
    plt.grid(True)
    plt.legend()
    # plt.show()

    return rad_avg3, peakposition3, properties3, bases1, bases2, bases3, spline_bg1, spline_bg2





