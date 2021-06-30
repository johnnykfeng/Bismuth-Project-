
# -*- coding: utf-8 -*-
"""
Bismuth_FindT0_v2.py
"""
import os
import glob
from skimage import io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.signal import find_peaks
from scipy import interpolate
from matplotlib.font_manager import FontProperties
import findCenters
import exponential_fit as exp_fit
import Linear_fit
from astropy.table import QTable, Table, Column
from astropy import units as u
from ScanDictionary import scanlist
import time
import csv
from BaseSubtraction import BaseSubtraction



#region MAIN SCAN LOOP
scan_index_length = 1
scan_start = 1
fit_var = []

for scan_index, scan_number in enumerate(np.arange(scan_start, scan_start + scan_index_length, 1)):
    print('+-+-+-+-+-+-+-+- scan index: '+ str(scan_index))
    print('+-+-+-+-+-+-+-+- scan number: '+ str(scan_number))

    #region Data Directory and dictionary structure
    # scan_number = 11
    dataDirec = 'D:\\Bismuth Project\\New Bismuth Data\\' + scanlist[scan_number]['date'] + '\\scans\\scan' + scanlist[scan_number]['scan']
    print('Data Directory: ' + dataDirec)
    plot_title = 'Scan #' +scanlist[scan_number]['scanindex'] + '  ' + \
                 'thickness= ' + scanlist[scan_number]['thickness'] + ' $nm$,  ' + \
                 'fluence= ' + scanlist[scan_number]['fluence'] + ' $mJ/cm^2$,  ' + \
                 'exposure= ' + scanlist[scan_number]['exposure'] + ' $s$'
    plot_title_short = 'thickness= ' + scanlist[scan_number]['thickness'] + ' $nm$,  ' + \
                       'fluence= ' + scanlist[scan_number]['fluence'] + ' $mJ/cm^2$'
    #endregion

    #region Toggle Plotting Variables
    #region pixel2Ghkl functions
    pixel2Ghkl = 2.126e-3   # converts CCD pixels to G_hkl = 1/d_hkl 'scattering vector'
    def pixel2scatt_vector(x):
        return x*pixel2Ghkl
    def scatt_vector2pixel(x):
        return x/pixel2Ghkl
    #endregion
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    show_centerfinder = 0
    show_all_diffImg =0  # toggles diff_Img plot
    show_vertical_lines = 1  # toggles vertical integration lines in diff_Img plot
    show_on_img = 0
    show_off_img = 0
    show_liquid_rise = 0
    show_expfits = 0
    normalize_intensity = 1  # normalize intensity I(t)/I_t0 in expfits
    show_flattened_on_img = 1
    plot_find_peaks = 0
    debye_waller_analysis = 0
    write_csv_file = 0

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #endregion

    #region Initializing variables for data collection
    onList = []  # for collecting ON images filenames
    offList = []  # for collecting OFF images filenames
    diffList = []  # for collecting DIFF images
    darkBg = np.zeros((1000, 1148))
    scatterBg = np.zeros((1000, 1148))
    onImg = np.zeros((1000, 1148))
    offImg = np.zeros((1000, 1148))
    sum_all_img = np.zeros((1000, 1148))
    t_list = []  # converts timepoint data to ps
    # for plotting the liquid rise
    lr_range = (148, 190)
    lr_data = []
    dby = []   # collecting y data for debye-waller analysis
    dbx = []   # collecting x data for debye-waller analysis
    y_peaks = []
    #endregion

    #region Image File Extraction - Computing On, Off, Diff Img

    #region Loading data images from Directory
    for file in glob.glob(os.path.join(dataDirec, '*.tiff')):
        # Iterate through file, find the ones we want.
        if (os.path.basename(file).split('_')[5] == 'On'):
            onList.append(file)
        if (os.path.basename(file).split('_')[5] == 'Off'):
            offList.append(file)
    print('loading on and off images successful')
    try:
        darkBg = io.imread(os.path.join(dataDirec, 'darkExp.TIF'))
        scatterBg = io.imread(os.path.join(dataDirec, 'pumpExp.TIF'))
        print('reading darkBg and scatterBg successful')
    except:
        print('dark images not okay!')
    print('Number of data points: ' + str(len(onList)))
    #endregion
    for j in range(len(onList)):
        # timepoint = os.path.basename(onList[j]).split('_')[3]  # reads time delay from file name
        # imgnum = os.path.basename(onList[j]).split('_')[1]  # not really used anymore
        onImg = io.imread(onList[j]) - scatterBg
        # onImg = io.imread(onList[j])
        offImg = io.imread(offList[j]) - darkBg
        # offImg = io.imread(offList[j])
        diffImg = (onImg - offImg)
        diffList.append(diffImg)
        sum_all_img = sum_all_img + offImg  # this variable sometimes isn't used
    #endregion

    #region FindCenter control
    start_timer = time.time()
    initCenter = (598, 543)
    step = 8;     radiusAvgMin = 140;     radiusAvgMax = 150
    center0 = findCenters.meanAroundRing(sum_all_img, initCenter, step, radiusAvgMin, radiusAvgMax, show_centerfinder)
    find_center_runtime = time.time() - start_timer
    print('Find Center Runtime = ' + str(find_center_runtime))
    #endregion

    #region Control plotting and integration variables
    plotting_index = np.arange(1, len(onList))  # skip certain timepoints in the scan, such as the first
    # peak_choices = (1, 2, 5, 6)
    peak_choices = (1, 2,3, 4,5,6)
    peak_step = 6
    peak_sum = np.zeros((len(peak_choices), len(plotting_index), scan_index_length))  # empty array for collecting the integrated peak data
    #endregion

    #region Plotting colors and fonts
    colors = cm.jet(np.linspace(0, 1, len(onList)))
    peak_colors = cm.turbo(np.linspace(0, 1, len(peak_choices)))
    colors_fluence = cm.Set1(np.linspace(0, 1, scan_index_length))
    fontP = FontProperties()
    fontP.set_size('x-small')
    #endregion

    #region FIND PEAKS_1  (not used anymore)
    FIND_PEAKS = 0
    start_timer = time.time()
    if FIND_PEAKS:
        image_index = 1  # first ON image, usually before t0
        timepoint = os.path.basename(onList[image_index]).split('_')[3]
        rad_avg_peakimg = findCenters.radialProfile(io.imread(onList[image_index]), center0)
        #----------- MAIN PEAK FINDING LINE OF CODE ------------ #
        peakposn, properties = find_peaks(rad_avg_peakimg, threshold=None, distance=20, prominence=1, width=2)
        bases = properties['left_bases']  # finds the bases for each peak
        bases = np.append(bases, [375,420])
        f = interpolate.interp1d(bases, rad_avg_peakimg[bases], kind='cubic')

        x_bases = np.arange(min(bases), max(bases), 1)
        spline_bg = f(x_bases)   # spline of bases, to be subtracted as background
        flattened_rad_avg = rad_avg_peakimg[min(bases): max(bases)]- spline_bg
        flattened_rad_avg = flattened_rad_avg -min(flattened_rad_avg[100:])
        peakposition, properties = find_peaks(flattened_rad_avg, threshold=None, distance=20, prominence=5, width=5)

        # plot_find_peaks = 1
        if plot_find_peaks:
            fig, ax = plt.subplots()
            ax.plot(rad_avg_peakimg, color=colors[-image_index], linewidth=1, linestyle='--', label= 'raw rad avg')
            ax.plot(peakposn, rad_avg_peakimg[peakposn], "x", markersize=5, color= 'blue')

            ax.plot(bases, rad_avg_peakimg[bases], 'o', color='red')   # show bases as red circle
            ax.plot(spline_bg, label= 'fitted background')

            # --- FIND PEAKS FOR THE FLATTENED RADIAL AVERAGE ----- #
            ax.plot(flattened_rad_avg, color=colors[image_index], linewidth=2, linestyle='-', label= 'flattened rad avg ')
            ax.plot(peakposition, flattened_rad_avg[peakposition], "x", markersize=5, color= 'green')
            # LABEL peaks in the plots
            for t in range(len(peakposition)):
                peakposn_text = str(peakposition[t])+ ', ' + str(np.round(pixel2Ghkl*peakposition[t],  3)) + '$A^{-1}$ '
                plt.text(peakposition[t] - 10, flattened_rad_avg[peakposition][t] + 300, peakposn_text)  # labels the peak positions
                plt.text(peakposition[t], flattened_rad_avg[peakposition][t] - 50, str(t), fontweight='bold') # labels the peak index number

            ax.set_xlabel('CCD pixels')
            secax = ax.secondary_xaxis('top', functions = (pixel2scatt_vector, scatt_vector2pixel))
            secax.set_xlabel('scattering vector (1/A)')
            plt.legend()
            plt.grid(True)
            plt.title(plot_title)
    find_peaks_runtime = time.time() - start_timer
    print('Find Peaks Runtime = ' + str(find_peaks_runtime))
    #endregion_1 #

    #region New Peak Finding using BaseSubtraction.py
    image_index = 1  # first ON image, usually before t0
    rad_avg_img = findCenters.radialProfile(io.imread(onList[image_index]), center0) # image
    flattened_rad_avg, peakposition, properties, bases1, bases2, bases3, spline_bg1, spline_bg2 = BaseSubtraction(rad_avg_img, 1, 0)
    plt.figure('Demonstrate peak and base finding, flattened_rad_avg')
    plt.plot(flattened_rad_avg)
    pp = peakposition[:-2]
    plt.plot(pp, flattened_rad_avg[pp], 'x')
    plt.plot(bases3, flattened_rad_avg[bases3], 'o')
    plt.grid(True)
    #endregion

    #region Main Loop for Computing and Plotting Radial Average
    start_timer = time.time()
    for i in plotting_index:
        # print(i)
        timepoint = os.path.basename(onList[i]).split('_')[3]
        imgnum = os.path.basename(onList[i]).split('_')[1]
        t_ps = int(timepoint) * 1e-3
        t_list.append(t_ps)  # accumulates the timepoint as ps
        # center0 = findCenters.meanAroundRing(io.imread(offList[i]), initCenter, 3, radiusAvgMin, radiusAvgMax, 0)
        rad_avg_diff = findCenters.radialProfile(diffList[i], center0)
        rad_avg_on = findCenters.radialProfile(io.imread(onList[i]), center0)
        rad_avg_off = findCenters.radialProfile(io.imread(offList[i]), center0)
        flattened_rad_avg_on = rad_avg_on[min(bases1):max(bases1)] - spline_bg1
        flattened_rad_avg_on = flattened_rad_avg_on[min(bases2):max(bases2)] - spline_bg2
        lr_data.append(sum(rad_avg_diff[lr_range[0]:lr_range[1]]))  # integrates pixels for liquid rise analysis

        if show_all_diffImg:
            plt.figure('[diffImg] ' + 'fluence:'+ scanlist[scan_number]['fluence'] + ', thickness:' + scanlist[scan_number]['thickness'] )
            if i % 2 == 1:
                # plt.plot(radial_distance*radius2dhkl, rad_avg_diff, color=colors[i], linewidth=2, linestyle='-', label=str(t_ps) + 'ps ')
                plt.plot(rad_avg_diff, color=colors[i], linewidth=2, linestyle='-', label=str(t_ps) + 'ps ')
            else:
                # plt.plot(radial_distance*radius2dhkl, rad_avg_diff, color=colors[i], linewidth=2, linestyle='--', label=str(t_ps) + 'ps ')
                plt.plot(rad_avg_diff, color=colors[i], linewidth=2, linestyle='--', label=str(t_ps) + 'ps ')
            if show_vertical_lines:
                for p in peak_choices:
                    plt.axvline(x=(peakposition[p]-peak_step), linestyle='--')
                    plt.axvline(x=(peakposition[p]+peak_step), linestyle='--')
            # if show_vertical_lines:
            #     for p in peak_choices:
            #         plt.axvline(x_bases=(peakposition[p]-peak_step), linestyle='--')
            #         plt.axvline(x_bases=(peakposition[p]+peak_step), linestyle='--')

            plt.legend(title='timepoints', bbox_to_anchor=(1, 1), loc='upper left', prop=fontP)
            plt.xlim(125, 600)
            plt.title(plot_title)
            plt.xlabel('Radius from center (pixel)')
            plt.ylabel('Radial average difference intensity')
            plt.grid(True)

        if show_off_img:
            plt.figure('OFF - ' + 'fluence:'+ scanlist[scan_number]['fluence'] + ', thickness:' + scanlist[scan_number]['thickness'] )
            if i % 2 == 1:
                plt.plot(rad_avg_off, color=colors[i], linewidth=2, linestyle='-', label=str(t_ps) + 'ps ')
            else:
                plt.plot(rad_avg_off, color=colors[i], linewidth=2, linestyle='--', label=str(t_ps) + 'ps ')

            if show_vertical_lines:
                for p in peak_choices:
                    plt.axvline(x=(peakposition[p]-peak_step), linestyle='--')
                    plt.axvline(x=(peakposition[p]+peak_step), linestyle='--')

            plt.legend(title='timepoints', bbox_to_anchor=(1, 1), loc='upper left', prop=fontP)
            plt.xlim(125, 400)
            plt.ylim(0, 12500)
            plt.title(plot_title)
            plt.xlabel('Radius from center (pixel)')
            plt.ylabel('Radial average difference intensity')
            plt.grid(True)  #

        if show_on_img:
            plt.figure('ON Images')
            if i % 2 == 1:
                plt.plot(rad_avg_on , color=colors[i], linewidth=2, linestyle='-', label=str(t_ps) + 'ps ')
            else:
                plt.plot(rad_avg_on, color=colors[i], linewidth=2, linestyle='--', label=str(t_ps) + 'ps ')

            plt.legend(title='timepoints', bbox_to_anchor=(1, 1), loc='upper left', prop=fontP)
            plt.xlim(min(bases), max(bases))

            plt.title(plot_title)
            plt.xlabel('Radius from center (pixel)')
            plt.ylabel('Radial average difference intensity')
            plt.grid(True)

        # flattened_rad_avg_on = rad_avg_on[min(bases):max(bases)] - spline_bg
        # show_flattened_on_img = 0
        if show_flattened_on_img:
            plt.figure('fluence:'+ scanlist[scan_number]['fluence'] + ', thickness:' + scanlist[scan_number]['thickness'] )

            if i % 2 == 1:
                # axs[scan_index].plot(flattened_rad_avg_on , color=colors[i], linewidth=2, linestyle='-', label=str(t_ps) + 'ps ')
                plt.plot(flattened_rad_avg_on , color=colors[i], linewidth=2, linestyle='-', label=str(t_ps) + 'ps ')
            else:
                # axs[scan_index].plot(flattened_rad_avg_on , color=colors[i], linewidth=2, linestyle='--', label=str(t_ps) + 'ps ')
                plt.plot(flattened_rad_avg_on, color=colors[i], linewidth=2, linestyle='--', label=str(t_ps) + 'ps ')

            if show_vertical_lines:
                for p in peak_choices:
                    plt.axvline(x=(peakposition[p]-peak_step), linestyle='--')
                    plt.axvline(x=(peakposition[p]+peak_step), linestyle='--')
                # LABEL peaks in the plots
            for t in range(len(peakposition)):
                peakposn_text = str(peakposition[t])
                plt.text(peakposition[t], flattened_rad_avg[peakposition][t] - 50, str(t),
                         fontweight='bold')  # labels the peak index number

            plt.legend(title='timepoints', bbox_to_anchor=(1, 1), loc='upper left', prop=fontP)
            plt.xlim(100, 400)
            plt.ylim(-500, 8500)
            plt.title(plot_title)
            plt.xlabel('Radius from center (pixel)')
            plt.ylabel('Radial average difference intensity')
            plt.grid(True)

        # Integrates the peaks in the radial average and accumulates them in the 'peak_sum' array
        for p_index, p in enumerate(peak_choices):
            peak_sum[p_index, i-1, scan_index] = sum(flattened_rad_avg_on[peakposition[p]-peak_step: peakposition[p]+peak_step])
    radialprofile_runtime = time.time() - start_timer
    print('radialprofile_runtime = ' + str(radialprofile_runtime))
    #endregion

    if show_expfits:
        fitted_variables = np.zeros((len(peak_choices), 5))
        fig, ax1 = plt.subplots()
        skipfirst = 0
        skiplast = 1
        x_peaks = t_list[skipfirst: -skiplast] #picks timepoints for the exponential fitting
        # Does the exponential fit for each peak

        for p_index, p in enumerate(peak_choices):

            a, c, tau, t0 = 8000, -8000, 2.0, 0
            if p == 3 or p == 4:
                a, c, tau, t0 = -5000, 5000, 2.0, 0

            y_peaks = (peak_sum[p_index, skipfirst:-skiplast, scan_index])
            x_fit, y_fit, popt, pcov = exp_fit.mainfitting(x_peaks, y_peaks, a, c, tau, t0, 'Peak '+ str(p))
            fitted_variables[p_index, :] = p, popt[0], popt[1], popt[2], popt[3]  # stick the fitted variables here to make a table

            if normalize_intensity:
                y_fit = np.true_divide(y_fit, popt[0]+popt[1])
                y_peaks = np.true_divide(y_peaks, popt[0]+popt[1])

            t_zero_correction = True
            if t_zero_correction:
                x_peaks = x_peaks - popt[3]
                x_fit = x_fit - popt[3]

            dby.append(-np.log( sum(y_peaks[-5:])/5.0 )) # ln of average of last 5 data points in the trace
            dbx.append((peakposition[p]*pixel2Ghkl)**2) # square of the scattering vector

            ax1.plot(x_peaks, y_peaks, '-o', color=peak_colors[p_index],  label = "Peak #" + str(p))
            ax1.plot(x_fit, y_fit, '--', color=peak_colors[p_index])

            ax1.tick_params(axis='y')
            ax1.set_ylabel('Normalized peak intensity')
            ax1.set_xlabel('Timepoint (ps)')
            ax1.legend(loc='center left')
            ax1.grid(True)
            ax1.set_title(plot_title)

            fluence_label = scanlist[scan_number]['fluence'] + ' $mJ/cm^2$'
            peakfigure = plt.figure('peak # ' + str(p))
            # peakfigure.set_title('Time Trace Peak # ' + str(p))
            plt.plot(x_peaks, y_peaks, '-o', color = colors_fluence[scan_index], label = fluence_label)
            plt.plot(x_fit, y_fit, '--',color = colors_fluence[scan_index])
            plt.xlabel('Timepoint (ps)')
            plt.ylabel('Normalized peak intensity')
            plt.grid(True)
            plt.title('peak # ' + str(p))
            plt.legend()

        ravel_fitted_variables = np.ravel(fitted_variables)
        # ravel_fitted_variables.insert(0, float(scanlist[scan_number]['fluence']))
        ravel_fitted_variables = np.insert(ravel_fitted_variables, 0, float(scanlist[scan_number]['fluence']))
        fit_var.append(ravel_fitted_variables)

        if show_liquid_rise:
            a, c, tau, t0 = -20000, 0, 0.5, 0
            y_peaks = lr_data[skipfirst:-skiplast]
            x_fit, y_fit, popt, pcov =exp_fit.mainfitting(x_peaks, y_peaks, a, c, tau, t0, 'Liquid Rise')
            ax2 = ax1.twinx()
            liquid_color = 'cadetblue'
            if normalize_intensity:
                y_fit = np.true_divide(y_fit, popt[0] + popt[1])
                y_peaks = np.true_divide(y_peaks, popt[0] + popt[1])
            ax2.plot(x_peaks, y_peaks, '-D', label="Liquid rise",color =liquid_color)
            ax2.plot(x_fit,y_fit,'--', color =liquid_color)
            ax2.set_ylabel('Liquid Rise scale', color =liquid_color)
            ax2.tick_params(axis='y_peaks', labelcolor=liquid_color)
            ax2.legend(bbox_to_anchor=(0, 0.25), loc='center left')

        output_qtable=1
        if output_qtable:
            qtable = QTable()
            qtable['Peak #'] = np.round(fitted_variables[:,0], 0)
            qtable['Tau'] = np.round(fitted_variables[:,3] , 2) *u.ps
            qtable['t0'] = np.round(fitted_variables[:,4] , 2)*u.ps
            qtable['a'] = np.round(fitted_variables[:, 1], 2)
            qtable['c'] = np.round(fitted_variables[:, 2], 2)
            print(qtable)

#region DEBYE WALLER ANALYSIS
    # debye_waller_analysis = 1
    if debye_waller_analysis:

        dbx_fit, dby_fit, popt, pcov = Linear_fit.mainfit(dbx, dby, 0.35, 0, True)
        dbx_fit2, dby_fit2, popt2, pcov2 = Linear_fit.mainfit(dbx, dby, 0.4, 0, False)
        print('slope_fit: ' + str(np.round(popt[0], 3)) + ', intercept_fit: ' + str(np.round(popt[1], 3))  )
        print('slope_fit no-intercept: ' + str(np.round(popt2, 3)) )

        plt.figure('Debye Waller')
        fluence_label = scanlist[scan_number]['fluence'] + ' $mJ/cm^2$'
        # plt.plot(dbx[:, scan_index], dby[:, scan_index], '-o', color = colors_fluence[scan_index], label = fluence_label)
        plt.plot(dbx, dby, '-o', color = colors_fluence[scan_index], label = fluence_label)
        # plt.plot(dbx_fit, dby_fit, color = colors_fluence[scan_number], linestyle = '--')  # plotting the fit line
        # plt.plot(dbx_fit2, dby_fit2, color = colors_fluence[scan_number], linestyle = '-.')  # plotting the fit line

        plt.ylabel('$-ln( I(t)/I_{t0} )$')
        plt.xlabel('$G_{hkl}^2$  ' +' ($A^{-2}$) ')
        plt.grid(True)
        plt.legend()
        plt.title('Debye waller analysis of peaks')
    #endregion

#region WRITE CSV FILE
if write_csv_file:
    print(fit_var)
    with open('fit_var_14nm.csv', 'w', newline='') as file:
        mywriter = csv.writer(file, delimiter = ',')
        mywriter.writerows(fit_var)
#endregion

plt.show()   # last line of the code

#endregion
