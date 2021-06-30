import numpy as np
from skimage import io
import os
import scipy.ndimage as ndi
import scipy.signal as sig
from skimage import measure as sk
import matplotlib.pyplot as plt
import glob
#import radialAvg
import datetime
import multiprocessing

def radialProfile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr  # does this divide by number of pixels? Then it is truly average
    return radialprofile

def valueToIndex(list, value):
    index = np.abs(np.array(list) - value).argmin()
    return index

def findCenter(img, xSlice, ySlice, blur, verbose): #original Slice values
#def findCenter(img, xSlice=0, ySlice=0, blur=30, verbose=True):
    # Apply median filter to get rid of noisy peaks
    gaySauce = img
    #img = ndi.gaussian_filter(img, sigma=blur)
    img = np.uint32(ndi.median_filter(img, blur))

    # Find image dimensions
    xSize = img.shape[1] - 1
    ySize = img.shape[0] - 1

    # Define slice start points and end points
    x1, y1 = [0, xSize], [ySize / 2 + ySlice, ySize / 2 + ySlice]
    x2, y2 = [xSize / 2 + xSlice, xSize / 2 + xSlice], [0, ySize]

    # Take profile of image along slice
    xProfile = sk.profile_line(img, (y1[0], x1[0]), (y1[1], x1[1]), mode='reflect')
    yProfile = sk.profile_line(img, (y2[0], x2[0]), (y2[1], x2[1]), mode='reflect')
    # Assign xCenter and yCenter to position of maximum of profile
    xCenter = valueToIndex(xProfile, max(xProfile[500:680]))
    yCenter = valueToIndex(yProfile, max(yProfile[450:615]))
    print('Center(x,y):')    
    print(xCenter, yCenter)

    # Plot stuff if verbose is True
    if verbose:
        # Get max value for scaling.

        plt.imshow(gaySauce, extent=[0, xSize, ySize, 0])
        plt.clim(0, 2000)

        plt.plot(x1, y1)
        plt.plot(x2, y2)
        plt.plot(xCenter, yCenter, marker='.')

        plt.gcf().gca().add_artist(plt.Circle((xCenter, yCenter), 150, fill=False))

        plt.plot(xProfile/max(xProfile)*500)
        plt.plot(yProfile/max(yProfile)*500)
        plt.show()

    return (xCenter, yCenter)

def cropImage(img, xCenter, yCenter, cropWidth):
    return img[yCenter - cropWidth:yCenter + cropWidth, xCenter - cropWidth:xCenter + cropWidth]

def meanAroundRing(img, initCenter, step, radiusAvgMin, radiusAvgMax, verbose):
    # Find initial coordinate guess.
    # Crop image arbitrarily around center of image.
    xSize = img.shape[1] - 1
    ySize = img.shape[0] - 1
    # croppedImg = cropImage(img, int(xSize/2), int(ySize/2), 400)

#    step = 8
    xRange = list(range(initCenter[0] - step, initCenter[0] + step))
    yRange = list(range(initCenter[1] - step, initCenter[1] + step))
    meanValList = []
    xList = []
    yList = []

    for x in xRange:
        for y in yRange:
            # Take radial average of point.
            # print(img.dtype)
            raTemp = radialProfile(img, (x,y))

            # Take mean between points.
            meanVal = np.mean(raTemp[radiusAvgMin:radiusAvgMax])
            # meanVal = np.mean()
            meanValList.append(meanVal)
            xList.append(x)
            yList.append(y)
#    print('datetime.now:')
#    print(datetime.datetime.now() - t0)

    # Return maximum of mean.
    xMean = xList[valueToIndex(meanValList, max(meanValList))]
    yMean = yList[valueToIndex(meanValList, max(meanValList))]

    # Print the two above results.
    print("Mean Center:")
    print(str(xMean)+','+str(yMean))

    if verbose:
        plt.figure('Mean Around Ring')
        plt.imshow(img, extent=[0, xSize, ySize, 0])
        plt.clim(700, 5000)
        plt.plot(xMean, yMean, marker='.')
        
        modify_colorlimits=0
        if modify_colorlimits:
            intensity_range = np.amax(img) - np.amin(img)
            intensity_scale = 0.3
            clim_mod = intensity_range* intensity_scale
            plt.clim(np.amin(img), np.amax(img)-clim_mod )
        
        plt.colorbar()
        plt.gcf().gca().add_artist(plt.Circle((xMean, yMean), radiusAvgMin, fill=False))
        plt.gcf().gca().add_artist(plt.Circle((xMean, yMean), radiusAvgMax, fill=False))
#        plt.gcf().gca().add_artist(plt.Circle((xMean, yMean), radiusAvgMax + radius_step, fill=False))
#        plt.gcf().gca().add_artist(plt.Circle((xMean, yMean), radiusAvgMax + 2*radius_step, fill=False))
        # plot crosshairs
        plt.axvline(x=xMean, color='r', linestyle='--')
        plt.axhline(y=yMean, color='r', linestyle='--')
        plt.title('meanAroundRing')
        plt.show()

    return xMean, yMean

    plt.show()
    
#showscript()
