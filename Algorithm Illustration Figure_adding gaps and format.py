# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 15:39:55 2022

Main: imports material properties and data from csv file, automatically calculates different smoothing parameters and associated loss function. Systematic comparison of different loss functions on smoothing.
Production of plots contrasting different smoothing paramters and loss functions.

Version: specific functions taken from Nico's code, Stuart reworked much of the structure and changed syntax to iterate through different loss functions as well as spectra and smoothing parameters
Current version iterates through a range of lowess smoothing windows, saving the outputs from each. The program then iterates through a range of gamma values in the loss function and selects the smoothed outputs according to the loss function. Design is to test possibility of a single, consistent gamma value in loss function.

@author: nicok
"""

import numpy as np
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
from tkinter.filedialog import askopenfilename
from matplotlib.widgets import Cursor
from scipy import interpolate
import scipy.optimize
from properties import *
import tkinter as tk
from tkinter import Tk
from GUI import material,methode
import timeit

AllPlots=False
Tk().withdraw()
Filename = askopenfilename()                                                #show an "Open" dialog box and return the path to the selected file      
algorithm="Nelder-Mead"                                                               #spektra excluded

short_name = Filename[Filename.rfind('/')+1:Filename.find('.csv')]          #extra the filename without path for saving files
print(short_name)
#%%

objekt=Material_Methode(material, methode)                                  #object with defined properties
try:
    data = np.genfromtxt(Filename,                                              #load data skip headlines
                 skip_header=2,delimiter=",")
except ValueError:
    data = np.genfromtxt(Filename,                                              #load data skip headlines
                 skip_header=2,delimiter=",",invalid_raise=False)
    print("Non-numerical values detected in csv file. Skipped parameter lines.")

wavelength=data[:,0] 
spectra=np.zeros([len(data[0]),len(data)])                                  #helparray for spectra data
count0=0
data_space = np.mean(wavelength[1:]-wavelength[:-1])

colour_list = ['#497AB9','violet','purple']

for i in range(len(data[0])):                                               #create spectra array out of data
    if sum(data[:,i])/len(data[:,i]) <= 200:     #NEW: since most spectrometers max at 10, redefined limit to that
        for j in range (len(data)):
            spectra[i][j]=data[j,i]
    if sum(spectra[i])==0.0:
        count0+=1
     
try:
    for i in range(len(spectra)-count0+1):                                      #delete empty lists in spectra-array
        if sum(spectra[i])== 0.0:
            spectra = np.delete(spectra,(i), axis = 0) 
except IndexError:
    for i in range(len(spectra)-count0):                                      #delete empty lists in spectra-array
        if sum(spectra[i])== 0.0:
            spectra = np.delete(spectra,(i), axis = 0)

labels = np.loadtxt(Filename,                                               #load header as str for labeling
                 delimiter=",",dtype=str,max_rows=1)

for i in range(len(labels)-count0):                                         #delete empty labels
    if labels[i]=="":
        labels = np.delete(labels,(i), axis = 0) 

#definition of left and right limit of analysis (for conventional measurement from high to low wavelength, left will be higher wavelengths). Note: smoothing and differentiation is applied to the entire spectrum to minimise boundary errors. Use of energy (function in properties) to correct to wavelength scale
lf_lim = ((objekt.index(objekt.wavelength_energy(objekt.E_bulk),wavelength)) - int(45/abs(data_space)))     #lower index, which means higher wavelength
rt_lim = ((objekt.index(objekt.wavelength_energy(objekt.E_ML),wavelength)) + int(45/abs(data_space)))     #higher index, which means lower wavelength
###Functions

#Smoothness parameter. Possibility of using the relative difference between adjacent points as per paper, or calculating the correlation coefficient of local points
def correlation_subset(x, y):
    r2_subset = []
    for p in np.arange(0,len(y)-subset_size,1):   
        #subset_x = x[p:p+subset_size]
        subset_y = y[p:p+subset_size]
        
        if np.std(subset_y) != 0:
            #r2_subset.append((np.corrcoef(subset_x,subset_y)[0][1])**2)
            r2_subset.append((np.std(subset_y)/np.mean(subset_y)))
        else:
            r2_subset.append(0.0)
    return r2_subset

#Root Mean Square Error of the smoothed and original data
def rmse(predicted,actual):
    summ = 0
    for p in range (len(predicted)):
        summ += (predicted[p]-actual[p])**2
    return np.sqrt(summ/len(predicted))

#Minimisation of the integral area either side of the central value. By minimising the difference, the center of mass of the peak can be calculated
def areamin (mini: int)-> float:
    arealeft=np.trapz(interpollist[interpolright:int(mini)],x_new[interpolright:int(mini)])
    arearight=np.trapz(interpollist[int(mini):interpolleft],x_new[int(mini):interpolleft])
    difference=abs(arealeft-arearight)
    return difference

### Atempt colour map
def hex_to_RGB(hex_str):
    """ #FFFFFF -> [255,255,255]"""
    #Pass 16 to the integer function for change of base
    return [int(hex_str[i:i+2], 16) for i in range(1,6,2)]

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#Colour gradient function can be used in plots to automatically form a pleasant blue-red gradient for a trend like gamma
def get_color_gradient(c1, c2, n):
    """
    Given two hex colors, returns a color gradient
    with n colors.
    """
    assert n > 1
    c1_rgb = np.array(hex_to_RGB(c1))/255
    c2_rgb = np.array(hex_to_RGB(c2))/255
    mix_pcts = [x/(n-1) for x in range(n)]
    rgb_colors = [((1-mix)*c1_rgb + (mix*c2_rgb)) for mix in mix_pcts]
    return ["#" + "".join([format(int(round(val*255)), "02x") for val in item]) for item in rgb_colors]

def round_up_to_odd(f):
    return int(np.ceil(f) // 2 * 2 + 1)

def right_intercept(interpollist):  #find the right crossing point of the x-axis and returns the index
    E_bulk_index = objekt.index(objekt.wavelength_energy(objekt.E_bulk),x_new)
    E_ml_index = objekt.index(objekt.wavelength_energy(objekt.E_ML),x_new)
    center = E_bulk_index + np.argmin(interpollist[E_bulk_index:E_ml_index])
    flag = False
    k = 1
    interpolright = np.argmax(interpollist[:center])
    while flag == False and k < center:
        if any(n > 0 for n in interpollist[:center]) and any(n < 0 for n in interpollist[:center]):     #check there are negative and values within the range of interest for there to be an intercept
            if interpollist[center - k] < 0:
                k += 1
            elif interpollist[center - k] > 0:
                interpolright = center - k
                flag = True
        else:
            flag = True
    return(interpolright)

def left_intercept(interpollist):
    E_bulk_index = objekt.index(objekt.wavelength_energy(objekt.E_bulk),x_new)
    E_ml_index = objekt.index(objekt.wavelength_energy(objekt.E_ML),x_new)
    center = E_bulk_index + np.argmin(interpollist[E_bulk_index:E_ml_index])
    flag = False
    k = 1
    interpolleft = np.argmax(interpollist[center:]) + center
    while flag == False and (center + k) < (int(len(interpollist))):
        if any(n > 0 for n in interpollist[center:]) and any(n < 0 for n in interpollist[center:]):  
            if interpollist[center + k] < 0:
                k += 1
            elif interpollist[center + k] > 0:
                interpolleft = center + k
                flag = True
        else:
            flag = True
    return(interpolleft)

def color_scheme(q, n):
    return([0.5*(1-q/n) , 0.5*(1-q/n) , 0.5+(q/(2*n)) , 0.7])

###

subset_size = round_up_to_odd(5/abs(data_space)**0.5)
print('Subset Size: ' + str(subset_size))

#fig,ax = plt.subplots(figsize=(6.2992,6.2992))

cm = 1/2.54

fig1 = plt.figure()
fig1.set_size_inches(16*cm,16*cm)


ax1 = plt.subplot2grid(shape=(15,2), loc=(0,0), rowspan=7)
ax2 = plt.subplot2grid(shape=(15,2), loc=(8,0), rowspan=7)
ax3 = plt.subplot2grid(shape=(15,2), loc=(0,1), rowspan=3)
ax4 = plt.subplot2grid(shape=(15,2), loc=(3,1), rowspan=3)
ax5 = plt.subplot2grid(shape=(15,2), loc=(6,1), rowspan=3)
ax6 = plt.subplot2grid(shape=(15,2), loc=(9,1), rowspan=3)
ax7 = plt.subplot2grid(shape=(15,2), loc=(12,1), rowspan=3)

fig1.subplots_adjust(hspace=0)
fig1.subplots_adjust(wspace=0.25)
fig1.subplots_adjust(left=0.125, right=0.975, top=0.975, bottom=0.1)

#Iterate the complete analysis process over each spectrum in the data file. Only the wavelength and smoothing window for each gamma value are saved globally. Every other data set is saved seperately for each spectra.
for i in range(0, len(spectra)):
#for i in range(0,1):
    lw_filtered=np.zeros([int(len(spectra[i])/8 + 1),len(spectra[0])])  #empty array to save the smoothed and differentiated spectra for error checking
    sec_dev=np.zeros([int(len(spectra[i])/8),len(spectra[0][lf_lim:rt_lim])])
    lw_filtered[0] = wavelength                                          #first entry to savelength the wavelength range
    sec_dev[0] = wavelength[lf_lim:rt_lim]

    rm_list = []
    r2_list = []
    windowlist= []
    fractionlist = []
    heading_string = str('Wavelength,')
    wavelength_a = []
    error_a = []
    diff_min = []
    y_diff_min = []
    rightzero = []
    leftzero = []
    exciton_energy = []
    for j in range(1, int(len(spectra[i])/8)):                                  #defines the lowess window range to be screened
        data_point_window = (j+0.5)*2
        windowlist.append(data_point_window)
        fractionlist.append(data_point_window/len(spectra[i]))
        heading_string = heading_string + 'wl:' + str('%.4f' %(data_point_window/len(spectra[i]))) + ','            #updates string to be used for headings in csv file
        smoothed = lowess(spectra[i],wavelength,frac=(data_point_window/len(spectra[i])),it=0,return_sorted=False)
        lw_filtered[j]= smoothed/(smoothed[objekt.index(objekt.normwavelength, wavelength)])
        sec_dev[j]=np.gradient(np.gradient(lw_filtered[j], wavelength),wavelength)[lf_lim:rt_lim]         #calculate second differential using x and y data points. Inclusion of x points crucial since spacing is NOT uniform.
        
        interpollist=np.zeros([len(spectra),1000])
        interpol= interpolate.interp1d(wavelength[lf_lim:rt_lim], sec_dev[j],kind="cubic")
        x_new=np.linspace(wavelength[lf_lim+1],wavelength[rt_lim-1],1000)
        interpollist=interpol(x_new)
        
        interpolright = right_intercept(interpollist)
        interpolleft = left_intercept(interpollist)
        rightzero.append(interpolright)
        leftzero.append(interpolleft)
        
        if interpolleft > interpolright:
            interpolminimum = (np.argmin(interpollist[interpolright:interpolleft])) + interpolright
        elif interpolright > interpolleft:
            interpolminimum = (np.argmin(interpollist[interpolleft:interpolright])) + interpolleft
        else:
            interpolminimum = 0
        diff_min.append(interpolminimum)
        y_diff_min.append(interpollist[interpolminimum])
            
        areas = []
        focused_wavelength = []
        
        for x in np.arange(start=interpolright,stop=interpolleft,dtype=int,step=1):
            areas.append(areamin(x)) 
        center = x_new[np.argmin(areas)+interpolright]
        five_percent = ((max(areas)-min(areas))*0.05)+min(areas)
        negative_error = x_new[interpolright + find_nearest(areas[np.argmin(areas):],five_percent) + np.argmin(areas)]
        positive_error = x_new[interpolright + find_nearest(areas[:np.argmin(areas)],five_percent)]
        wavelength_a.append(center)        #this is being calculated only for each smoothing window, more efficient.  
        error_a.append(max(positive_error-center,center-negative_error))
    
    right_smoothness = correlation_subset(windowlist, rightzero)
    left_smoothness = correlation_subset(windowlist, leftzero)
    center_smoothness = correlation_subset(windowlist, diff_min)
       
    
    solved_flag = False
    n = 0
    marker_size = 48
    
    while solved_flag != True and n < len(windowlist):
        if (x_new[rightzero[n]] - x_new[leftzero[n]]) > 15:
            if x_new[rightzero[n]] - x_new[diff_min[n]] > 5 and x_new[diff_min[n]] - x_new[leftzero[n]] > 5 :
                if right_smoothness[n] <= 0.02 and left_smoothness[n] <= 0.02 and center_smoothness[n] <= 0.02:
                    exciton_energy.append(wavelength_a[n])
                    solved_flag = True
                else: n += 1
            else:
                n += 1
        else:
            n += 1 
    
    def plot_data(ax, q, x_ticks=False, y_ticks=False):
        ax.plot(wavelength[lf_lim:rt_lim], sec_dev[q + 1], color=color_scheme(q, n), label=('Window: ' + str(windowlist[q])))
        ax.axhline(0, color='k', linewidth=1, linestyle=':')
        ax.scatter(x_new[diff_min[q]], y_diff_min[q], color='#FFA500', marker='o', edgecolors='black', s=marker_size)
        ax.scatter(x_new[rightzero[q]], 0, color='#C0504D', marker="^", edgecolors='black', s=marker_size)
        ax.scatter(x_new[leftzero[q]], 0, color='#96B954', marker="v", edgecolors='black', s=marker_size)
        ax.axvline(wavelength_a[q], color='#497AB9', linestyle='--', alpha=0.6)
        ax.fill_between(wavelength[find_nearest(wavelength, x_new[rightzero[q]]):find_nearest(wavelength, wavelength_a[q]) + 1],
                        y1=0, y2=sec_dev[q + 1][find_nearest(wavelength, x_new[rightzero[q]]) - lf_lim:find_nearest(wavelength, wavelength_a[q]) - lf_lim + 1], color=color_scheme(q, n), alpha=0.2)
    
        ax.set_yticks([]) if y_ticks else None
        ax.set_xticks([]) if x_ticks else None
    
        ax.set_ylabel("$\\frac{\partial^2 Ext.}{\partial \lambda^2}$", fontsize=12)
        ax.legend(frameon=False, fontsize=8)
    
    if i == 0:
        
        ax2.plot(windowlist, x_new[diff_min], color='#FFA500', linestyle=':', marker='.', label='Minimum', markersize=4)
        ax2.plot(windowlist, x_new[rightzero], color='#C0504D', linestyle=':', marker='^', label='Right Intercept', markersize=4)
        ax2.plot(windowlist, x_new[leftzero], color='#96B954', linestyle=':', marker='v', label='Left Intercept', markersize=4)
        ax2.plot(windowlist, wavelength_a, color=colour_list[i], linestyle=':', marker='x', label='Center', markersize=4)
        ax2.axhline(np.mean(exciton_energy), color='k', linestyle='--', alpha=0.8)
        ax2.axvline(windowlist[n], color=colour_list[i], linestyle='-.', alpha=0.8)
        ax2.set_xlabel('Smoothing Window', fontsize=10)
        ax2.set_ylabel('Wavelength / nm', fontsize=10)
        ax2.tick_params(axis='both', which='major', labelsize=8)
        ax2.set_xlim(0,50)
        ax2.set_ylim(640,690)
        ax2.legend(frameon=False, fontsize=8)
        ax2.text(-0.15, 0, 'b)', fontsize=12)
        
        ax7.plot(wavelength[lf_lim:rt_lim], sec_dev[n+1], color=color_scheme(n,n), label=('Window: ' + str(windowlist[n])))
        ax7.axhline(0, color='k', linewidth=1, linestyle=':')
        ax7.scatter(x_new[diff_min[n]], y_diff_min[n], color='#FFA500', marker='o', edgecolors='black', s=marker_size)
        ax7.scatter(x_new[rightzero[n]], 0, color='#C0504D', marker="^", edgecolors='black', s=marker_size)
        ax7.scatter(x_new[leftzero[n]], 0, color='#96B954', marker="v", edgecolor='black', s=marker_size)
        ax7.axvline(wavelength_a[n], color='#497AB9')
        ax7.fill_between(wavelength[find_nearest(wavelength,x_new[rightzero[n]]):find_nearest(wavelength,wavelength_a[n])+1],y1=0,y2=sec_dev[n+1][find_nearest(wavelength,x_new[rightzero[n]])-lf_lim:find_nearest(wavelength,wavelength_a[n])-lf_lim+1], color=color_scheme(n,n), alpha=0.2)
        ax7.set_yticks([])
        ax7.tick_params(axis='x', which='major', labelsize=8)
        ax7.set_xlabel('Wavelength / nm', fontsize=10)
        ax7.set_ylabel("$\\frac{\partial^2 Ext.}{\partial \lambda^2}$", fontsize=12)
        ax7.legend(frameon=False, fontsize=8)
        
        plot_data(ax3, 0, x_ticks=True, y_ticks=True)
        ax3.text(-0.15, 0, 'c)', fontsize=12)
        plot_data(ax4, n//5, x_ticks=True, y_ticks=True)
        plot_data(ax5, 2*n//5, x_ticks=True, y_ticks=True)
        plot_data(ax6, 3*n//5, x_ticks=True, y_ticks=True)
        
    ax1.errorbar(windowlist, wavelength_a, yerr=error_a, xerr=None, label=(labels[i]), color=colour_list[i])
    ax1.axvline(windowlist[n], color=colour_list[i], linestyle='-.', alpha=0.8)
    
ax1.axhline(np.mean(exciton_energy), color='k', linestyle='--', alpha=0.8)
ax1.set_ylabel('Exciton Position / nm', fontsize=10)
ax1.tick_params(axis='both', which='major', labelsize=8)
ax1.legend(frameon=False,fontsize=8)
ax1.set_xlim(0,200)
ax1.set_ylim(650,680)

plt.figtext(0.02, 0.96, 'a)', fontsize=14)
plt.figtext(0.02, 0.49, 'b)', fontsize=14)
plt.figtext(0.55, 0.96, 'c)', fontsize=14)

fig1.savefig('Combined Algorithm Illustration.svg')