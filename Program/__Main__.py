# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 18:27:19 2023

Main: imports material properties and data from csv file, automatically calculates different smoothing parameters and associated loss function. Systematic comparison of different loss functions on smoothing.
Production of plots contrasting different smoothing paramters and loss functions.

Version: 1.2.2 - algorithm to identify the minimum smoothing required using the variability of the x-axis intercepts and peak minimum in second derivative space, as well as minimum peak width and distance between features.
Once minimum smoothing identified, center of integral area used to determine the peak position in nanometers which can be converted into energy and used for metric calculation.
Current metrics: old ln relation from Kevin's thesis, updated fitting required
Estimate of error in wavelength updated, new experimental form to use the gradient of the data outside the intercepts to estimate noise and the peak intensity to estimate signal.

@author: Stuart Goldie, Nico Kubetschek
"""

"""------Import of the used modules--------------------------------------------"""
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
from tkinter.filedialog import askopenfilename
from scipy import interpolate
from properties import *
from tkinter import Tk
import GUI as gui
"""----------------------------------------------------------------------------"""

import_flag = False

"""------Import and sort data with format and range checks for error handling-------"""
while import_flag != True:
    Tk().withdraw()
    Filename = askopenfilename()    
    
    if '.csv' in Filename:
        short_name = Filename[Filename.rfind('/')+1:Filename.find('.csv')]
    elif '.txt' in Filename:
        short_name = Filename[Filename.rfind('/')+1:Filename.find('.txt')]
    
    print("File loaded:", short_name)
    #%%
    
    objekt=Material_Methode(gui.material, gui.methode)
    format_error = False
    
    try:
        file_found = False
        rows_skipped = 0

        while file_found != True:
            if '.csv' in Filename:
                trial_in = np.genfromtxt(Filename, skip_header=rows_skipped, max_rows=1, delimiter=",")
                if np.isnan(trial_in[0]):
                    rows_skipped += 1
                else:    
                    file_found = True
                    data = np.genfromtxt(Filename, skip_header=rows_skipped, delimiter=",")
            
            elif '.txt' in Filename:
                try:
                    trial_in = np.genfromtxt(Filename, skip_header=rows_skipped, max_rows=1, delimiter=None)
                    if np.isnan(trial_in[0]):
                        rows_skipped += 1
                    else:    
                        file_found = True
                        data = np.genfromtxt(Filename, skip_header=rows_skipped, delimiter=None)
                except:
                    print('File format could not be recognised')
                    break
        
        #Preparation and sorting of the data
        wavelength=data[:,0]       
        spectra=np.zeros([data.shape[1],data.shape[0]])
                
        count0=0
        data_space = abs(np.mean(wavelength[1:]-wavelength[:-1]))  
        
        for i in range(data.shape[1]):
            if sum(data[:,i])/len(data[:,i]) <= 200:
                for j in range (data.shape[0]):
                    spectra[i][j] = data[j,i]
            if sum(spectra[i])==0.0:         
                count0+=1
                   
        try:
            for i in range(len(spectra)-count0+1):
                if sum(spectra[i])== 0.0:
                    spectra = np.delete(spectra,(i), axis = 0) 
        except IndexError:
            for i in range(len(spectra)-count0): 
                if sum(spectra[i])== 0.0:
                    spectra = np.delete(spectra,(i), axis = 0)
        
        
        #Confirm spectra range is valid
        if max(wavelength) > (objekt.wavelength_energy(objekt.E_bulk) + 45) and min(wavelength) < (objekt.wavelength_energy(objekt.E_ML) - 45):
            neededminwavelength=min([objekt.R1,objekt.R2,objekt.normwavelength,objekt.coeff[0]])
            minwavelength = min(wavelength)
            if minwavelength > neededminwavelength:
                print("Caution! Spectrum energy range does not fit to required energy for concentration or flake length determination")
                gui.open_second_window()
                if gui.overide_flag == True:
                    length_flag =False
                    import_flag = True
                elif gui.overide_flag == False:
                    import_flag = False
            else:
                length_flag = True
                import_flag = True        
        
        else:
            print("Caution! Spectrum energy range does not include resonance energy")
          
    except UnicodeDecodeError:
        print("Non-data table file type selected")
    
    except IndexError:
        print("Data format not recognised. Please see GitHub page for information on correct file format")
    
    except:
        print("The file could not be read. \n Ensure .csv format and correct layout")
    

"""----------------------------------------------------------------------------"""

"""-----create color gradient in dependence to subset size and load labels-----"""
colors=[]
for i in range(1,spectra.shape[0]+1):
    colors.append([1-i/(spectra.shape[0]),0.1,i/(spectra.shape[0]),0.7])


try:
    if '.csv' in Filename:
        labels = np.loadtxt(Filename, delimiter=",", dtype=str, max_rows=1)
    elif '.txt' in Filename:
        labels = np.loadtxt(Filename, delimiter=None, dtype=str, max_rows=1)

    for i in range(len(labels)-count0):
        if labels[i]=="":
            labels = np.delete(labels,(i), axis = 0) 
except:
    labels = range(1, len(spectra))
"""----------------------------------------------------------------------------"""


"""
Definition of left and right limit of analysis (for conventional measurement 
from high to low wavelength, left will be higher wavelengths). 
Note: smoothing and differentiation is applied to the entire spectrum to minimise boundary errors.
Use of energy (function in properties) to correct to wavelength scale.

The borders confine a regime in the wavelengthspace with a size of 90 nm around 
the wavelength which corresponds to the estimated bulk A exciton transition energy.
"""
    
lf_lim = ((objekt.index(objekt.wavelength_energy(objekt.E_bulk),wavelength)) - int(45/abs(data_space)))      #lower index,  which means higher wavelength
rt_lim = ((objekt.index(objekt.wavelength_energy(objekt.E_ML),wavelength)) + int(45/abs(data_space)))        #higher index, which means lower wavelength


"""------------------------Functions-------------------------------------------"""
def correlation_subset(j, y):
    """
    Calculates a subset size in which the variation after smoothing is calculated
    """
    subset_y = y[j:]
    if np.std(subset_y) != 0:
        return(float((np.std(subset_y)/np.mean(subset_y))))
    else:
        return(float(0.0))

def areamin (mini, interpollist, interpolright, interpolleft):
    """
    Minimisation of the integral area either side of the central value. By minimising the difference, the center of mass of the peak can be calculated
    """
    arealeft=np.trapz(interpollist[interpolright:int(mini)],x_new[interpolright:int(mini)])
    arearight=np.trapz(interpollist[int(mini):interpolleft],x_new[int(mini):interpolleft])
    difference=abs(arealeft-arearight)
    return difference

def find_nearest(array, value):
    """This function returns the index of an value from an array"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def round_up_to_odd(f):
    """
    The function  takes a floating-point number f as input and returns the nearest odd integer that is greater than or equal to f.
    """
    return int(np.ceil(f) // 2 * 2 + 1)

def right_intercept(interpollist):
    """ This function finds the right crossing point of the x-axis and returns the index"""
    E_bulk_index = objekt.index(objekt.wavelength_energy(objekt.E_bulk),x_new)
    E_ml_index = objekt.index(objekt.wavelength_energy(objekt.E_ML),x_new)
    center = E_bulk_index + np.argmin(interpollist[E_bulk_index:E_ml_index])
    flag = False
    m = 1
    interpolright = np.argmax(interpollist[:center])
    while flag == False and m < center:
        if any(n > 0 for n in interpollist[:center]) and any(n < 0 for n in interpollist[:center]):     #check there are negative and values within the range of interest for there to be an intercept
            if interpollist[center - m] < 0:
                m += 1
            elif interpollist[center - m] > 0:
                interpolright = center - m
                flag = True
        else:
            flag = True
    return(interpolright)

def left_intercept(interpollist):
    """ This function finds the left crossing point of the x-axis and returns the index"""
    E_bulk_index = objekt.index(objekt.wavelength_energy(objekt.E_bulk),x_new)
    E_ml_index = objekt.index(objekt.wavelength_energy(objekt.E_ML),x_new)
    center = E_bulk_index + np.argmin(interpollist[E_bulk_index:E_ml_index])
    flag = False
    m = 1
    interpolleft = np.argmax(interpollist[center:]) + center
    while flag == False and (center + m) < (int(len(interpollist))):
        if any(n > 0 for n in interpollist[center:]) and any(n < 0 for n in interpollist[center:]):  
            if interpollist[center + m] < 0:
                m += 1
            elif interpollist[center + m] > 0:
                interpolleft = center + m
                flag = True
        else:
            flag = True
    return(interpolleft)

"""----------------------------------------------------------------------------"""

subset_size = (round_up_to_odd(5/abs(data_space)**0.5))                        #automatically calculate the subset size based on the data spacing. This was empirically observed by studying the behaviour of the function for many data sets.
print('Subset Size: ' + str(subset_size))

"""-----Creation of helplists to contain output values from each spectrum------"""
final_wavelength_a   =   []
final_error_a        =   []
final_window_list    =   []
final_fraction_list  =   []
final_flake_length   =   []
"""----------------------------------------------------------------------------"""


"""-------------Helparrays for datastorage-------------------------------------"""
lw_filtered = np.zeros([spectra.shape[0],spectra.shape[1],spectra.shape[1]])   #used to save final processed spectra for illustrative plotting
sec_diff = np.zeros([spectra.shape[0],spectra.shape[1],spectra.shape[1]])      #array of differentiated data, [i,j,x] where i is the spectrum index, j is the smoothness index, x is the wavelength index
interpol_array = np.zeros([spectra.shape[0],spectra.shape[1],1000])            #array for interpolated data storage, [i,j,x] where i is the spectrum index, j is the smoothness index, x is the interpolated wavelength index
x_new=np.linspace(wavelength[lf_lim+1],wavelength[rt_lim-1],1000)

processed_spectra = np.zeros([int(1+2*spectra.shape[0]), int(len(wavelength))])#empty array to save the exciton wavelength and optimal window for each gamma value
processed_spectra[0] = wavelength
"""----------------------------------------------------------------------------"""

for i in range(spectra.shape[0]):                                              #Iterate over each spectra in turn for analysis
    print("Processing spectrum:", i+1)
    diff_min = []
    rightzero = []
    leftzero = []
    local_window_list = []
    local_fraction_list = []
    
    for j in range(0,subset_size):                                             #Iterate over the sub-set to produce enough data points for varation analysis
        data_point_window = (j+1.5)*2                                          #j is an index, here we increase smoothing window length by 2, starting at 3 points when j=0
        local_window_list.append(data_point_window)
        local_fraction_list.append(data_point_window/len(spectra[i]))
        smoothed = lowess(spectra[i], wavelength, frac=(data_point_window/len(spectra[i])),it=0,return_sorted=False)
        if length_flag == True:
            lw_filtered[i,j]= smoothed/(smoothed[objekt.index(objekt.normwavelength, wavelength)])
        elif length_flag == False:
            lw_filtered[i,j]=smoothed/max(smoothed)
        sec_diff[i,j] = np.gradient(np.gradient(lw_filtered[i,j], wavelength),wavelength)
        interpol= interpolate.interp1d(wavelength[lf_lim:rt_lim], sec_diff[i,j][lf_lim:rt_lim],kind="cubic")
        interpol_array[i,j] = interpol(x_new)
        
        rightzero.append(right_intercept(interpol_array[i,j]))  #use functions to find the left and right intercept - currently causes error for some samples wth a wide possible peak range because the initial centre can be positive to start.
        leftzero.append(left_intercept(interpol_array[i,j]))
        diff_min.append(np.argmin(interpol_array[i,j,rightzero[j]:leftzero[j]]) + rightzero[j])
        
    #First smoothing points calculated. At least the initial sub-set must be calculated for the analysis
    #at this point j is equal to subset size. Analysis to start with the first point, and extra smoothing and differntial spectra are only calculated as needed
    convergence_flag = False
    k = 0
    plottable_right_smoothness = []
    plottable_left_smoothness = []
    plottable_center_smoothness = []
    
    try:
        while convergence_flag != True:
            if k == 999:
                print("Caution: " + labels[i] + " could not be analysed.")
                final_wavelength_a.append("--")
                final_error_a.append("--")
                final_window_list.append(local_window_list[k-1])
                final_fraction_list.append(local_fraction_list[k-1])
                if length_flag == True:
                    final_flake_length.append("--")
                elif length_flag == False:
                    final_flake_length.append('Insufficient Data')
                break
            right_smoothness = correlation_subset(k, rightzero)
            plottable_right_smoothness.append(right_smoothness)
            left_smoothness = correlation_subset(k, leftzero)
            plottable_left_smoothness.append(left_smoothness)
            center_smoothness = correlation_subset(k, diff_min)
            plottable_center_smoothness.append(center_smoothness)
            if (x_new[rightzero[k]] - x_new[leftzero[k]]) > 15:     #a peak narrower than 15nm is assumed to be noise crossing the x-axis
                if x_new[rightzero[k]] - x_new[diff_min[k]] > 5 and x_new[diff_min[k]] - x_new[leftzero[k]] > 5 :      #the intercepts must be either side and ruther than 5nm from the minimal point. (Using the minimum point like this can cause errors with spectra that have discontinuities at localised points due to source or detector changes)
                    if right_smoothness <= 0.025 and left_smoothness <= 0.025 and center_smoothness <= 0.025:         #all metrics must show less than 2.5% variation within a sub-set size defined from the data interval.
                        areas = []
                        focussed_wavelength = []
                        for x in np.arange(start=rightzero[k],stop=leftzero[k],dtype=int,step=1):
                            areas.append(areamin(x, interpol_array[i,k], rightzero[k], leftzero[k]))
                        center = x_new[np.argmin(areas)+rightzero[k]]
                        
                        noise_list = np.gradient(sec_diff[i,k][objekt.index(x_new[leftzero[k]],wavelength):rt_lim])      #calculate the gradient as an approximation of the noise in the peak region, to the left of the x-intercept of the peak
                        np.append(noise_list, np.gradient(sec_diff[i,k][lf_lim:objekt.index(x_new[rightzero[k]],wavelength)]))
                        noise = np.mean(np.absolute(noise_list))
                        signal = np.min(sec_diff[i,k][objekt.index(x_new[rightzero[k]],wavelength):objekt.index(x_new[leftzero[k]],wavelength)])
                        sn = np.abs(signal / noise)
                        #print('Signal to Noise: {:.2f}'.format(sn))
                        five_percent = ((max(areas)-min(areas))*0.05)+min(areas)
                        negative_error = x_new[rightzero[k] + find_nearest(areas[np.argmin(areas):],five_percent) + np.argmin(areas)]
                        positive_error = x_new[rightzero[k] + find_nearest(areas[:np.argmin(areas)],five_percent)]
                        final_wavelength_a.append(center)        #this is being calculated only for each smoothing window, more efficient.  
                        final_error_a.append(max(positive_error-center,center-negative_error) * (1+(1/sn)))
                        #print('5% Peak Minimum Error: {:.2f}'.format(max(positive_error-center,center-negative_error)))
                        #print('Combined Error: {:.2f}'.format(final_error_a[i]))
                        
                        final_window_list.append(local_window_list[k])
                        final_fraction_list.append(local_fraction_list[k])
                        if length_flag == True:
                            final_flake_length.append(objekt.length(lw_filtered[i,k,:], wavelength))
                        elif length_flag == False:
                            final_flake_length.append('Insufficient Data')
                        processed_spectra[1+i] = lw_filtered[i,k,:]
                        processed_spectra[1+len(spectra)+i] = sec_diff[i,k,:]
                        convergence_flag = True
                    else:
                        l = k + subset_size
                        data_point_window = (l+1.5)*2
                        local_window_list.append(data_point_window)
                        local_fraction_list.append(data_point_window/len(spectra[i]))
                        smoothed = lowess(spectra[i], wavelength, frac=(data_point_window/len(spectra[i])),it=0,return_sorted=False)
                        lw_filtered[i,l]= smoothed/(smoothed[objekt.index(objekt.normwavelength, wavelength)])
                        sec_diff[i,l] = np.gradient(np.gradient(lw_filtered[i,l], wavelength),wavelength)
                        interpol= interpolate.interp1d(wavelength[lf_lim:rt_lim], sec_diff[i,l][lf_lim:rt_lim],kind="cubic")
                        interpol_array[i,l] = interpol(x_new)
                        
                        rightzero.append(right_intercept(interpol_array[i,l]))
                        leftzero.append(left_intercept(interpol_array[i,l]))
                        diff_min.append(np.argmin(interpol_array[i,l,rightzero[l]:leftzero[l]]) + rightzero[l])
                        k += 1
                        
                else:
                    l = k + subset_size
                    data_point_window = (l+1.5)*2
                    local_window_list.append(data_point_window)
                    local_fraction_list.append(data_point_window/len(spectra[i]))
                    smoothed = lowess(spectra[i], wavelength, frac=(data_point_window/len(spectra[i])),it=0,return_sorted=False)
                    if length_flag == True:
                        lw_filtered[i,l]= smoothed/(smoothed[objekt.index(objekt.normwavelength, wavelength)])
                    elif length_flag == False:
                        lw_filtered[i,l]= smoothed/max(smoothed)
                    sec_diff[i,l] = np.gradient(np.gradient(lw_filtered[i,l], wavelength),wavelength)
                    interpol= interpolate.interp1d(wavelength[lf_lim:rt_lim], sec_diff[i,l][lf_lim:rt_lim],kind="cubic")
                    interpol_array[i,l] = interpol(x_new)
                    
                    rightzero.append(right_intercept(interpol_array[i,l]))
                    leftzero.append(left_intercept(interpol_array[i,l]))
                    diff_min.append(np.argmin(interpol_array[i,l,rightzero[l]:leftzero[l]]) + rightzero[l])
                    k += 1
            else:
                l = k + subset_size
                data_point_window = (l+1.5)*2
                local_window_list.append(data_point_window)
                local_fraction_list.append(data_point_window/len(spectra[i]))
                smoothed = lowess(spectra[i], wavelength, frac=(data_point_window/len(spectra[i])),it=0,return_sorted=False)
                if length_flag == True:
                    lw_filtered[i,l]= smoothed/(smoothed[objekt.index(objekt.normwavelength, wavelength)])
                elif length_flag == False:
                    lw_filtered[i,l]= smoothed/max(smoothed)
                sec_diff[i,l] = np.gradient(np.gradient(lw_filtered[i,l], wavelength),wavelength)
                interpol= interpolate.interp1d(wavelength[lf_lim:rt_lim], sec_diff[i,l][lf_lim:rt_lim],kind="cubic")
                interpol_array[i,l] = interpol(x_new)
                
                rightzero.append(right_intercept(interpol_array[i,l]))
                leftzero.append(left_intercept(interpol_array[i,l]))
                diff_min.append(np.argmin(interpol_array[i,l,rightzero[l]:leftzero[l]]) + rightzero[l])
                k += 1
    
    except ValueError:
        print("Caution: " + labels[i] + " could not be analysed.")
        final_wavelength_a.append("--") 
        final_error_a.append("--")
        final_window_list.append(local_window_list[k])
        final_fraction_list.append(local_fraction_list[k])
        if length_flag == True:
            final_flake_length.append("--")
        elif length_flag == False:
            final_flake_length.append('Insufficient Data')


#%% At this point, smooth parameters have been identified and recorded in small lists denoted final.
#Intermediate analysis of smoothing and differentation are available in 3D arrays [i,j,k] where i is the spectra number, j the smoothing index and k the wavelength index.
#Plot the final smoothed spectra and second derivative with the A-exciton position highlighted.
   
plt.figure("Spectral_Analyis Smoothness",figsize=(10,5))
   
for i in range(len(spectra)):
    plt.subplot(1,2,1)
    plt.title('Spectrum', fontsize = 14)
    
    plt.plot(wavelength, lw_filtered[i,int((final_window_list[i]/2)-1.5),:], color=colors[i], label=labels[i])
    if final_wavelength_a[i] != '--':
        plt.axvline(final_wavelength_a[i], color=colors[i])
    plt.xlabel("$\lambda / nm$",fontsize=14)
    plt.ylabel("$Ext$",fontsize=14)
    plt.xlim(wavelength[-1],wavelength[0])   
    plt.legend()
    
    plt.subplot(1,2,2)
    plt.title('Differentiated',fontsize=14)
    
    if final_wavelength_a[i] != 0:
        plt.axvline(final_wavelength_a[i],color=colors[i])
        plt.plot(x_new,interpol_array[i,int((final_window_list[i]/2)-1.5),:],color=colors[i], linestyle='-', label=("$\lambda_A$: {:.2f} nm".format(final_wavelength_a[i])),linewidth=2,alpha=0.4)
    
    plt.xlim(wavelength[rt_lim+5],wavelength[lf_lim-5])
    plt.xlabel("$\lambda / nm$",fontsize=14)
    plt.ylabel("$\\frac{\partial^2 Ext.}{\partial \lambda^2}$",fontsize=14)
    plt.axhline(0,color="k")
    plt.legend()


names = []
for i in range(len(spectra)):
    names.append(labels[i])
    
    
final_energy_a = []
final_energy_a_error = []
final_thickness = []
final_concentration = []
final_Nv_error = []
     
for m in np.arange(len(final_wavelength_a)):
    if final_wavelength_a[m] == "--":
        final_energy_a.append("--")
        final_energy_a_error.append("--")
        final_thickness.append("--")
        final_Nv_error.append("--")     
    else:
        final_energy_a.append(objekt.energy_wavelength(final_wavelength_a[m]))
        final_energy_a_error.append(final_energy_a[m] * final_error_a[m] / final_wavelength_a[m])
        final_thickness.append(objekt.thickness(final_wavelength_a[m]))
        final_Nv_error.append(objekt.thickness_error(final_energy_a[m],final_energy_a_error[m]))
    if objekt.coeff == 0:
        final_concentration.append("Not Available")
    elif length_flag == False:
        final_concentration.append('Insufficient Data')
    elif length_flag == True:
        final_concentration.append((spectra[m][objekt.index(objekt.coeff[0],wavelength)])/objekt.coeff[1])
output_headings = str('Sample,Exciton Wavelength / nm,Error / nm,Exciton Energy / eV,Error / eV,<N>vf,<N>vf error,Flake Length / nm,Smoothing Window,Smoothing Fraction,Concentration g/(L*cm)')

try:
    np.savetxt(Filename.removesuffix('.csv') + '_Metrics.csv', np.column_stack((names, final_wavelength_a, final_error_a, final_energy_a, final_energy_a_error, final_thickness,final_Nv_error, final_flake_length, final_window_list, final_fraction_list,final_concentration)), delimiter=',', header=output_headings, fmt='%s')
except AttributeError:
    np.savetxt(Filename[:-4] + '_Metrics.csv', np.column_stack((names, final_wavelength_a, final_error_a, final_energy_a, final_energy_a_error, final_thickness, final_Nv_error, final_flake_length, final_window_list, final_fraction_list,final_concentration)), delimiter=',', header=output_headings, fmt='%s')
processed_headings = "Wavelength / nm,"
for i in range(len(spectra)):
    processed_headings = processed_headings + str(labels[i]) + '-Smoothed,'
for i in range(len(spectra)):
    processed_headings = processed_headings + str(labels[i]) + '-Differentiated,'
try:
    processed_headings.removesuffix(',')        #cleans format by removing final comma at end of string
except:
    print('Upgrade to Python3.9 for string cleaning. \n The output file may contain additional empty columns')
    
processed_output = [list(n) for n in zip(*processed_spectra)]
try:
    np.savetxt(Filename.removesuffix('.csv') + '_spectra.csv', (processed_output), delimiter=',', header=processed_headings)
except AttributeError:
    np.savetxt(Filename[:-4] + '_spectra.csv', (processed_output), delimiter=',', header=processed_headings)

print('Data saved')
plt.show()
input('Press Enter or close window to exit program')
