# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 11:58:33 2022

Properties: contains all the constants and material dependent properties to be imported by the "main" UV-vis analysis code
Constants use conventional units for discussing relative properties, e.g wavelengths are in nm, exiton energies in eV, 

@author: nicok
"""
from scipy.constants import h,c,e
import numpy as np


face="#738595"                                         #facecolor for figures
back="lightgrey"#"#363737"                             #backgroundcolor for figures
colors=['mediumorchid', 'blueviolet', 'navy', 'royalblue', 'darkslategrey',    #colorlist for plotting
        'limegreen', 'darkgreen', 'beige', 'khaki', 'orange', 'darkorange', "peru", "maroon"]



class Material_Methode:
    def __init__(self, material: str , methode: str):
        self._material = material
        self._methode = methode
    '''params for length metrics'''
        
    @property 
    def R1(self):
        if self._material == "WS2":
            return 235
        
        elif self._material == "MoS2":
            return 270
        
        elif self._material == "WSe2":
            return 235
        
        elif self._material == "MoSe2":
            return 280
        
        elif self._material == "NiPS3":
            return 370
        
        elif self._material == "RuCl3":
            return 370
        
        elif self._material == "PtSe2":
            return 800
        
    @property 
    def R2(self):
        if self._material == "WS2":
            return 295
        
        elif self._material == "MoS2":
            return 345
        
        elif self._material == "WSe2":
            return 347
        
        elif self._material == "MoSe2":
            return 390
        
        elif self._material == "NiPS3":
            return 460
        
        elif self._material == "RuCl3":
            return 555
        
        elif self._material == "PtSe2":
            return 255
        
    @property 
    def A1(self):
        if self._material == "WS2":
            if self._methode == "Ext":
                return 0.0159
            elif self._methode == "Abs":
                return 0.0252
        
        elif self._material == "MoS2":
            if self._methode == "Ext":
                return 0.0144
            elif self._methode == "Abs":
                return 0.0044
        
        elif self._material == "WSe2":
            if self._methode == "Ext":
                return 0.0143
            elif self._methode == "Abs":
                return 0.0219
        
        elif self._material == "MoSe2":
            if self._methode == "Ext":
                return 0.0117
            elif self._methode == "Abs":
                return 0.0064
        
        elif self._material == "NiPS3":
            if self._methode == "Ext":
                return -0.0185
            elif self._methode == "Abs":
                return -0.0009

        elif self._material == "RuCl3":
            if self._methode == "Ext":
                return -0.0019
            elif self._methode == "Abs":
                return -0.0018

        elif self._material == "PtSe2":
            if self._methode == "Ext":
                return 0.0134
            elif self._methode == "Abs":
                return 0.0141
            
    @property 
    def B1(self):
        if self._material == "WS2":
            if self._methode == "Ext":
                return 2.20
            elif self._methode == "Abs":
                return 2.66
        
        elif self._material == "MoS2":
            if self._methode == "Ext":
                return 1.97
            elif self._methode == "Abs":
                return 2.04
        
        elif self._material == "WSe2":
            if self._methode == "Ext":
                return 2.20
            elif self._methode == "Abs":
                return 2.64
        
        elif self._material == "MoSe2":
            if self._methode == "Ext":
                return 1.83
            elif self._methode == "Abs":
                return 1.79

        elif self._material == "NiPS3":
            if self._methode == "Ext":
                return 1.72
            elif self._methode == "Abs":
                return 1.69

        elif self._material == "RuCl3":
            if self._methode == "Ext":
                return 1.68
            elif self._methode == "Abs":
                return 1.75

        elif self._material == "PtSe2":
            if self._methode == "Ext":
                return -0.06
            elif self._methode == "Abs":
                return -0.11
            
    @property 
    def A2(self):
        if self._material == "WS2":
            if self._methode == "Ext":
                return 0.0166
            elif self._methode == "Abs":
                return 0.0270
        
        elif self._material == "MoS2":
            if self._methode == "Ext":
                return 0.0160
            elif self._methode == "Abs":
                return 0.0083
        
        elif self._material == "WSe2":
            if self._methode == "Ext":
                return 0.0200
            elif self._methode == "Abs":
                return 0.0308
        
        elif self._material == "MoSe2":
            if self._methode == "Ext":
                return 0.0154
            elif self._methode == "Abs":
                return 0.0103
        
        elif self._material == "NiPS3":
            if self._methode == "Ext":
                return 0.0003
            elif self._methode == "Abs":
                return 0.0011

        elif self._material == "RuCl3":
            if self._methode == "Ext":
                return -0.001
            elif self._methode == "Abs":
                return -0.0009

        elif self._material == "PtSe2":
            if self._methode == "Ext":
                return 0.0035
            elif self._methode == "Abs":
                return 0.0042
    
    @property 
    def B2(self):
        return 1
    
    """params for Thickness"""
      
    @property 
    def E_ML(self):
        if self._material == "WS2":
            if self._methode == "Ext":
                return 2.034
            elif self._methode == "Abs":
                return 2.033
        
        elif self._material == "MoS2":
            if self._methode == "Ext":
                return 1.896
            elif self._methode == "Abs":
                return 1.895
        
        elif self._material == "WSe2":
            if self._methode == "Ext":
                return 1.686
            elif self._methode == "Abs":
                return 1.692
        
        elif self._material == "MoSe2":
            if self._methode == "Ext":
                return 1.581
            elif self._methode == "Abs":
                return 1.599
        
        elif self._material == "NiPS3":
            if self._methode == "Ext":
                return 3.456
            elif self._methode == "Abs":
                return 3.414
        
        elif self._material == "RuCl3":
            if self._methode == "Ext":
                return 3.35
            elif self._methode == "Abs":
                return 3.305
        
        elif self._material == "PtSe2":
            if self._methode == "Ext":
                return 1.97
            elif self._methode == "Abs":
                return 3.508
            
    @property 
    def E_bulk(self):
        if self._material == "WS2":
            if self._methode == "Ext":
                return 1.957
            elif self._methode == "Abs":
                return 1.966
        
        elif self._material == "MoS2":
            if self._methode == "Ext":
                return 1.828
            elif self._methode == "Abs":
                return 1.846
        
        elif self._material == "WSe2":
            if self._methode == "Ext":
                return 1.610
            elif self._methode == "Abs":
                return 1.626
        
        elif self._material == "MoSe2":
            if self._methode == "Ext":
                return 1.516
            elif self._methode == "Abs":
                return 1.544
        
        elif self._material == "NiPS3":
            if self._methode == "Ext":
                return 3.139
            elif self._methode == "Abs":
                return 3.243
        
        elif self._material == "RuCl3":
            if self._methode == "Ext":
                return 3.2
            elif self._methode == "Abs":
                return 3.262
        
        elif self._material == "PtSe2":
            if self._methode == "Ext":
                return 1.81
            elif self._methode == "Abs":
                return 3.206
            
    @property 
    def N0(self):
        if self._material == "WS2":
            if self._methode == "Ext":
                return 3.69
            elif self._methode == "Abs":
                return 3.72
        
        elif self._material == "MoS2":
            if self._methode == "Ext":
                return 5.72
            elif self._methode == "Abs":
                return 3.37
        
        elif self._material == "WSe2":
            if self._methode == "Ext":
                return 5.15
            elif self._methode == "Abs":
                return 3.29
        
        elif self._material == "MoSe2":
            if self._methode == "Ext":
                return 8.59
            elif self._methode == "Abs":
                return 2.80
            
        elif self._material == "NiPS3":
            if self._methode == "Ext":
                return 5.8
            elif self._methode == "Abs":
                return 2.92
            
        elif self._material == "RuCl3":
            if self._methode == "Ext":
                return 5.30
            elif self._methode == "Abs":
                return 6.47

        elif self._material == "PtSe2":
            if self._methode == "Ext":
                return 3.42
            elif self._methode == "Abs":
                return 5.08
       
    @property 
    def normwavelength(self):
        if self._material == "WS2":
            return 295
            
        elif self._material == "MoS2":
            return 345
        
        elif self._material == "WSe2":
            return 347
        
        elif self._material == "MoSe2":
            return 390
        
        elif self._material == "NiPS3":
            return 460
        
        elif self._material == "RuCl3":
            return 555
        
        elif self._material == "PtSe2":
            return 1149
            
    
    def index(self, lamb:float, liste: list)->int:                       #function to get index from given value out of list
        self._lamb = lamb
        self._liste = liste
        return int(np.argmin(abs(self._liste - self._lamb)))


    def energy_wavelength(self,lamb : float)->float:                      #calculate energy in eV to given wavelength in nm           
        return (h*c)/(lamb*10**-9*e)

    
    def wavelength_energy(self,lamb : float)->float:                    #calculate wavelength in nm given energy in eV
        return (h*c*10**9)/(lamb*e)


    def length (self, spectrum_y: list, wavelength_x : list)->float: 
        self._spectrum = spectrum_y
        self._wavelength = wavelength_x
        
        Ex1 = self._spectrum[self.index(self.R1,self._wavelength)]
        Ex2 = self._spectrum[self.index(self.R2,self._wavelength)]
        length = (self.B1*Ex2-self.B2*Ex1)/(self.A2*Ex1-self.A1*Ex2)
        return length
    

    def thickness(self, exitonwavelength : float)->float: 
        self._exitonwavelength = exitonwavelength
        exitonenergy=self.energy_wavelength(self._exitonwavelength)
        
        dif=exitonenergy-self.E_bulk
        if self._material == "WS2":
            return 0.88 + 2.311998E31*np.exp(-35.62*exitonenergy)
                    
        else:
            if dif >0:
                return abs(1+self.N0*np.log(((self.E_ML-self.E_bulk)/(dif))))     
            else: 
                return 50
    


    