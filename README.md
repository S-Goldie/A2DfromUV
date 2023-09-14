# A2DfromUV #
*Analysis of UV/VIS spectra from semiconducting 2D dispersions*

Python script for applying an automated analysis protocol to UV/VIS spectra of 2D nanoflake dispersions to extract thickness and length metrics. Automated smoothing and second derivative determination allows easy and reproducible analysis of typical spectra.

Included in the repository are the following:
* Instructions for use
* Source code in the `programm` folder
* Compiled program files suitable for download and use without any python installation (instructions below)
* Example data sets in `examples` folder

## Dependencies ##

The current version is verified to work on Python >=3.9.
Standard anaconda installation includes all required packages [click here for installation](https://docs.anaconda.com/free/anaconda/install/)

#### Required Packages ####

* numpy
* matplotlib
* tkinter
* scipy
*  statsmodels

## Installation ##

* Clone the repository `git clone https://github.com/SGoldie4/A2DfromUV/programm.git`
* Install all required packages (see section *Dependencies > Required Packages*, above).
* Launch the script from `main` ensuring the `GUI` and `properties` files are in the same folder

## Usage ##

Example data sets are available within the `example` folder

Ensure datasets are saved in `.csv` format, with **Wavelength** data in the leftmost column. Select the material and experiment parameters in the pop-up window, select the data file and the analysis will be completed automatically outputting the metrics and smoothed spectra in the same directory as the target file.

## References ##

If you use A2DfromUV in your work, please cite: *manuscript in preparation*
