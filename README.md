# A2DfromUV #
*Analysis of UV/VIS spectra from semiconducting 2D dispersions*

Python script for applying an automated analysis protocol to UV/VIS spectra of 2D nanoflake dispersions to extract thickness and length metrics. Automated smoothing and second derivative determination allows easy and reproducible analysis of typical spectra.

__Disclaimer:__ please note the metrics developed and output by this program are subject to limitations detailed in our published paper [doi coming soon]. Thickness and concentration values were fitted from dispersions in aqueous dispersion. We cannot guarantee their reliability in different solvents, especially those with different dielectric permittivity's. Length metrics were fitted for small nanosheets typically produced by LPE, but larger aspect ratio sheets may deviate from this trend.

Included in the repository are the following:
* Instructions for use
* Source code in the `program` folder
* Compiled program files suitable for download and use without any python installation (instructions below)
* Example data sets in `examples` folder

## Dependencies ##

The current version is verified to work on Python >=3.9.
A standard anaconda installation includes all required packages [click here for installation](https://docs.anaconda.com/free/anaconda/install/)

#### Required Packages ####

* numpy
* matplotlib
* tkinter
* scipy
* statsmodels

## Installation ##

#### Python Code ####

* Clone the repository `git clone https://github.com/SGoldie4/A2DfromUV/program.git`
* Alternatively, copy all files in the `Program` folder into the same folder on your machine and launch __main__ in python
* Install all required packages (see section *Dependencies > Required Packages*, above)
* Launch the script from `__main__` ensuring the `GUI` and `properties` files are in the same folder

#### Compiled .exe files ####
*Coming soon*
* Download the relevant program folder for your operating system
* Extract the files into a new folder, this will contain all the program files
* Launch the application from the .exe file within this folder (for convinience you may wish to make a Desktop shortcut to this .exe file)

Please be aware of the large file size of the program when downloaded and saved in this format; this is a limitation of compiling Python code into machine code. We have attempted to minimise this file size where possible, but be aware that the most up-to-date versions and best performance will be achieved by using the .py files within your own Python environment. When the .exe file is opened for the first time in a given session (i.e after restarting) it may take up to 30 seconds to load - within this time it may appear the program has frozen on the black command window.

## Usage ##

Example data sets are available within the `example` folder

Ensure datasets are saved with **Wavelength** data, **in nm**, in the leftmost column and extinction data in subsequent columns. For best results, use `.csv` files with sample labels in the top row above Wavelength columns. Best efforts have been made for flexibility in data format and `.txt` files with common delimiters are possible, but sample labels may become lost.

For analysis select the required material in the pop-up window, select the data file and the analysis will be completed automatically outputting the metrics and smoothed spectra in the same directory as the target file.
* Metrics: this file contains the calculated metrics as discussed in our paper. Concentration values are returned assuming a 1 cm pathlength, if a different path length is used correct according to Beer-Lambert.
* Spectra: this file contains the smoothed spectra and second derivative used for the analysis in `.csv` format.

If spectra are imported that do not include the required spectral range, this is will be indicated in the command window:
* If the required exciton transition is absent the program will terminate with an error warning.
* If the length or concentration range is absent this will be indicated in the `Metrics` file, but thickness determination will still be attempted.

## References ##

If you use A2DfromUV in your work, please cite: *manuscript in preparation*
