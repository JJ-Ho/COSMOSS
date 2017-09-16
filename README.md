# COSMOSS 
COSMOSS = *Coupled OScillator MOdel Spectrum Simulator*

### DESCRIPTION
Given a molecule structure, this Matlab code can simulate vibrational spectrum base on coupled oscillator model (Frankle Exiton Model).  

COSMOSS can be used to simulate:
* Fourier Transform Infrared spectrum (FTIR)
* Sum-Frequency Generation spectrum (SFG)
* Two-dimensional Infrared spectrum (2DIR)
* Two-dimensional Sum-Frequency Generation spectrum (2DSFG)

For more information check [Wiki page](../../wiki)

Currently, COSMOSS support the following structure construction models:
* Two coupled oscillators
* Extract Amide-I mode from PDB file
* Build 2D-Grid from formatted quantum simulation outputs, such as Gaussian 09. 
* Amide-I mode of an ideal beta-sheet
* Combination of any two structures above

This code is designed to be used for other molecules as well. You can generate your interested molecule and their spectra. Check [Wiki page](../../wiki) for more information.

### Installation
To use COSMOSS, please clone or download it from GitHub and place the whole folder in your Matlab path.  Executing the "COSMOSS.m" file will bring up the main GUI.

You can also download the application package(maybe outdated) from Matlab file exchange. 

https://www.mathworks.com/matlabcentral/fileexchange/64433-cosmoss

The package is integrated with the Matlab App Installation process. Once installed, you can directly access COSMOSS from your App banner of Matlab.



Note:
The GUI is built upon the "GUI Layout Toolbox",  written by David Sampson. Further information about GUI Layout Toolbox can be found in: [GUI Layout Tool Box][GUILayoutToolbox]

### License

GPL2



[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)
   [wiki-Home]: https://gitlab.com/jjho/COSMOSS/wikis/home
   [GUILayoutToolbox]: http://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox
   
