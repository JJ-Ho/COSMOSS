# COSMOSS 
COSMOSS = *Coupled OScillator MOdel Spectrum Simulator*

### DESCRIPTION
Given a molecule structure, this Matlab code can simulate vibrational spectrum base on coupled oscillator model (Frankle Exiton Model).  
For now COSMOSS can be used to simulate:
* Fourier Transform Infared spectrum (FTIR)
* Sum-Frequency Generation spectrum (SFG)
* Two dimentional Infared spectrum (2DIR)
* Two dimentional Sum-Frequency Generation spectrum (2DSFG)

For more information check [Wiki page](wiki-Home)

Currently COSMOSS support the following structure construction models:
* Two coupled oscillators
* Extract Amide-I mode from PDB file
* Build 2D-Grid from formatted quantum simulation ouput
* Amide-I mode of an ideal betasheet
* Combination of any two structures above

This code is designed to be used for other molecules as well. You can generate your interested molecule and their spectra. Check [structure mode wiki page](wiki-Structure) for more infomation.

### Installation
To use this code, please download COSMOSS into your Matlab path and excute COSMOSS.m.

Note:
The GUI portion is created base on the "GUI Layout Toolbox" wrote by David Sampson. Further information about GUI Layout Tool box can be found in: [GUI Layout Tool Box](GUILayoutToolbox)

### License

GPL2



[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)
   [wiki-Home]: <https://gitlab.com/jjho/COSMOSS/wikis/home>
   [wiki-Structure]: <https://gitlab.com/jjho/COSMOSS/wikis/Structure-Model>
   [GUILayoutToolbox]: <http://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox>
   
