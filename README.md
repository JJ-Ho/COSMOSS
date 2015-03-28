# COSMOSS
==========================================================================
COSMOSS = Coupled OScillator MOdel Spectrum Simulator  

Given molecule structure, this code can generates different kinds of vibrational spectrum 
base on coupled oscillator model. For now COSMOSS can be used to simulate:
1. Fourier Transform Infared spectrum (FTIR)
2. Sum-Frequency Generation spectrum (SFG)
3. Two dimentional Infared spectrum (2DIR)
4. Two dimentional Sum-Frequency Generation spectrum (2DSFG)

The supported structure generation functions include:
1. Two coupled oscillators
2. PDB_AmideI
This code is design to be used for other molecules as well. You can generate your interest molecule
and their spectra. Please check sub-function in "MoleculeConstruction" folder for more infomation. 

To use this code, please download it into your Matlab path and excute COSMOSS.m.
==========================================================================
Note: 
The GUI portion is created base on the "GUI Layout Toolbox" wrote by David Sampson. 
For this version (v0.6) you have to install this tool box to use COSMOSS. In near future, 
the Matlab default GUI layout will be supported, too.

Further information about GUI Layout Tool box can be found in: 
  http://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox
==========================================================================
jjho 15/03/28
