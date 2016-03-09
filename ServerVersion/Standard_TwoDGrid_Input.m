function Output = Standard_TwoDGrid_Input

% For Monomer -------------------------------------------------------------
Output.Ang_Phi     =  0; %NUMBER, Angle Degrees%
Output.Ang_Psi     =  0; %NUMBER, Angle Degrees%
Output.Ang_Theta   =  0; %NUMBER, Angle Degrees%
Output.Delta_Phi   =  0; %NUMBER, Angle Degrees%
Output.Delta_Psi   =  0; %NUMBER, Angle Degrees%
Output.Delta_Theta =  0; %NUMBER, Angle Degrees%
% -------------------------------------------------------------------------

% For 2D Grid -------------------------------------------------------------
Output.Vec_1 = [7,0,0]; %ARRAY, Angstron%
Output.Vec_2 = [0,4,0]; %ARRAY, Angstron%
Output.N_1   = 3; %NUMBER%
Output.N_2   = 3; %NUMBER%
% -------------------------------------------------------------------------

% For Local Mode frequencies ----------------------------------------------
Output.NLFreq  = 1720; %NUMBER, cm^-1%
Output.LFreq   = 1750; %NUMBER, cm^-1%
Output.Anharm  =   20; %NUMBER, cm^-1%
Output.L_Index =   []; %ARRAY%
% -------------------------------------------------------------------------