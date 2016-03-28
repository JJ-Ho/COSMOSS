function Output = Standard_Main_Input
% Frequency Range ---------------------------------------------------------
Output.F_Min     = 1550; %NUMBER, cm^-1%
Output.F_Max     = 1800; %NUMBER, cm^-1%
Output.FreqRange = Output.F_Min:Output.F_Max;
% ------------------------------------------------------------------------

% Coupling model ----------------------------------------------------------
CouplingModelIndex = 1; %VALUE, popmenu selection%
[~,CouplingList]   = Coupling('List','None');

Output.CouplingType = CouplingList{CouplingModelIndex};
Output.Beta_NN      = 0.8; %NUMBER, cm^-1%
% ------------------------------------------------------------------------

% For Diagonal Disorder ---------------------------------------------------
Output.Sampling   =   0; %VALUE, toggles%
Output.Sample_Num = 100; %NUMBER%
Output.FWHM       =  10; %ARRAY, cm^-1%
Output.P_FlucCorr = 100; %NUMBER, percentage %
% ------------------------------------------------------------------------

% For 1D ------------------------------------------------------------------
Output.A_IR    = 90; %NUMBER, incident angle, Degrees%
Output.A_Vis1D = 90; %NUMBER, incident angle, Degrees%
Output.A_Sig1D = 90; %NUMBER, incident angle, Degrees%
Output.P_IR    =  0; %NUMBER, polarization angle 0=> P; 90=> S%
Output.P_Vis1D =  0; %NUMBER, polarization angle 0=> P; 90=> S%
Output.P_Sig1D =  0; %NUMBER, polarization angle 0=> P; 90=> S%
% ------------------------------------------------------------------------

% For 2D ------------------------------------------------------------------
Output.A_Pump  = 90; %NUMBER, incident angle, Degrees% 
Output.A_Probe = 90; %NUMBER, incident angle, Degrees% 
Output.A_Vis2D = 90; %NUMBER, incident angle, Degrees% 
Output.A_Sig2D = 90; %NUMBER, incident angle, Degrees% 
Output.P_Pump1 =  0; %NUMBER, polarization angle 0=> P; 90=> S% 
Output.P_Pump2 =  0; %NUMBER, polarization angle 0=> P; 90=> S% 
Output.P_Probe =  0; %NUMBER, polarization angle 0=> P; 90=> S% 
Output.P_Vis2D =  0; %NUMBER, polarization angle 0=> P; 90=> S% 
Output.P_Sig2D =  0; %NUMBER, polarization angle 0=> P; 90=> S% 
% ------------------------------------------------------------------------

% For MF to LF ------------------------------------------------------------
Output.Avg_Phi    =  0; %NUMBER, Angle, Degrees%
Output.Avg_Theta  =  0; %NUMBER, Angle, Degrees%  
Output.Avg_Psi    =  0; %NUMBER, Angle, Degrees%
Output.Avg_Rot    =  1; %VALUE, popmenu selection%
Output.Avg_Mirror =  1; %VALUE, popmenu selection%
% ------------------------------------------------------------------------



