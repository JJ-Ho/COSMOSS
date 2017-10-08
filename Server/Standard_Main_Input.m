function O = Standard_Main_Input
% For Sample Symmetry -----------------------------------------------------
O.Avg_Rot      = 1; %'Phi' C_Inf
O.Avg_Mirror   = 1; % no mirror plane
% -------------------------------------------------------------------------

% For 1D ------------------------------------------------------------------
O.A_IR    = 90; %NUMBER, incident angle, Degrees%
O.A_Vis1D = 90; %NUMBER, incident angle, Degrees%
O.A_Sig1D = 90; %NUMBER, incident angle, Degrees%
O.P_IR    =  0; %NUMBER, polarization angle 0=> P; 90=> S%
O.P_Vis1D =  0; %NUMBER, polarization angle 0=> P; 90=> S%
O.P_Sig1D =  0; %NUMBER, polarization angle 0=> P; 90=> S%
% -------------------------------------------------------------------------

% For 2D ------------------------------------------------------------------
O.A_Pump  = 90; %NUMBER, incident angle, Degrees% 
O.A_Probe = 90; %NUMBER, incident angle, Degrees% 
O.A_Vis2D = 90; %NUMBER, incident angle, Degrees% 
O.A_Sig2D = 90; %NUMBER, incident angle, Degrees% 
O.P_Pump1 =  0; %NUMBER, polarization angle 0=> P; 90=> S% 
O.P_Pump2 =  0; %NUMBER, polarization angle 0=> P; 90=> S% 
O.P_Probe =  0; %NUMBER, polarization angle 0=> P; 90=> S% 
O.P_Vis2D =  0; %NUMBER, polarization angle 0=> P; 90=> S% 
O.P_Sig2D =  0; %NUMBER, polarization angle 0=> P; 90=> S% 
% -------------------------------------------------------------------------

% For Diagonal Disorder ---------------------------------------------------
O.Sampling      =   0; %VALUE, toggles%
O.Sample_Num    = 100; %NUMBER%
O.DD_FWHM       = 10;
O.ODD_FWHM      = 5;
O.FWHM          =  10; %ARRAY, cm^-1%
O.P_FlucCorr    = 100; %NUMBER, percentage %
O.DynamicUpdate = 0;
O.UpdateStatus  = 0;
% -------------------------------------------------------------------------

% For Coupling/Signal -----------------------------------------------------
O.LocFreqType   = 1;

% Coupling model
CouplingModelIndex = 3;
[~,CouplingList]   = Coupling('List','None','None');
O.CouplingType     = CouplingList{CouplingModelIndex};

O.Beta_NN       = 0.8;
O.PCutOff       = 1E-5;
% -------------------------------------------------------------------------


% Frequency Range ---------------------------------------------------------
O.F_Min     = 1550; %NUMBER, cm^-1%
O.F_Max     = 1800; %NUMBER, cm^-1%
O.FreqRange = O.F_Min:O.F_Max;
% ------------------------------------------------------------------------




