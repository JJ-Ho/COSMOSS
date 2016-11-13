function O = ParseGUI_Main(hGUIs)
%% Get GUI inputs

O.Label_Index   =    str2num(get(hGUIs.LIndex       ,'String'));
O.Label_Freq    = str2double(get(hGUIs.LFreq        ,'String'));
O.Beta_NN       = str2double(get(hGUIs.Beta_NN      ,'String'));
O.LineWidth     = str2double(get(hGUIs.LineWidth    ,'String'));
O.Sampling      =            get(hGUIs.Sampling     ,'Value');
O.Sample_Num    = str2double(get(hGUIs.Sample_Num   ,'String'));
O.FWHM          =    str2num(get(hGUIs.FWHM         ,'String')); 
O.P_FlucCorr    = str2double(get(hGUIs.P_FlucCorr   ,'String')); 
O.DynamicUpdate =            get(hGUIs.DynamicUpdate,'Value');
O.UpdateStatus  =            get(hGUIs.UpdateStatus ,'Value');

% For 1D ------------------------------------------------------------------
O.P_IR         = str2double(get(hGUIs.P_IR       ,'String'));
O.P_Vis1D      = str2double(get(hGUIs.P_Vis1D    ,'String'));
O.P_Sig1D      = str2double(get(hGUIs.P_Sig1D    ,'String'));
O.A_IR         = str2double(get(hGUIs.A_IR       ,'String')); % fix it with update JonesTensor later
O.A_Vis1D      = str2double(get(hGUIs.A_Vis1D    ,'String'));
O.A_Sig1D      = str2double(get(hGUIs.A_Sig1D    ,'String'));
% ------------------------------------------------------------------------

% For 2D ------------------------------------------------------------------
O.A_Pump       = str2double(get(hGUIs.A_Pump1    ,'String')); % fix it with update JonesTensor later
O.A_Probe      = str2double(get(hGUIs.A_Probe    ,'String'));
O.A_Vis2D      = str2double(get(hGUIs.A_Vis2D    ,'String'));
O.A_Sig2D      = str2double(get(hGUIs.A_Sig2D    ,'String'));
O.P_Pump1      = str2double(get(hGUIs.P_Pump1    ,'String'));
O.P_Pump2      = str2double(get(hGUIs.P_Pump2    ,'String'));
O.P_Probe      = str2double(get(hGUIs.P_Probe    ,'String'));
O.P_Vis2D      = str2double(get(hGUIs.P_Vis2D    ,'String'));
O.P_Sig2D      = str2double(get(hGUIs.P_Sig2D    ,'String'));
% ------------------------------------------------------------------------

% For MF to LF ------------------------------------------------------------
O.Avg_Phi      = str2double(get(hGUIs.Avg_Phi    ,'String'));
O.Avg_Theta    = str2double(get(hGUIs.Avg_Theta  ,'String'));
O.Avg_Psi      = str2double(get(hGUIs.Avg_Psi    ,'String'));
O.Avg_Rot      =            get(hGUIs.Avg_Rot    ,'Value');
O.Avg_Mirror   =            get(hGUIs.Avg_Mirror ,'Value');
% ------------------------------------------------------------------------

% For Figures ------------------------------------------------------------
O.Num_Contour  = str2double(get(hGUIs.Num_Contour ,'String'));
O.PlotStick    =            get(hGUIs.PlotStick   ,'Value');
O.PlotNorm     =            get(hGUIs.PlotNorm    ,'Value');
O.PlotCursor   =            get(hGUIs.PlotCursor  ,'Value');
O.IntegralArea =            get(hGUIs.IntegralArea,'Value');
O.CMAP_Index   =            get(hGUIs.CMAP_Index  ,'Value');
% ------------------------------------------------------------------------

%% Variable Need to be processed 
O.F_Min     = str2double(get(hGUIs.X_Min,'String'));
O.F_Max     = str2double(get(hGUIs.X_Max,'String'));
O.FreqRange = O.F_Min:O.F_Max;

% Coupling model
CouplingModelIndex = get(hGUIs.CouplingModelIndex,'Value');
[~,CouplingList]   = Coupling('List','None');
O.CouplingType     = CouplingList{CouplingModelIndex};

LineShape = get(hGUIs.LineShape,'Value');
switch LineShape
    case 1
        O.LineShape = 'G';
    case 2
        O.LineShape = 'L';
    case 3
        O.LineShape = 'KK';
    case 4
        O.LineShape = 'None'; 
end

SpecType = get(hGUIs.SpecType, 'Value');
switch SpecType
    case 1
        O.SpecType = 'Abs';
    case 2
        O.SpecType = 'R';
    case 3
        O.SpecType = 'NR';
end

Pathway = get(hGUIs.Pathway,'Value');
switch Pathway
    case 1
        O.Pathway = 'All';
    case 2
        O.Pathway = 'GB';
    case 3
        O.Pathway = 'SE';
    case 4
        O.Pathway = 'EA';
end

Signal_Type = get(hGUIs.Sig_Type,'Value');
switch Signal_Type
    case 1 % heterodyne
        O.Signal_Type = 'Hetero';
    case 2 % homodyne
        O.Signal_Type = 'Homo';
end