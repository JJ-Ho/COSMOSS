function O = ParseGUI_Main(handles)
%% Get GUI inputs
GUI_Main = handles.GUI_Main;

O.Label_Index  =    str2num(get(GUI_Main.LIndex     ,'String'));
O.Label_Freq   = str2double(get(GUI_Main.LFreq      ,'String'));
O.Beta_NN      = str2double(get(GUI_Main.Beta_NN    ,'String'));
O.LineWidth    = str2double(get(GUI_Main.LineWidth  ,'String'));
O.Sampling     =            get(GUI_Main.Sampling   ,'Value');
O.Sample_Num   = str2double(get(GUI_Main.Sample_Num ,'String'));
O.FWHM         = str2double(get(GUI_Main.FWHM       ,'String')); 
O.P_FlucCorr   = str2double(get(GUI_Main.P_FlucCorr ,'String')); 

% For 1D ------------------------------------------------------------------
O.P_IR         = str2double(get(GUI_Main.P_IR       ,'String'));
O.P_Vis1D      = str2double(get(GUI_Main.P_Vis1D    ,'String'));
O.P_Sig1D      = str2double(get(GUI_Main.P_Sig1D    ,'String'));
O.A_IR         = str2double(get(GUI_Main.A_IR       ,'String')); % fix it with update JonesTensor later
O.A_Vis1D      = str2double(get(GUI_Main.A_Vis1D    ,'String'));
O.A_Sig1D      = str2double(get(GUI_Main.A_Sig1D    ,'String'));
% ------------------------------------------------------------------------

% For 2D ------------------------------------------------------------------
O.A_Pump       = str2double(get(GUI_Main.A_Pump1    ,'String')); % fix it with update JonesTensor later
O.A_Probe      = str2double(get(GUI_Main.A_Probe    ,'String'));
O.A_Vis2D      = str2double(get(GUI_Main.A_Vis2D    ,'String'));
O.A_Sig2D      = str2double(get(GUI_Main.A_Sig2D    ,'String'));
O.P_Pump1      = str2double(get(GUI_Main.P_Pump1    ,'String'));
O.P_Pump2      = str2double(get(GUI_Main.P_Pump2    ,'String'));
O.P_Probe      = str2double(get(GUI_Main.P_Probe    ,'String'));
O.P_Vis2D      = str2double(get(GUI_Main.P_Vis2D    ,'String'));
O.P_Sig2D      = str2double(get(GUI_Main.P_Sig2D    ,'String'));
% ------------------------------------------------------------------------
  
O.Avg_Rot      =            get(GUI_Main.Avg_Rot    ,'Value');
O.Avg_Mirror   =            get(GUI_Main.Avg_Mirror ,'Value');

O.Num_Contour  = str2double(get(GUI_Main.Num_Contour,'String'));
O.PlotStick    =            get(GUI_Main.PlotStick  ,'Value');
O.PlotNorm     =            get(GUI_Main.PlotNorm   ,'Value');
O.PlotCursor   =            get(GUI_Main.PlotCursor ,'Value');

%% Variable Need to be processed 
O.F_Min     = str2double(get(GUI_Main.X_Min,'String'));
O.F_Max     = str2double(get(GUI_Main.X_Max,'String'));
O.FreqRange = O.F_Min:O.F_Max;

Coupling = get(GUI_Main.Coupling,'Value');
switch Coupling
    case 1
        O.Coupling = 'NN_Mix_TDC';
    case 2
        O.Coupling = 'TDC';
    case 3
        O.Coupling = 'Cho_PB';
    case 4 
        O.Coupling = 'Cho_APB';
end

LineShape = get(GUI_Main.LineShape,'Value');
switch LineShape
    case 1
        O.LineShape = 'G';
    case 2
        O.LineShape = 'L';
    case 3
        O.LineShape = 'KK';
end

SpecType = get(GUI_Main.SpecType, 'Value');
switch SpecType
    case 1
        O.SpecType = 'Abs';
    case 2
        O.SpecType = 'R';
    case 3
        O.SpecType = 'NR';
end

Pathway = get(GUI_Main.Pathway,'Value');
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

O.Signal_Type = get(GUI_Main.Sig_Type,'Value');