function [O,T] = ParseGUI_TwoDGrid(hGUIs)

% For Monomer -------------------------------------------------------------
O.Frame_Type   =               get(hGUIs.Frame_Type  ,'Value');
O.MF_Center    =       str2num(get(hGUIs.MF_Center   ,'String'));
O.MF_Zi        =    str2double(get(hGUIs.MF_Zi       ,'String'));
O.MF_Zf        =    str2double(get(hGUIs.MF_Zf       ,'String'));
O.MF_XZi       =    str2double(get(hGUIs.MF_XZi      ,'String'));
O.MF_XZf       =    str2double(get(hGUIs.MF_XZf      ,'String'));
O.BondAvg      =               get(hGUIs.BondAvg     ,'Value');

O.LF_Phi   =    str2double(get(hGUIs.LF_Phi         ,'String'));
O.LF_Psi   =    str2double(get(hGUIs.LF_Psi         ,'String'));
O.LF_Theta =    str2double(get(hGUIs.LF_Theta       ,'String'));
% -------------------------------------------------------------------------

% For 2D Grid -------------------------------------------------------------
O.Delta_Phi   =    str2double(get(hGUIs.Delta_Phi   ,'String'));
O.Delta_Psi   =    str2double(get(hGUIs.Delta_Psi   ,'String'));
O.Delta_Theta =    str2double(get(hGUIs.Delta_Theta ,'String'));
O.Vec_1       =       str2num(get(hGUIs.Vec_1 ,'String'));
O.Vec_2       =       str2num(get(hGUIs.Vec_2 ,'String'));
O.N_1         =    str2double(get(hGUIs.N_1   ,'String'));
O.N_2         =    str2double(get(hGUIs.N_2   ,'String'));
% -------------------------------------------------------------------------

% For Local Mode frequencies ----------------------------------------------
O.NLFreq  = str2double(get(hGUIs.NLFreq ,'String'));
O.LFreq   = str2double(get(hGUIs.LFreq  ,'String'));
O.Anharm  = str2double(get(hGUIs.Anharm ,'String'));
O.L_Index =    str2num(get(hGUIs.L_Index,'String'));
% -------------------------------------------------------------------------

% For molecule ploting ----------------------------------------------------
O.Plot_Atoms      = get(hGUIs.Plot_Atoms     ,'Value' );
O.Plot_Bonds      = get(hGUIs.Plot_Bonds     ,'Value' );
O.Plot_Axis       = get(hGUIs.Plot_Axis      ,'Value' );
O.Plot_Lattice    = get(hGUIs.Plot_Lattice   ,'Value' );
O.Plot_Atom_Index = get(hGUIs.Plot_Atom_Index,'Value' );
% -------------------------------------------------------------------------

%% Export field name for each parameter
% For Monomer -------------------------------------------------------------
T.Frame_Type   =  'Value'; 
T.MF_Center    = 'String';
T.MF_Zi        = 'String';
T.MF_Zf        = 'String';
T.MF_XYi       = 'String';
T.MF_XYf       = 'String';
T.BondAvg      =  'Value';
T.LF_Phi       = 'String';
T.LF_Psi       = 'String';
T.LF_Theta     = 'String';
% -------------------------------------------------------------------------

% For 2D Grid -------------------------------------------------------------
T.Delta_Phi   = 'String';
T.Delta_Psi   = 'String';
T.Delta_Theta = 'String';
T.Vec_1       = 'String';
T.Vec_2       = 'String';
T.N_1         = 'String';
T.N_2         = 'String';
% -------------------------------------------------------------------------

% For Local Mode frequencies ----------------------------------------------
T.NLFreq  = 'String';
T.LFreq   = 'String';
T.Anharm  = 'String';
T.L_Index = 'String';
% -------------------------------------------------------------------------

% For molecule ploting ----------------------------------------------------
T.Plot_Atoms      = 'Value';
T.Plot_Bonds      = 'Value';
T.Plot_Axis       = 'Value';
T.Plot_Lattice    = 'Value';
T.Plot_Atom_Index = 'Value';
% -------------------------------------------------------------------------

