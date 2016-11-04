function O = ParseGUI_TwoDGrid(hGUI)

% For Monomer -------------------------------------------------------------
O.Frame_Type   =               get(hGUI.Frame_Type  ,'Value');
O.MF_Center    =       str2num(get(hGUI.MF_Center   ,'String'));
O.MF_Zi        =    str2double(get(hGUI.MF_Zi       ,'String'));
O.MF_Zf        =    str2double(get(hGUI.MF_Zf       ,'String'));
O.MF_XYi       =    str2double(get(hGUI.MF_XYi      ,'String'));
O.MF_XYf       =    str2double(get(hGUI.MF_XYf      ,'String'));
O.BondAvg      =               get(hGUI.BondAvg     ,'Value');

O.LF_Phi   =    str2double(get(hGUI.LF_Phi         ,'String'));
O.LF_Psi   =    str2double(get(hGUI.LF_Psi         ,'String'));
O.LF_Theta =    str2double(get(hGUI.LF_Theta       ,'String'));
% -------------------------------------------------------------------------

% For 2D Grid -------------------------------------------------------------
O.Delta_Phi   =    str2double(get(hGUI.Delta_Phi   ,'String'));
O.Delta_Psi   =    str2double(get(hGUI.Delta_Psi   ,'String'));
O.Delta_Theta =    str2double(get(hGUI.Delta_Theta ,'String'));
O.Vec_1       =       str2num(get(hGUI.Vec_1 ,'String'));
O.Vec_2       =       str2num(get(hGUI.Vec_2 ,'String'));
O.N_1         =    str2double(get(hGUI.N_1   ,'String'));
O.N_2         =    str2double(get(hGUI.N_2   ,'String'));
% -------------------------------------------------------------------------

% For Local Mode frequencies ----------------------------------------------
O.NLFreq  = str2double(get(hGUI.NLFreq ,'String'));
O.LFreq   = str2double(get(hGUI.LFreq  ,'String'));
O.Anharm  = str2double(get(hGUI.Anharm ,'String'));
O.L_Index =    str2num(get(hGUI.L_Index,'String'));
% -------------------------------------------------------------------------

% For molecule ploting ----------------------------------------------------
O.Plot_Atoms      = get(hGUI.Plot_Atoms     ,'Value' );
O.Plot_Bonds      = get(hGUI.Plot_Bonds     ,'Value' );
O.Plot_Axis       = get(hGUI.Plot_Axis      ,'Value' );
O.Plot_Lattice    = get(hGUI.Plot_Lattice   ,'Value' );
O.Plot_Atom_Index = get(hGUI.Plot_Atom_Index,'Value' );
% -------------------------------------------------------------------------