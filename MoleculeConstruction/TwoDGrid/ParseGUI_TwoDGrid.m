function O = ParseGUI_TwoDGrid(hGUI)

% For Monomer -------------------------------------------------------------
O.Ang_Phi     =    str2double(get(hGUI.Ang_Phi     ,'String'));
O.Ang_Psi     =    str2double(get(hGUI.Ang_Psi     ,'String'));
O.Ang_Theta   =    str2double(get(hGUI.Ang_Theta   ,'String'));
O.Delta_Phi   =    str2double(get(hGUI.Delta_Phi   ,'String'));
O.Delta_Psi   =    str2double(get(hGUI.Delta_Psi   ,'String'));
O.Delta_Theta =    str2double(get(hGUI.Delta_Theta ,'String'));
% -------------------------------------------------------------------------

% For 2D Grid -------------------------------------------------------------
O.Vec_1 =    str2num(get(hGUI.Vec_1 ,'String'));
O.Vec_2 =    str2num(get(hGUI.Vec_2 ,'String'));
O.N_1   = str2double(get(hGUI.N_1   ,'String'));
O.N_2   = str2double(get(hGUI.N_2   ,'String'));
% -------------------------------------------------------------------------

% For Local Mode frequencies ----------------------------------------------
O.NLFreq  = str2double(get(hGUI.NLFreq ,'String'));
O.LFreq   = str2double(get(hGUI.LFreq  ,'String'));
O.Anharm  = str2double(get(hGUI.Anharm ,'String'));
O.L_Index =    str2num(get(hGUI.L_Index,'String'));