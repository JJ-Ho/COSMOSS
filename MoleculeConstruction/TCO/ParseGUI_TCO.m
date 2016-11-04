function O = ParseGUI_TCO(hGUI)

O.Phi_D1        = str2double(get(hGUI.Phi1  ,'String'));
O.Psi_D1        = str2double(get(hGUI.Psi1  ,'String'));
O.Theta_D1      = str2double(get(hGUI.Theta1,'String'));
O.Phi_D2        = str2double(get(hGUI.Phi2  ,'String'));
O.Psi_D2        = str2double(get(hGUI.Psi2  ,'String'));
O.Theta_D2      = str2double(get(hGUI.Theta2,'String'));
O.Displacement  =    str2num(get(hGUI.Trans ,'String'));

O.Rot_X        = str2double(get(hGUI.Rot_X  ,'String'));
O.Rot_Y        = str2double(get(hGUI.Rot_Y  ,'String'));
O.Rot_Z        = str2double(get(hGUI.Rot_Z  ,'String'));

% For Local Mode frequencies ----------------------------------------------
O.NLFreq  = str2double(get(hGUI.NLFreq ,'String'));
O.LFreq   = str2double(get(hGUI.LFreq  ,'String'));
O.Anharm  = str2double(get(hGUI.Anharm ,'String'));
O.L_Index =    str2num(get(hGUI.L_Index,'String'));

% For molecule ploting ----------------------------------------------------
O.Plot_Atoms  = get(hGUI.Plot_Atoms ,'Value' );
O.Plot_Bonds  = get(hGUI.Plot_Bonds ,'Value' );
O.Plot_Axis   = get(hGUI.Plot_Axis  ,'Value' );
