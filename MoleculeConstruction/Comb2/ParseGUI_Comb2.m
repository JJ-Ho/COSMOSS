function O = ParseGUI_Comb2(hGUI)

% For relative orientation ------------------------------------------------
O.Trans_X     =    str2double(get(hGUI.Trans_X     ,'String'));
O.Trans_Y     =    str2double(get(hGUI.Trans_Y     ,'String'));
O.Trans_Z     =    str2double(get(hGUI.Trans_Z     ,'String'));
O.Rot_Phi     =    str2double(get(hGUI.Rot_Phi     ,'String'));
O.Rot_Psi     =    str2double(get(hGUI.Rot_Psi     ,'String'));
O.Rot_Theta   =    str2double(get(hGUI.Rot_Theta   ,'String'));
% -------------------------------------------------------------------------
