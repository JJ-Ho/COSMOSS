function [O,T] = ParseGUI_Comb2(hGUI)

% For relative orientation ------------------------------------------------
O.Conc_Scaling =    str2double(get(hGUI.Conc_Scaling ,'String'));
O.Trans_X      =    str2double(get(hGUI.Trans_X      ,'String'));
O.Trans_Y      =    str2double(get(hGUI.Trans_Y      ,'String'));
O.Trans_Z      =    str2double(get(hGUI.Trans_Z      ,'String'));
O.Rot_Phi      =    str2double(get(hGUI.Rot_Phi      ,'String'));
O.Rot_Psi      =    str2double(get(hGUI.Rot_Psi      ,'String'));
O.Rot_Theta    =    str2double(get(hGUI.Rot_Theta    ,'String'));
% -------------------------------------------------------------------------

%% Export field name for each parameter
% For relative orientation ------------------------------------------------
T.Conc_Scaling = 'String';
T.Trans_X      = 'String';
T.Trans_Y      = 'String';
T.Trans_Z      = 'String';
T.Rot_Phi      = 'String';
T.Rot_Psi      = 'String';
T.Rot_Theta    = 'String';
% -------------------------------------------------------------------------
