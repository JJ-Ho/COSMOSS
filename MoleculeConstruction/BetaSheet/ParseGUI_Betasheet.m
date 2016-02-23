function O = ParseGUI_Betasheet(hGUI)

% Orientation -------------------------------------------------------------
O.SheetType  =            get(hGUI.SheetType ,'Value' ) ;
O.N_Residue  = str2double(get(hGUI.N_Residue ,'String'));
O.N_Strand   = str2double(get(hGUI.N_Strand  ,'String'));
O.Trans_X    = str2double(get(hGUI.Trans_X   ,'String'));
O.Trans_Y    = str2double(get(hGUI.Trans_Y   ,'String'));
O.Trans_Z    = str2double(get(hGUI.Trans_Z   ,'String'));
O.Twist_X    = str2double(get(hGUI.Twist_X   ,'String'));
O.Twist_Y    = str2double(get(hGUI.Twist_Y   ,'String'));
O.Twist_Z    = str2double(get(hGUI.Twist_Z   ,'String'));
O.Phi_D      = str2double(get(hGUI.Phi_D     ,'String'));
O.Psi_D      = str2double(get(hGUI.Psi_D     ,'String'));
O.Theta_D    = str2double(get(hGUI.Theta_D   ,'String'));

% For Local Mode frequencies ----------------------------------------------
O.NLFreq  = str2double(get(hGUI.NLFreq ,'String'));
O.LFreq   = str2double(get(hGUI.LFreq  ,'String'));
O.Anharm  = str2double(get(hGUI.Anharm ,'String'));
O.L_Index =    str2num(get(hGUI.L_Index,'String'));

% For molecule ploting ----------------------------------------------------
O.Plot_Atoms  = get(hGUI.Plot_Atoms ,'Value' );
O.Plot_Bonds  = get(hGUI.Plot_Bonds ,'Value' );
