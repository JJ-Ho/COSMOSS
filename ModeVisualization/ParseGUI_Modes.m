function O = ParseGUI_Modes(handles)
%% Get GUI inputs
GUI_Modes = handles.GUI_Modes;

O.Mu_Alpha_Plot  =            get(GUI_Modes.Mu_Alpha_Plot , 'Value' ) ; 
O.Mu_Alpha_Type  =            get(GUI_Modes.Mu_Alpha_Type , 'Value' ) ; 
O.Mu_Alpha_Ind   =    str2num(get(GUI_Modes.Mu_Alpha_Ind  , 'String'));

O.TDV_Plot       =            get(GUI_Modes.TDV_Plot      , 'Value' ) ;
O.TDV_Scale      = str2double(get(GUI_Modes.TDV_Scale     , 'String'));
O.Raman_Plot     =            get(GUI_Modes.Raman_Plot    , 'Value' ) ;
O.Raman_Scale    = str2double(get(GUI_Modes.Raman_Scale   , 'String'));
O.Raman_Type     =            get(GUI_Modes.Raman_Type    , 'Value' ) ;

O.Normalize      =            get(GUI_Modes.Normalize     , 'Value' ) ;

O.Plot_EigenVec  =            get(GUI_Modes.Plot_EigenVec , 'Value' ) ;
O.EigneVec_Ind   = str2double(get(GUI_Modes.EigneVec_Ind  , 'String'));

