function O = ParseGUI_Modes(handles)
%% Get GUI inputs
GUI_Modes = handles.GUI_Modes;

O.Plot_Loc       =            get(GUI_Modes.Plot_Loc      , 'Value' ) ; 
O.Loc_Ind        =    str2num(get(GUI_Modes.Loc_Ind       , 'String'));

O.Plot_Ex        =            get(GUI_Modes.Plot_Ex       , 'Value' ) ; 
O.Ex_Ind         =    str2num(get(GUI_Modes.Ex_Ind        , 'String'));

O.Plot_TDV       =            get(GUI_Modes.Plot_TDV      , 'Value' ) ;
O.Scale_TDV      = str2double(get(GUI_Modes.Scale_TDV     , 'String'));

O.Plot_Raman     =            get(GUI_Modes.Plot_Raman    , 'Value' ) ;
O.Scale_Raman    = str2double(get(GUI_Modes.Scale_Raman   , 'String'));

O.Normalize      =            get(GUI_Modes.Normalize     , 'Value' ) ;

O.Plot_EigenVec  =            get(GUI_Modes.Plot_EigenVec , 'Value' ) ;
O.EigneVec_Ind   = str2double(get(GUI_Modes.EigneVec_Ind  , 'String'));

