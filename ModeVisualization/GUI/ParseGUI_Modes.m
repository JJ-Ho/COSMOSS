function O = ParseGUI_Modes(hGUIs)
%% Get GUI inputs
O.SpecType       =            get(hGUIs.SpecType      , 'Value' ) ;

O.Mu_Alpha_Plot  =            get(hGUIs.Mu_Alpha_Plot , 'Value' ) ; 
O.Mu_Alpha_Type  =            get(hGUIs.Mu_Alpha_Type , 'Value' ) ; 
O.Mu_Alpha_Ind   =    str2num(get(hGUIs.Mu_Alpha_Ind  , 'String'));

O.TDV_Plot       =            get(hGUIs.TDV_Plot      , 'Value' ) ;
O.TDV_Scale      = str2double(get(hGUIs.TDV_Scale     , 'String'));
O.Raman_Plot     =            get(hGUIs.Raman_Plot    , 'Value' ) ;
O.Raman_Scale    = str2double(get(hGUIs.Raman_Scale   , 'String'));
O.Raman_Type     =            get(hGUIs.Raman_Type    , 'Value' ) ;

O.Normalize      =            get(hGUIs.Normalize     , 'Value' ) ;

O.Plot_EigenVec  =            get(hGUIs.Plot_EigenVec , 'Value' ) ;
O.EigneVec_Ind   = str2double(get(hGUIs.EigneVec_Ind  , 'String'));
O.EigneVec_Conv  =            get(hGUIs.EigneVec_Conv , 'Value' ) ;
