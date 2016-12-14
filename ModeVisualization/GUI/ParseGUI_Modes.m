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

O.Plot_EigVec    =            get(hGUIs.Plot_EigVec   , 'Value' ) ;
O.EigVec_Ind     =    str2num(get(hGUIs.EigVec_Ind    , 'String'));
O.EigVec_Conv    =            get(hGUIs.EigVec_Conv   , 'Value' ) ;
O.EigVec_Scale   = str2double(get(hGUIs.EigVec_Scale  , 'String'));

O.Sig_Scale      = str2double(get(hGUIs.Sig_Scale     , 'String'));
O.Sig_NGrid      = str2double(get(hGUIs.Sig_NGrid     , 'String'));
O.Sig_Plot3D     =            get(hGUIs.Sig_Plot3D    , 'Value' ) ;
O.Sig_PlotCT     =            get(hGUIs.Sig_PlotCT    , 'Value' ) ;
O.Sig_PlotCT_R   =            get(hGUIs.Sig_PlotCT_R  , 'Value' ) ;