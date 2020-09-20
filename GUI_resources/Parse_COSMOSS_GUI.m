function O = Parse_COSMOSS_GUI(app)
% Molecular response ------------------------------------------------------
O.PCutOff    = app.EditField_FeynmannCutoff.Value/100; % percentage to ratio
% -------------------------------------------------------------------------

% For Sample Symmetry -----------------------------------------------------
O.Avg_Rot    = app.DropDown_RotAvg.Value;
O.Avg_Mirror = app.DropDown_MirrorPlane.Value;
% -------------------------------------------------------------------------

% Esemble average ---------------------------------------------------------
O.Sampling         = app.CheckBox_Sampling.Value;
O.Sample_Num       = app.EditField_SampleN.Value;
O.DynamicUpdate    = app.CheckBox_DynamicFigUpdate.Value;
O.UpdateStatus     = app.CheckBox_Continue.Value;
O.ViewSampling     = app.CheckBox_ViewSampling.Value;
O.ViewSamplingMode = app.DropDown_ViewSamplingMode.Value;
O.isotopeDilution  = app.isotopeDilutionCheckBox.Value;
O.P_FlucCorr       = app.EditField_FluctuationCorrelation.Value;
O.DD_FWHM          = app.EditField_DD.Value;
O.ODD_FWHM         = app.EditField_ODD.Value;
% -------------------------------------------------------------------------

% Exp parameters 1D -------------------------------------------------------
Exp1D = app.UITable_1D.Data;
O.A_IR         = Exp1D(1,1);
O.A_Vis1D      = Exp1D(1,2);
O.A_Sig1D      = Exp1D(1,3);
O.P_IR         = Exp1D(2,1);
O.P_Vis1D      = Exp1D(2,2);
O.P_Sig1D      = Exp1D(2,3);
O.Exp1D        = Exp1D;
% -------------------------------------------------------------------------

% Exp parameters 2D -------------------------------------------------------
Exp2D = app.UITable_2D.Data;
O.A_Pump1      = Exp2D(1,1);
O.A_Pump2      = Exp2D(1,2);
O.A_Probe      = Exp2D(1,3);
O.A_Vis2D      = Exp2D(1,4);
O.A_Sig2D      = Exp2D(1,5);
O.P_Pump1      = Exp2D(2,1);
O.P_Pump2      = Exp2D(2,2);
O.P_Probe      = Exp2D(2,3);
O.P_Vis2D      = Exp2D(2,4);
O.P_Sig2D      = Exp2D(2,5);
O.Exp2D        = Exp2D;
% -------------------------------------------------------------------------

% For Figures -------------------------------------------------------------
O.SaveFig      = app.CheckBox_SaveFig.Value;
O.SavePath     = app.EditField_SavePath.Value;
O.existFig     = app.ExisitingfigureCheckBox.Value;
O.hFig         = app.hFEditField.Value;

% 1D
O.PlotStick_1D    = app.CheckBox_Sticks_1D.Value;
O.PlotNorm_1D     = app.CheckBox_Normalize_1D.Value;
O.PlotCursor_1D   = app.CheckBox_Cursor_1D.Value;
O.IntArea_1D      = app.CheckBox_IntArea_1D.Value;
O.YSym_1D         = app.CheckBox_YSym_1D.Value;
O.F_Min_1D        = app.EditField_FreqMin_1D.Value;
O.F_Max_1D        = app.EditField_FreqMax_1D.Value;
O.FreqRange_1D    = O.F_Min_1D:O.F_Max_1D;
O.LineShape_1D    = app.DropDown_LineShape_1D.Value;
O.LineWidth_1D    = app.EditField_Width_1D.Value;
O.Signal_Type_1D  = app.DropDown_SignalType_1D.Value;

% 2D
O.PlotNorm_2D     = app.CheckBox_Normalize_2D.Value;
O.PlotCursor_2D   = app.CheckBox_Cursor_2D.Value;
O.Num_Contour_2D  = app.EditField_Contour_2D.Value;
O.F_Min_2D        = app.EditField_FreqMin_2D.Value;
O.F_Max_2D        = app.EditField_FreqMax_2D.Value;
O.FreqRange_2D    = O.F_Min_2D:O.F_Max_2D;
O.LineShape_2D    = app.DropDown_LineShape_2D.Value;
O.LineWidth_2D    = app.EditField_Width_2D.Value;
O.CMAP_Index_2D   = app.DropDown_ColorMap_2D.Value;
O.Pathway_2D      = app.DropDown_Pathway_2D.Value;
O.SpecType_2D     = app.DropDown_SpectralType_2D.Value;
%O.Signal_Type_2D  = app.DropDown_SignalType_2D.Value; % is not used in plot2D
% -------------------------------------------------------------------------

% For Analysis Tools ------------------------------------------------------

% -------------------------------------------------------------------------
