function O = Parse_COSMOSS(app)
% Molecular response ------------------------------------------------------
O.LocFreqType = app.DropDown_LocFreq.Value;

% Coupling model
O.CouplingType = app.DropDown_Coupling.Value;
O.Beta_NN      = app.EditField_NN.Value;
O.PCutOff      = app.EditField_FeynmannCutoff.Value;
% -------------------------------------------------------------------------

% For Sample Symmetry -----------------------------------------------------
O.Avg_Rot    = app.DropDown_RotAvg.Value;
O.Avg_Mirror = app.DropDown_MirrorPlane.Value;
% -------------------------------------------------------------------------

% Esemble average ---------------------------------------------------------
O.Sampling      = app.CheckBox_Sampling.Value;
O.Sample_Num    = app.EditField_SampleN.Value;
O.DynamicUpdate = app.CheckBox_DynamicFigUpdate.Value;
O.UpdateStatus  = app.CheckBox_Continue.Value;
O.P_FlucCorr    = app.EditField_FluctuationCorrelation.Value;
O.DD_FWHM       = app.EditField_DD.Value;
O.ODD_FWHM      = app.EditField_ODD.Value;
% -------------------------------------------------------------------------

% Exp parameters 1D -------------------------------------------------------
Exp1D = app.UITable_1D.Data;
O.A_IR         = Exp1D(1,1);
O.A_Vis1D      = Exp1D(1,2);
O.A_Sig1D      = Exp1D(1,3);
O.P_IR         = Exp1D(2,1);
O.P_Vis1D      = Exp1D(2,2);
O.P_Sig1D      = Exp1D(2,3);
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
% -------------------------------------------------------------------------

% For Figures -------------------------------------------------------------
O.SaveFig      = app.CheckBox_SaveFig.Value;
O.SavePath     = app.EditField_SavePath.Value;
O.PlotStick    = app.CheckBox_Sticks.Value;
O.PlotNorm     = app.CheckBox_Normalize.Value;
O.PlotCursor   = app.CheckBox_Cursor.Value;
O.F_Min        = app.EditField_FreqMin.Value;
O.F_Max        = app.EditField_FreqMax.Value;
O.FreqRange    = O.F_Min:O.F_Max;
O.LineShape    = app.DropDown_LineShape.Value; % Need to debug
O.LineWidth    = app.EditField_Width.Value;
O.Num_Contour  = app.EditField_Contour.Value;
O.CMAP_Index   = app.DropDown_ColorMap.Value; % Need to debug
O.SpecType     = app.DropDown_SignalType.Value; % Need to debug
O.Signal_Type  = app.DropDown_SpectralType; % Need to debug
O.Pathway      = app.DropDown_Pathway.Value; % Need to debug
% -------------------------------------------------------------------------

% For Analysis Tools ------------------------------------------------------

% -------------------------------------------------------------------------




