% Template of Server_2DSFG running script
% This template scan the structural pameters on 2D_Grid of MBA 
%% setup path
% % Initialization Script
COSMOSS_path = '/home/jho/Projects/2DSFG_FGAIL_SAM/COSMOSS';
addpath(COSMOSS_path);

Initialization
Save_Output = 1;
Debug_FTIR  = 0;
Debug_SFG   = 0;
Debug_2DIR  = 0;
Debug_2DSFG = 0;

%% Input parameters
% APB Relative orientation
PHI    = #PHI#;
PSI    = #PSI#;
THETA  = #THETA#;
TransV = [0,0,#Z#];

% % APB Relative orientation
% PHI    = 0;
% PSI    = 0;
% THETA  = 0;
% TransV = [0,0,15];

% Name the outputs
Rot_Vec_Str   = sprintf('%03.0f_%03.0f_%03.0f',PHI,PSI,THETA);
Trans_Vec_Str = sprintf('%01.0f_%01.0f_%02.0f',TransV(1),TransV(2),TransV(3));
SaveName = ['APB_R5S3_MBA_4x4','_R_',Rot_Vec_Str,'_V_',Trans_Vec_Str];

% COSMOSS Inputs
COSMOSS_Input = Standard_Main_Input;
COSMOSS_Input.LocFreqType  = 2;% Use Jensen Locl mode freuqency type
COSMOSS_Input.CouplingType = 'Jansen_TDC';
COSMOSS_Input.Sampling     = 1;
COSMOSS_Input.Sample_Num   = 1000;
COSMOSS_Input.DD_FWHM      = 20; % cm-1
COSMOSS_Input.ODD_FWHM     = 5;  % cm-1
COSMOSS_Input.P_FlucCorr   = 100; % percentage
COSMOSS_Input.MEM_CutOff   = 1; % memmore cutoffGB
COSMOSS_Input.PCutOff      = 1E-5; % Pathway min intensity cutoff
COSMOSS_Input.F_Min        = 1500; % cm-1
COSMOSS_Input.F_Max        = 1900; % cm-1
COSMOSS_Input.FreqRange    = 1500:1900;
COSMOSS_Input.LineShape    = 'L';
COSMOSS_Input.LineWidth    = 5;
COSMOSS_Input.Pathway      = 'All';
COSMOSS_Input.SpecType     = 'Abs';

% Mimicing GUI input for ensemble average part
hGUIs.DynamicUpdate.Value = 0;
hGUIs.UpdateStatus.Value  = 1;

%% Construct the MBA 2D Grid
% Monomer
G09_Path = which('131029_MBA_Reverse_TDV.txt');
S_G09_Monomer = ReadG09Input(G09_Path);
S_G09_Monomer.hPlotFunc = @PlotXYZ_Grid; % Add function handles for plotting

% Grid 1
M1.MF_Center  = 13;
M1.MF_Zi      = 13;
M1.MF_Zf      = 3;
M1.MF_XZi     = 2;
M1.MF_XZf     = 12;
M1.Frame_Type = 1;
M1.LF_Phi     = -130;
M1.LF_Psi     =  270;  
M1.LF_Theta   = 37.4;  

S_M1 = R_MF2LF(S_G09_Monomer,M1);
S_M1.hPlotFunc = @PlotXYZ_Grid;% Add function handles for plotting

G.Vec_1 = 5.1*[1,0,0];
G1.Vec_2 = 3.2*[-sqrt(3),sqrt(13),0];
G1.N_1   = 4;
G1.N_2   = 2;
S_Grid1 = ConstructGrid(S_M1,G1);
S_Grid1.hPlotFunc = @PlotXYZ_Grid;% Add function handles for plotting

% Grid 2
M2 = M1;
M2.LF_Phi   =  -65; 
M2.LF_Theta = 37.4;  

S_M2 = R_MF2LF(S_G09_Monomer,M2);
S_M2.hPlotFunc = @PlotXYZ_Grid;% Add function handles for plotting

G2 = G1;
S_Grid2 = ConstructGrid(S_M2,G2);
S_Grid2.hPlotFunc = @PlotXYZ_Grid;% Add function handles for plotting

% Generate the full grid
Grid2TransV = 3.2*[-sqrt(3),sqrt(13),0]/2;

% move S_Grid1&2's origin to the first "S" atom
S_Grid1_0 = SD_Trans(S_Grid1, -S_Grid1.XYZ(13,:));
S_Grid2_0 = SD_Trans(S_Grid2, -S_Grid2.XYZ(13,:));
S_Grid2_T = SD_Trans(S_Grid2_0, Grid2TransV );

S_Grid_All = SD_Comb2(S_Grid1_0,S_Grid2_T);
S_Grid_All.hPlotFunc = @PlotXYZ_Grid;

%% Construct a ideal betasheet
BSheet.SheetTypeV = 2; % APB
BSheet.N_Residue  = 5;
BSheet.N_Strand   = 3;
BSheet.Trans_X    = 0;
BSheet.Trans_Y    = 0;
BSheet.Trans_Z    = 4.75;
BSheet.Twist_X    = 0;
BSheet.Twist_Y    = 0;
BSheet.Twist_Z    = 0;

S_BSheet = ConstuctBetaSheet(BSheet);
S_BSheet = SD_GetAmideI(S_BSheet);
S_BSheet.LocAnharm = ones(S_BSheet.Nmodes,1).* 20;
S_BSheet.LocFreq   = ones(S_BSheet.Nmodes,1).* 1665;

%% Combine the Betasheet with the 2D Grid
% MBA Grid
S_Grid_All_0 = SD_Trans(S_Grid_All,-S_Grid_All.CoM);
S_Grid_All_0.hPlotFunc = @PlotXYZ_Grid;

% Betasheet
PHI_r   = PHI/180*pi;
PSI_r   = PSI/180*pi;
THETA_r = THETA/180*pi;
R_Matrix = R1_ZYZ_0(PHI_r,PSI_r,THETA_r);
S_BSheet_0    = SD_Trans(S_BSheet,-S_BSheet.CoM);
S_BSheet_R    = SD_Rot(S_BSheet_0,R_Matrix);
S_BSheet_RT   = SD_Trans(S_BSheet_R,TransV);
S_BSheet_RT.hPlotFunc = @Plot_Betasheet_AmideI;
S_BSheet_RT.Extra.RotV = [PHI,PSI,THETA];

% Combine
S_APB_MBA = SD_Comb2(S_Grid_All_0,S_BSheet_RT);
S_APB_MBA.hPlotFunc = @PlotComb2;% Add function handles for plotting

%% Test run a 2D SFG to determine the scaling factor so that the MBA peak is about the same height of Betasheet peak
COSMOSS_Input_Test = COSMOSS_Input;
COSMOSS_Input_Test.Sampling = 0;
[SpectraGrid0,~] = TwoD_Iteration(@TwoDSFG_Main_Sparse,S_APB_MBA,COSMOSS_Input_Test,hGUIs);
CVL_Test = Conv2D(SpectraGrid0,COSMOSS_Input_Test);

if Debug_2DSFG
    hF  = figure;
    hAx = axes('Parent',hF);
    CVL_Test.FilesName = 'Scaling Test'; % pass filesname for figure title
    Plot2D(hAx,CVL_Test,COSMOSS_Input_Test,'2DSFG');
end

% Define Region og interest
Range_BSheet = (1600:1680) - COSMOSS_Input_Test.FreqRange(1) + 1;
Range_MBA    = (1710:1750) - COSMOSS_Input_Test.FreqRange(1) + 1;

M_All = real(CVL_Test.selected);
M_BSheet = M_All(Range_BSheet,Range_BSheet);
M_MBA    = M_All(   Range_MBA,   Range_MBA);

Max_BSheet = max(M_BSheet(:));
Max_MBA    = max(M_MBA(:));

ScalingF = (Max_MBA/Max_BSheet)^(1/4);
% Apply the scaling factor
S_BSheet_RT_S = SD_ScaleTransitions(S_BSheet_RT,ScalingF);
S_BSheet_RT_S.hPlotFunc  = @Plot_Betasheet_AmideI;% Add function handles for plotting
S_BSheet_RT_S.Extra.RotV = [PHI,PSI,THETA];

% Combine
S_APB_MBA_Scaled = SD_Comb2(S_Grid_All_0,S_BSheet_RT_S);
S_APB_MBA_Scaled.hPlotFunc = @PlotComb2;% Add function handles for plotting
Structure = S_APB_MBA_Scaled;

%% FTIR
FTIR = OneD_Iteration(@FTIR_Main,Structure,COSMOSS_Input,hGUIs);
if Debug_FTIR
    hF  = figure;
    hAx = axes('Parent',hF);
    Plot1D(hAx,FTIR,COSMOSS_Input);
end

%% SFG
OneDSFG = OneD_Iteration(@OneDSFG_Main,Structure,COSMOSS_Input,hGUIs);
if Debug_SFG
    hF  = figure;
    hAx = axes('Parent',hF);
    Plot1D(hAx,OneDSFG,COSMOSS_Input);
end

%% 2DIR
[SpectraGrid,Response] = TwoD_Iteration(@TwoDIR_Main_Sparse,Structure,COSMOSS_Input,hGUIs);
TwoDIR                 = Response;
TwoDIR.SpectraGrid     = SpectraGrid;

if Debug_2DIR
    hF  = figure;
    hAx = axes('Parent',hF);
    CVL = Conv2D(SpectraGrid,COSMOSS_Input);
    CVL.FilesName = Structure.FilesName; % pass filesname for figure title
    Plot2D(hAx,CVL,COSMOSS_Input,Response.SpecType);
end

%% 2D SFG
[SpectraGrid,Response] = TwoD_Iteration(@TwoDSFG_Main_Sparse,Structure,COSMOSS_Input,hGUIs);
TwoDSFG                = Response;
TwoDSFG.SpectraGrid    = SpectraGrid;

if Debug_2DSFG
    hF  = figure;
    hAx = axes('Parent',hF);
    CVL = Conv2D(SpectraGrid,COSMOSS_Input);
    CVL.FilesName = Structure.FilesName; % pass filesname for figure title
    Plot2D(hAx,CVL,COSMOSS_Input,Response.SpecType);
end

%% Outputs
Output.Structure = Structure;
Output.FTIR      = FTIR;
Output.OneDSFG   = OneDSFG;
Output.TwoDIR    = TwoDIR;
Output.TwoDSFG   = TwoDSFG;
Output.ScalingF  = ScalingF;
if Save_Output
    save(SaveName,'Output','SpectraGrid','Structure','ScalingF')
end

exit