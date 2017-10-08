% Template of Server_2DSFG running script
% This template scan the structural pameters on 2D_Grid of MBA 
%% Input parameters
% setup path
Initialization

% APB Relative orientation
PHI    = 0;
PSI    = 0;
THETA  = 0;
TransV = [0,0,15];

% Name the outputs
Rot_Vec_Str   = sprintf('%03.0f_%03.0f_%03.0f',PHI,PSI,THETA);
Trans_Vec_Str = sprintf('%01.0f_%01.0f_%02.0f',TransV(1),TransV(2),TransV(3));
SaveName = ['APB_R5S3_MBA_4x4','_R_',Rot_Vec_Str,'_V_',Trans_Vec_Str];

% COSMOSS Inputs
COSMOSS_Input = Standard_Main_Input;
COSMOSS_Input.Sampling = 1;
COSMOSS_Input.Sample_Num = 500;
COSMOSS_Input.MEM_CutOff = 1; % GB

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
R_Matrix = R1_ZYZ_0(PHI,PSI,THETA);
S_BSheet_0  = SD_Trans(S_BSheet,-S_BSheet.CoM);
S_BSheet_R  = SD_Rot(S_BSheet_0,R_Matrix);
S_BSheet_RT = SD_Trans(S_BSheet_R,TransV);
S_BSheet_RT.hPlotFunc = @Plot_Betasheet_AmideI;% Add function handles for plotting
S_BSheet_RT.Extra.RotV = [PHI,PSI,THETA];

% Combine
S_APB_MBA = SD_Comb2(S_Grid_All_0,S_BSheet_RT);
S_APB_MBA.hPlotFunc = @PlotComb2;% Add function handles for plotting
Structure = S_APB_MBA;

%% FTIR
FTIR = OneD_Iteration(@FTIR_Main,Structure,COSMOSS_Input,hGUIs);
% hF  = figure;
% hAx = axes('Parent',hF);
% Plot1D(hAx,FTIR,COSMOSS_Input);

%% SFG
OneDSFG = OneD_Iteration(@OneDSFG_Main,Structure,COSMOSS_Input,hGUIs);
% hF  = figure;
% hAx = axes('Parent',hF);
% Plot1D(hAx,SFG,COSMOSS_Input);

%% 2DIR
[SpectraGrid,Response] = TwoD_Iteration(@TwoDIR_Main_Sparse,Structure,COSMOSS_Input,hGUIs);
TwoDIR                 = Response;
TwoDIR.SpectraGrid     = SpectraGrid;

% hF  = figure;
% hAx = axes('Parent',hF);
% CVL = Conv2D(SpectraGrid,COSMOSS_Input);
% CVL.FilesName = Structure.FilesName; % pass filesname for figure title
% Plot2D(hAx,CVL,COSMOSS_Input,Response.SpecType);

%% 2D SFG
[SpectraGrid,Response] = TwoD_Iteration(@TwoDSFG_Main_Sparse,Structure,COSMOSS_Input,hGUIs);
TwoDSFG                = Response;
TwoDSFG.SpectraGrid    = SpectraGrid;

% hF  = figure;
% hAx = axes('Parent',hF);
% CVL = Conv2D(SpectraGrid,COSMOSS_Input);
% CVL.FilesName = Structure.FilesName; % pass filesname for figure title
% Plot2D(hAx,CVL,COSMOSS_Input,Response.SpecType);

%% Outputs
Output.Structure = Structure;
Output.FTIR      = FTIR;
Output.OneDSFG   = OneDSFG;
Output.TwoDIR    = TwoDIR;
Output.TwoDSFG   = TwoDSFG;
save(SaveName,'Output')