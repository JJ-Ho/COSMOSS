% Template of Server_2DSFG running script
% This template scan the structural pameters on 2D_Grid of MBA 

%% Construct 2D Grid
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
S_Grid_All.hPlotFunc = @PlotComb2;% Add function handles for plotting

%% Construct ideal betasheet


%% COSMOSS Inputs
COSMOSS_Input = Standard_Main_Input;

%% Call Server_2DSFG
hF  = figure;
hAx = axes('Parent',hF);

FTIR = FTIR_Main(S_Grid_All,COSMOSS_Input);
Plot1D(hAx,FTIR,COSMOSS_Input);

% [SpectraGrid,Response] = Server_2DSFG(Structure,COSMOSS_Input);


%% Ouputs
save(SaveName,'Output')