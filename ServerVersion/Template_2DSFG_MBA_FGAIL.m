% Template of Server_2DSFG running script
% This template scan the vertical distance between MBA and R5S3 betasheet

%% Initialization Script
% Initialization Script
COSMOSS_path = '/home/jho/Projects/2DSFG_FGAIL_SAM/COSMOSS';
addpath(COSMOSS_path);

% Add path to SubFunctions
SpectalFunctions = genpath([COSMOSS_path '/SpecFunc']);
addpath(SpectalFunctions);

MoleculeConstruction = genpath([COSMOSS_path '/MoleculeConstruction']);
addpath(MoleculeConstruction);

StructureFiles = genpath([COSMOSS_path '/StructureFiles']);
addpath(StructureFiles);

ModeVisualization = genpath([COSMOSS_path '/ModeVisualization']);
addpath(ModeVisualization);

Main_GUIs = genpath([COSMOSS_path '/Main_GUI']);
addpath(Main_GUIs);

ServerVersion = genpath([COSMOSS_path '/ServerVersion']);
addpath(ServerVersion);


%% Manual Inputs 
% Relative orientaiton
Conc_Scaling = 1;
Trans_X      = 10.5;
Trans_Y      = #SCAN_Y#; % scan from -3:1:3
Trans_Z      = #SCAN_Z#; % Scan from 10:1:20
Rot_Phi      = 0;
Rot_Psi      = 0;
Rot_Theta    = 0;

% Main part
Coupling   = 'TDC+Cho_APB';
Sampling   = 1;
Sample_Num = 200;
FWHM       = 30;

%% Load default Inputs and Change default Inputs with Manual Inputs
Main_Input     = Standard_Main_Input;

Main_Input.Coupling   = Coupling;
Main_Input.Sampling   = Sampling;
Main_Input.Sample_Num = Sample_Num;
Main_Input.FWHM       = FWHM;

%% Structural loading and construction
Tmp1 = load('MBA_XPS_Grid_5x4');
StrucData1 = Tmp1.Structure;

Tmp2 = load('APB_R5S3');
StrucData2 = Tmp2.Structure;

GUI_Inputs.Conc_Scaling = Conc_Scaling;
GUI_Inputs.Trans_X      = Trans_X;
GUI_Inputs.Trans_Y      = Trans_Y;
GUI_Inputs.Trans_Z      = Trans_Z;
GUI_Inputs.Rot_Phi      = Rot_Phi/180*pi;
GUI_Inputs.Rot_Psi      = Rot_Psi/180*pi;
GUI_Inputs.Rot_Theta    = Rot_Theta/180*pi;

Structure = Comb2(StrucData1,StrucData2,GUI_Inputs);

%% For parfor version used
matlabpool close force local
matlabpool open local 8

%% Call Server_2DSFG
[SpectraGrid,Response] = Server_2DSFG(Structure,Main_Input);
%% Ouputs
Output.SpectraGrid = SpectraGrid;
Output.Response    = Response;

SaveName = ['MBA_FGAIL_Y_',num2str(Trans_Y),'_Z_',num2str(Trans_Z)];
save(SaveName,'Output')

matlabpool close


