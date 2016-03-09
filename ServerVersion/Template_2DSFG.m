% Template of Server_2DSFG running script
% This template scan the structural pameters on 2D_Grid of MBA 

%% Manual Inputs 
% Main part
Avg_Theta  = 30; % MF to LF in Degrees
Coupling   = 'TDC';
Sampling   = 1;
Sample_Num = 100;

% Structural part
G09_FileName = '131029_MBA.txt';
D1 = 12; % Grid vector length 1, Ang
D2 = 12; % Grid vector length 2, Ang

MBA_Vec_1 = D1.*[cosd(Avg_Theta),0,sind(Avg_Theta)]; % X dimension in Lab Frame
MBA_Vec_2 = D2.*[0,1,0];  % Y dimension in Lab Frame
MBA_L_Index = 5;


%% Load default Inputs
Main_Input     = Standard_Main_Input;
TwoDGrid_Input = Standard_TwoDGrid_Input;

%% Change default Inputs with Manual Inputs
% Main Part
Main_Input.Avg_Theta  = Avg_Theta;
Main_Input.Coupling   = Coupling;
Main_Input.Sampling   = Sampling;
Main_Input.Sample_Num = Sample_Num;

% Structural part
TwoDGrid_Input.Vec_1 = MBA_Vec_1;
TwoDGrid_Input.Vec_2 = MBA_Vec_2;
TwoDGrid_Input.L_Index = MBA_L_Index;

%% Construct Molecule
G09_Path = which(G09_FileName);
Eular_MF_D = [TwoDGrid_Input.Ang_Phi,TwoDGrid_Input.Ang_Psi,TwoDGrid_Input.Ang_Theta];
Eular_MF_R = Eular_MF_D./180*pi;

G09_Input = ReadG09Input(G09_Path,Eular_MF_R);
Structure = ConstructGrid(G09_Input,TwoDGrid_Input);

%% Call Server_2DSFG
[SpectraGrid,Response] = Server_2DSFG(Structure,Main_Input);
%% Ouputs
save(SaveName,'Output')


