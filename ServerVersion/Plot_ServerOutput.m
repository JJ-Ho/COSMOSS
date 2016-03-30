%% Load ouput 
Default_folder = '/Users/jjho/Dropbox/Wisconsin/LAB/Data/Raw_data/Project5_FGAIL-SAM/160329_COSMOSS_MBA-FGAIL/';

[FilesName,PathName,~] = uigetfile({'*.mat','MAT files'; ...
                                    '*,*','All Files'},...
                                    'Select inputs',Default_folder);
    
MAT_Path = [PathName FilesName];                             
O = load(MAT_Path);
Output = O.Output;

%% Figure parameters
SG = Output.SpectraGrid;
GUI_Inputs = Output.Main_Input;
GUI_Inputs.LineWidth = 10;
GUI_Inputs.Num_Contour = 100;

%% make 2D figure
CVL = Conv2D(SG,GUI_Inputs);
CVL.FilesName = FilesName;
hF = figure;
Plot2DSFG(hF,CVL,GUI_Inputs);

%% make Molecule figure
XYZ = Output.Structure.XYZ;
Conn = Connectivity(XYZ);
figure
gplot3(Conn,XYZ)
axis equal