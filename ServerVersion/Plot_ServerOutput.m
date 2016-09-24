%% Load ouput 

Fig_Save = 1;

Default_folder = '/Users/jjho/Dropbox/Wisconsin/LAB/Data/Raw_data/Project5_FGAIL-SAM/';

[FilesName,PathName,~] = uigetfile({'*.mat','MAT files'; ...
                                    '*,*','All Files'},...
                                    'Select inputs',Default_folder,'MultiSelect','on');

NFile = length(FilesName);                                
for i = 1:NFile                                
    MAT_Path = [PathName FilesName{i}];                             
    O = load(MAT_Path);
    Output = O.Output;

    %% Figure parameters
    SG = Output.SpectraGrid;
    GUI_Inputs = Output.Main_Input;
    GUI_Inputs.LineWidth = 10;
    GUI_Inputs.Num_Contour = 100;

    %% make 2D figure
    CVL = Conv2D(SG,GUI_Inputs);
    CVL.FilesName = FilesName{i};
    hF_2D = figure;
    Plot2DSFG(hF_2D,CVL,GUI_Inputs);

    %% make Molecule figure
%     XYZ = Output.Structure.XYZ;
%     Conn = Connectivity(XYZ);
%     figure
%     gplot3(Conn,XYZ)
%     axis equal

    %% Save figures
    if Fig_Save
        SavePath = PathName;
        FigName  = regexprep(FilesName{i},'\.mat','');
        SaveName = [SavePath, FigName];
        savefig(hF_2D,SaveName)
        saveas (hF_2D,SaveName,'png')
        disp('2DSFG spectrums saved to:')
        disp(SavePath)
    end
    
    close all
end