% batch script for 1DSFG/2DSFG generation of ideal betasheet
%% Inputs
%- Major options  ---------------------------------------------------------
Scan_Theta = 0:10:170;

% SheetType = 'Anti';
SheetType = 'Para';
N_Residue = 6;
N_Strand  = 3;
O.Coupling = 'Cho_PB';
% O.Coupling = 'Cho_APB';
Fig_Save = 1;
BaseFileName = 'PB_R6S3_';
PathName = pwd;

%- Structural part --------------------------------------------------------

TransV    = [0,0,4.75];
TwistV    = [0,0,0];
Phi_D     = 0;
Psi_D     = 0;
% Theta_D   = 0;
NLFreq    = 1644;
Anharm    = 20;

%- Labeling part ----------------------------------------------------------
O.Label_Index  = 'Not Labeled';
O.Label_Freq   = 1604;
O.LineWidth    = 5;

% For 1D ------------------------------------------------------------------
O.P_IR         = 0;
O.P_Vis1D      = 0;
O.P_Sig1D      = 0;
O.A_IR         = 90;
O.A_Vis1D      = 90;
O.A_Sig1D      = 90;
% ------------------------------------------------------------------------

% For 2D ------------------------------------------------------------------
O.A_Pump       = 90;
O.A_Probe      = 90;
O.A_Vis2D      = 90;
O.A_Sig2D      = 90;
O.P_Pump1      = 0;
O.P_Pump2      = 0;
O.P_Probe      = 0;
O.P_Vis2D      = 0;
O.P_Sig2D      = 0;
% ------------------------------------------------------------------------
  
O.Avg_Rot      = 1; %'Phi' C_Inf
O.Avg_Mirror   = 1; % no mirror plane

O.Num_Contour  = 20;
O.PlotStick    = 1;
O.PlotNorm     = 1;
O.PlotCursor   = 0;

O.F_Min     = 1550;
O.F_Max     = 1700;
O.FreqRange = O.F_Min:O.F_Max;

O.LineShape = 'L';
O.SpecType = 'Abs';
O.Pathway = 'All';
O.Signal_Type = 'Hetero';


%% Scanning 
 

for i = 1:length(Scan_Theta)

    %% Structure 
    BB        = ConstuctBetaSheet(SheetType,N_Residue,N_Strand,TransV,TwistV);
    Structure = GetAmideI(BB.Num_Atoms,...
                          BB.XYZ,...
                          BB.AtomName,...
                          BB.FilesName,...
                          'Phi_D',Phi_D,...
                          'Psi_D',Psi_D,...
                          'Theta_D',Scan_Theta(i),...
                          'NLFreq',NLFreq,...
                          'Anharm',Anharm);

    Structure.N_Residue         = N_Residue;
    Structure.N_Strand          = N_Strand;
    Structure.N_Mode_per_Starnd = N_Residue-1;

    % C terminus Index
    Structure.Ind_H = BB.Ind_H;
    Structure.Ind_O = BB.Ind_O;

    % Betasheet orientation info export
    Structure.TransV = TransV;
    Structure.TwistV = TwistV;
    Structure.RotV   = [Phi_D,Psi_D,Scan_Theta(i)];

    hF_Struc = Plot_Betasheet_AmideI(Structure);
    
    %% OneDSFG
    OneDSFG    = OneDSFG_Main(Structure,O);
    hF_OneDSFG = PlotOneDSFG(OneDSFG,O);

    %% Two DSFG
    [SpectraGrid,Response] = TwoDSFG_Main(Structure,O);

    CVL = Conv2D(SpectraGrid,O);
    CVL.FilesName = Structure.FilesName; % pass filesname for figure title
    hF_TwoDSFG = Plot2DSFG(CVL,O);
    
    %% Save figure
    if Fig_Save
        SavePath = PathName;
        FigName1  = [BaseFileName,'Th',num2str(Scan_Theta(i)),'_Structure'];
        SaveName1 = [SavePath,'/', FigName1];
        FigName2  = [BaseFileName,'Th',num2str(Scan_Theta(i)),'_1DSFG'];
        SaveName2 = [SavePath,'/', FigName2];
        FigName3  = [BaseFileName,'Th',num2str(Scan_Theta(i)),'_2DSFG'];
        SaveName3 = [SavePath,'/', FigName3];
        
        saveas(hF_Struc  ,SaveName1,'png') 
        saveas(hF_OneDSFG,SaveName2,'png')
        saveas(hF_TwoDSFG,SaveName3,'png')
        close all
    end
end