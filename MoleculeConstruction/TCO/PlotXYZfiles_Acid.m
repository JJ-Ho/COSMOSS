function hF = PlotXYZfiles_Acid(Struc_Data,GUI_Inputs)
%% PlotXYZfiles
%  
% Given the Stuctural infomation that generated by GetAmideI.m this
% function can plot the molecule and its' amide I modes.
% 
% ------- Version log -----------------------------------------------------
% 
% Ver. 1.5  140723  Add hydrogen atom
% 
% Ver. 1.4  140605  Readuce the process for reading XYZ data
% 
% Ver. 1.3  131107  Add imput "Rot_Mat" to rotate molecue as user assigned
%                   This is the correction follow by v1.2 of
%                   OneDSFG_AmideI.m
% 
% Ver. 1.2  130930  Move gplot3 from Molecular_Plot to here to accomodate
%                   the update of molecular_plot
% 
% Ver. 1.1  130814  Bug in atom position index fixed 
% 
% Ver. 1.0  130729  Isolated from TwoDSFG_Simulation
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2013

%% Debug
% PDB_Data    = GetAcid('Debug');
% Orientation = [0,0,0];

%% Input parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultAvg_Phi    = 0;
defaultAvg_Theta  = 0;
defaultAvg_Psi    = 0;
defaultPlot_Atoms = 1;
defaultPlot_Bonds = 1;
defaultPlot_Axis  = 1;

% Add optional inputs to inputparser object
addOptional(INPUT,'Avg_Phi'   ,defaultAvg_Phi);
addOptional(INPUT,'Avg_Theta' ,defaultAvg_Theta);
addOptional(INPUT,'Avg_Psi'   ,defaultAvg_Psi);
addOptional(INPUT,'Plot_Atoms',defaultPlot_Atoms);
addOptional(INPUT,'Plot_Bonds',defaultPlot_Bonds);
addOptional(INPUT,'Plot_Axis' ,defaultPlot_Axis);

parse(INPUT,GUI_Inputs_C{:});

Avg_Phi    = INPUT.Results.Avg_Phi;
Avg_Theta  = INPUT.Results.Avg_Theta;
Avg_Psi    = INPUT.Results.Avg_Psi;
Plot_Atoms = INPUT.Results.Plot_Atoms;
Plot_Bonds = INPUT.Results.Plot_Bonds;
Plot_Axis  = INPUT.Results.Plot_Axis;

%% Rotate the molecule to Lab frame

XYZ_MF          = Struc_Data.XYZ;
Center_MF       = Struc_Data.center;
Displacement_MF = Struc_Data.Displacement;

% Orientation = Orientation/180*pi; % turn to radius unit
Avg_Phi_R   =   Avg_Phi/180*pi;
Avg_Psi_R   =   Avg_Psi/180*pi;
Avg_Theta_R = Avg_Theta/180*pi;

R_MF_LF = R1_ZYZ_0(Avg_Phi_R,Avg_Psi_R,Avg_Theta_R);

XYZ_LF = (R_MF_LF*XYZ_MF')';
Center_LF = (R_MF_LF*Center_MF')';
Displacement_LF = (R_MF_LF*Displacement_MF')';

%% Define C,OD,OS atom position
Carbon_Pos   = XYZ_LF(Struc_Data.AtomSerNo(:,1),:);
OxygenD_Pos  = XYZ_LF(Struc_Data.AtomSerNo(:,2),:);
OxygenS_Pos  = XYZ_LF(Struc_Data.AtomSerNo(:,3),:);
Hydrogen_Pos = XYZ_LF(Struc_Data.AtomSerNo(:,4),:);

hF = figure; 
hold on

    %% draw bonds
    if Plot_Bonds
        Conn = Connectivity(XYZ_LF);
            % add line between the two Carbon
            Conn(1,5) = 1;
            Conn(5,1) = 1;

        gplot3(Conn,XYZ_LF);
    end
    %%  Draw atoms
    if Plot_Atoms
        plot3(Center_LF(:,1)  ,Center_LF(:,2)  ,Center_LF(:,3)  ,'LineStyle','none','Marker','d','MarkerFaceColor','w','MarkerSize',10)
        plot3(Carbon_Pos(:,1)  ,Carbon_Pos(:,2)  ,Carbon_Pos(:,3)  ,'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerSize',10)
        plot3(OxygenD_Pos(:,1) ,OxygenD_Pos(:,2) ,OxygenD_Pos(:,3) ,'LineStyle','none','Marker','o','MarkerFaceColor',[1,0,0],'MarkerSize',10)
        plot3(OxygenS_Pos(:,1) ,OxygenS_Pos(:,2) ,OxygenS_Pos(:,3) ,'LineStyle','none','Marker','o','MarkerFaceColor',[1,0,0],'MarkerSize',10)
        plot3(Hydrogen_Pos(:,1),Hydrogen_Pos(:,2),Hydrogen_Pos(:,3),'LineStyle','none','Marker','o','MarkerFaceColor',[1,1,1],'MarkerSize',10)
    end
    
    %% Draw molecular and Lab frame
    if Plot_Axis
        hAx = findobj(hF,'type','axes');
        Lab_Frame = [1,0,0;
                     0,1,0;
                     0,0,1 ];
        COA = [-3,-3,-3];
        PlotRotMolFrame(hAx,Lab_Frame,R_MF_LF,COA)
    end
hold off

%% figure setting 
hAx = findobj(hF,'type','axes');

LL = 5;
Orig_Box = [-LL,LL,-LL,LL,-LL,LL];

XMin = min( Orig_Box(1),Orig_Box(1)+ Displacement_LF(1));
XMax = max( Orig_Box(2),Orig_Box(2)+ Displacement_LF(1));
YMin = min( Orig_Box(3),Orig_Box(3)+ Displacement_LF(2));
YMax = max( Orig_Box(4),Orig_Box(4)+ Displacement_LF(2));
ZMin = min( Orig_Box(5),Orig_Box(5)+ Displacement_LF(3));
ZMax = max( Orig_Box(6),Orig_Box(6)+ Displacement_LF(3));

hAx.XLim = [XMin,XMax];
hAx.YLim = [YMin,YMax];
hAx.ZLim = [ZMin,ZMax];

rotate3d on
grid on


xlabel('X')
ylabel('Y')
zlabel('Z')

view([0,0])
