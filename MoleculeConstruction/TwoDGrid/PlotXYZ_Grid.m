function hF = PlotXYZ_Grid(hAx,SData,GUI_Inputs)
%% Input parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultPlot_Atoms      = 1;
defaultPlot_Bonds      = 1;
defaultPlot_Axis       = 1;
defaultPlot_Lattice    = 1;
defaultPlot_Atom_Index = 0;

% Add optional inputs to inputparser object
addOptional(INPUT,'Plot_Atoms'     ,defaultPlot_Atoms);
addOptional(INPUT,'Plot_Bonds'     ,defaultPlot_Bonds);
addOptional(INPUT,'Plot_Axis'      ,defaultPlot_Axis);
addOptional(INPUT,'Plot_Lattice'   ,defaultPlot_Lattice);
addOptional(INPUT,'Plot_Atom_Index',defaultPlot_Atom_Index);

parse(INPUT,GUI_Inputs_C{:});

Plot_Atoms      = INPUT.Results.Plot_Atoms;
Plot_Bonds      = INPUT.Results.Plot_Bonds;
Plot_Axis       = INPUT.Results.Plot_Axis;
Plot_Lattice    = INPUT.Results.Plot_Lattice;
Plot_Atom_Index = INPUT.Results.Plot_Atom_Index;

%% Rotate the molecule to Lab frame
XYZ      = SData.XYZ;
AtomName = SData.AtomName;
Center   = SData.LocCenter;
CoM      = SData.CoM;
Nmodes   = SData.Nmodes;
Atom_Num = SData.NAtoms;

% Grid info
N_Vec1 = SData.Extra.N_Vec1;
N_Vec2 = SData.Extra.N_Vec2;
Vec_1  = SData.Extra.Vec_1;
Vec_2  = SData.Extra.Vec_2;

%% Decide atoms types
C_Ind = strcmp(AtomName,'C');
O_Ind = strcmp(AtomName,'O');
N_Ind = strcmp(AtomName,'N');
S_Ind = strcmp(AtomName,'S');
H_Ind = strcmp(AtomName,'H');


%% Make figure
% Prep hAx 
if ~ishandle(hAx)
    hF = figure; 
    hAx = axes('Parent',hF);
else
    hF = hAx.Parent;
end

hold(hAx,'on')
    %% draw bonds
    if Plot_Bonds
        Conn = Connectivity(AtomName,XYZ);
        gplot3(Conn,XYZ,'Parent',hAx);
    end
    
    %% Add atom indexes
    if Plot_Atom_Index
        if eq(N_Vec1*N_Vec2,1)
            Atom_Ind_Str = strsplit(num2str(1:Atom_Num));
            XYZ_Atom_Ind = XYZ + 0.1;
            text(hAx,XYZ_Atom_Ind(:,1),XYZ_Atom_Ind(:,2),XYZ_Atom_Ind(:,3),Atom_Ind_Str)
        else
            disp('Atom Index only work for single molecule now...')
        end
    end
    
    %% Draw mode index
    Plot_Mode_Index = 1;
    if Plot_Mode_Index
        Mode_Ind_Str = strsplit(num2str(1:Nmodes));
        XYZ_Mode_Ind = Center + 0.1;
        text(hAx,XYZ_Mode_Ind(:,1),XYZ_Mode_Ind(:,2),XYZ_Mode_Ind(:,3),Mode_Ind_Str)
    end
    
    %% Draw atoms
    if Plot_Atoms
        plot3(hAx,Center(:,1)  ,Center(:,2)  ,Center(:,3)  ,'LineStyle','none','Marker','d','MarkerFaceColor','w')

        plot3(hAx,XYZ(C_Ind,1),XYZ(C_Ind,2),XYZ(C_Ind,3),'LineStyle','none','Marker','o','MarkerFaceColor',[0,0,0],'MarkerSize',10)
        plot3(hAx,XYZ(O_Ind,1),XYZ(O_Ind,2),XYZ(O_Ind,3),'LineStyle','none','Marker','o','MarkerFaceColor',[1,0,0],'MarkerSize',10)
        plot3(hAx,XYZ(N_Ind,1),XYZ(N_Ind,2),XYZ(N_Ind,3),'LineStyle','none','Marker','o','MarkerFaceColor',[0,0,1],'MarkerSize',10)        
        plot3(hAx,XYZ(S_Ind,1),XYZ(S_Ind,2),XYZ(S_Ind,3),'LineStyle','none','Marker','o','MarkerFaceColor',[1,1,0],'MarkerSize',10)        
        plot3(hAx,XYZ(H_Ind,1),XYZ(H_Ind,2),XYZ(H_Ind,3),'LineStyle','none','Marker','o','MarkerFaceColor',[1,1,1],'MarkerSize',5)        

    end
    
    %% Draw molecular and Lab frame
    if Plot_Axis
        Lab_Frame = [1,0,0;
                     0,1,0;
                     0,0,1 ];
                 
        PlotRotMolFrame(hAx,Lab_Frame,Lab_Frame,CoM)
    end
    
    %% Draw lattice, for MBA now
    if Plot_Lattice
        L_Vec1 = norm(Vec_1);
        L_Vec2 = norm(Vec_2);
        
        if eq(N_Vec1*N_Vec2,1)
            disp('No grid lines for single molecule...')
        else
            Conn_V1 = Connectivity('None',Center,'BondLength',L_Vec1);
            Conn_V2 = Connectivity('None',Center,'BondLength',L_Vec2);
            gplot3(or(Conn_V1,Conn_V2),Center,'Color',[0,1,0],'LineWidth',2,'Parent',hAx)
        end


    end

hold(hAx,'off')

%% Figure options
axis(hAx,'image'); %tight
rotate3d(hAx,'on')
grid(hAx,'on')
box(hAx,'on')
view(hAx,[54,12])
hAx.XLabel.String = 'X';
hAx.YLabel.String = 'Y';
hAx.ZLabel.String = 'Z';


ExtScal = 0.2;
hAx.XLim = hAx.XLim + ExtScal*sum(hAx.XLim.*[-1,1])*[-1,1];
hAx.YLim = hAx.YLim + ExtScal*sum(hAx.YLim.*[-1,1])*[-1,1];
hAx.ZLim = hAx.ZLim + ExtScal*sum(hAx.ZLim.*[-1,1])*[-1,1];
hAx.DataAspectRatio = [1,1,1];

