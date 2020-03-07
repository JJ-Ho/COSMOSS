function hF = PlotXYZ_Grid(hAx,SData)
%% Input parser
GUI_Inputs = SData.GUI_Inputs;
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultPlot_Atoms     = 1;
defaultPlot_Bonds     = 1;
defaultPlot_Axis      = 1;
defaultPlot_Lattice   = 0;
defaultPlot_AtomIndex = 0;
defaultMF_Center      = 1;
defaultPlot_Label     = 0;
defaultL_Index        = [];

% Add optional inputs to inputparser object
addOptional(INPUT,'Plot_Atoms'     ,defaultPlot_Atoms);
addOptional(INPUT,'Plot_Bonds'     ,defaultPlot_Bonds);
addOptional(INPUT,'Plot_Axis'      ,defaultPlot_Axis);
addOptional(INPUT,'Plot_Lattice'   ,defaultPlot_Lattice);
addOptional(INPUT,'Plot_AtomIndex' ,defaultPlot_AtomIndex);
addOptional(INPUT,'MF_Center'      ,defaultMF_Center);
addOptional(INPUT,'Plot_Label'     ,defaultPlot_Label);
addOptional(INPUT,'L_Index'        ,defaultL_Index);

parse(INPUT,GUI_Inputs_C{:});

Plot_Atoms     = INPUT.Results.Plot_Atoms;
Plot_Bonds     = INPUT.Results.Plot_Bonds;
Plot_Axis      = INPUT.Results.Plot_Axis;
Plot_Lattice   = INPUT.Results.Plot_Lattice;
Plot_AtomIndex = INPUT.Results.Plot_AtomIndex;
MF_Center      = INPUT.Results.MF_Center;
Plot_Label     = INPUT.Results.Plot_Label;
L_Index        = INPUT.Results.L_Index;

%% Rotate the molecule to Lab frame
XYZ      = SData.XYZ;
AtomName = SData.AtomName;
Center   = SData.LocCenter;
CoM      = SData.CoM;
Nmodes   = SData.Nmodes;
Atom_Num = SData.NAtoms;

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
        gplot3(Conn,XYZ,'Parent',hAx,'Color','b');     
    end
    
    %% Add atom indexes
    if Plot_AtomIndex
        N_Vec1 = GUI_Inputs.N_1;
        N_Vec2 = GUI_Inputs.N_2;
        
        if eq(N_Vec1*N_Vec2,1)
            Atom_Ind_Str = strsplit(num2str(1:Atom_Num));
            XYZ_Atom_Ind = XYZ + 0.1;
            text(hAx,XYZ_Atom_Ind(:,1),XYZ_Atom_Ind(:,2),XYZ_Atom_Ind(:,3),Atom_Ind_Str,...
                'Color','red','FontSize',14)
        else
            disp('Atom Index only work for single molecule now...')
        end
    end
    
    %% Draw mode index
    Plot_Mode_Index = 0;
    if Plot_Mode_Index
        Mode_Ind_Str = strsplit(num2str(1:Nmodes));
        XYZ_Mode_Ind = Center + 0.1;
        text(hAx,XYZ_Mode_Ind(:,1),XYZ_Mode_Ind(:,2),XYZ_Mode_Ind(:,3),Mode_Ind_Str)
    end
    
    %% Draw atoms
    if Plot_Atoms
        %plot3(hAx,Center(:,1)  ,Center(:,2)  ,Center(:,3)  ,'LineStyle','none','Marker','d','MarkerFaceColor','w')

        PlotAtom(hAx,'C',XYZ(C_Ind,:));
        PlotAtom(hAx,'O',XYZ(O_Ind,:));
        PlotAtom(hAx,'N',XYZ(N_Ind,:));
        PlotAtom(hAx,'S',XYZ(S_Ind,:));
        PlotAtom(hAx,'H',XYZ(H_Ind,:));

    end
    
    %% Draw molecular and Lab frame
    if Plot_Axis
        Lab_Frame = [1,0,0;
                     0,1,0;
                     0,0,1 ];
                 
        PlotRotMolFrame(hAx,Lab_Frame,Lab_Frame,CoM)
    end
    
    %% Draw lattice
    if Plot_Lattice
        % find the anchor points
        XYZ_N = reshape(SData.XYZ,[],SData.Nmodes,3);
        AtomCenter = reshape(XYZ_N(MF_Center,:,:),[],3);
        
        % Grid info
        N_Vec1 = GUI_Inputs.N_1;
        N_Vec2 = GUI_Inputs.N_2;
        Vec_1  = GUI_Inputs.Vec_1;
        Vec_2  = GUI_Inputs.Vec_2;
        
        L_Vec1 = norm(Vec_1);
        L_Vec2 = norm(Vec_2);
        
        if eq(N_Vec1*N_Vec2,1)
            disp('No grid lines for single molecule...')
        else
            Conn_V1_Upper = Connectivity('None',AtomCenter,'BondLength',L_Vec1);
            Conn_V1_Lower = Connectivity('None',AtomCenter,'BondLength',L_Vec1*0.9);
            Conn_V1 = and(Conn_V1_Upper,~Conn_V1_Lower);
               
            Conn_V2_Upper = Connectivity('None',AtomCenter,'BondLength',L_Vec2);
            Conn_V2_Lower = Connectivity('None',AtomCenter,'BondLength',L_Vec2*0.9);
            Conn_V2 = and(Conn_V2_Upper,~Conn_V2_Lower);         
            gplot3(or(Conn_V1,Conn_V2),AtomCenter,'Color',[0,0.5,0],'LineWidth',2,'Parent',hAx)
        end


    end

    %% Draw labeled atoms
    if Plot_Label
        if ~isempty(L_Index)
            NAtoms = SData.Children.NAtoms;
            XYZ_Monomer_Array = reshape(XYZ,NAtoms,[],3);
            Center_Ind = SData.Extra.Mol_Frame.Center_Ind;
            A1 = reshape(XYZ_Monomer_Array(Center_Ind,:,:),[],3);

            PlotAtom(hAx,'Label',A1(L_Index,:));
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

camlight
daspect([1 1 1]);

