function hF = PlotXYZ_Grid(Structure,GUI_Inputs)
%% Input parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultAvg_Phi         = 0;
defaultAvg_Theta       = 0;
defaultAvg_Psi         = 0;
defaultPlot_Atoms      = 1;
defaultPlot_Bonds      = 1;
defaultPlot_Axis       = 1;
defaultPlot_Lattice    = 1;
defaultPlot_Atom_Index = 0;

% Add optional inputs to inputparser object
addOptional(INPUT,'Avg_Phi'        ,defaultAvg_Phi);
addOptional(INPUT,'Avg_Theta'      ,defaultAvg_Theta);
addOptional(INPUT,'Avg_Psi'        ,defaultAvg_Psi);
addOptional(INPUT,'Plot_Atoms'     ,defaultPlot_Atoms);
addOptional(INPUT,'Plot_Bonds'     ,defaultPlot_Bonds);
addOptional(INPUT,'Plot_Axis'      ,defaultPlot_Axis);
addOptional(INPUT,'Plot_Lattice'   ,defaultPlot_Lattice);
addOptional(INPUT,'Plot_Atom_Index',defaultPlot_Atom_Index);

parse(INPUT,GUI_Inputs_C{:});

Avg_Phi         = INPUT.Results.Avg_Phi;
Avg_Theta       = INPUT.Results.Avg_Theta;
Avg_Psi         = INPUT.Results.Avg_Psi;
Plot_Atoms      = INPUT.Results.Plot_Atoms;
Plot_Bonds      = INPUT.Results.Plot_Bonds;
Plot_Axis       = INPUT.Results.Plot_Axis;
Plot_Lattice    = INPUT.Results.Plot_Lattice;
Plot_Atom_Index = INPUT.Results.Plot_Atom_Index;

%% Rotate the molecule to Lab frame
XYZ_MF    = Structure.XYZ;
Center_MF = Structure.center;

% Orientation = Orientation/180*pi; % turn to radius unit
Avg_Phi_R   =   Avg_Phi/180*pi;
Avg_Psi_R   =   Avg_Psi/180*pi;
Avg_Theta_R = Avg_Theta/180*pi;

R_MF_LF = R1_ZYZ_0(Avg_Phi_R,Avg_Psi_R,Avg_Theta_R);

XYZ_LF    = (R_MF_LF*XYZ_MF')';
Center_LF = (R_MF_LF*Center_MF')';

% number of replicas
N_Vec1 = Structure.N_Vec1;
N_Vec2 = Structure.N_Vec2;

%% Decide connectivity
Num_Modes = Structure.Num_Modes;
Atom_Num  = Structure.Monomer.Atom_Num;

Carbon_Pos   = XYZ_LF(Structure.AtomSerNo(:,1),:);
O_Double_Pos = XYZ_LF(Structure.AtomSerNo(:,2),:);
O_Single_Pos = XYZ_LF(Structure.AtomSerNo(:,3),:);

hF = figure;
hold on
    %% draw bonds
    if Plot_Bonds
        Conn      = Connectivity(XYZ_LF);
        % Add S-C connection for Ester
        for k = 1:Num_Modes
            Conn(13+(k-1)*Atom_Num,3+(k-1)*Atom_Num) = 1; 
        end
        gplot3(Conn,XYZ_LF);
    end
    
    %% Add atom indexes
    if Plot_Atom_Index
        if eq(N_Vec1*N_Vec2,1)
            Index_Str = strsplit(num2str(1:Atom_Num));
            XYZ_Text = XYZ_LF + 0.1;
            text(XYZ_Text(:,1),XYZ_Text(:,2),XYZ_Text(:,3),Index_Str)
        else
            disp('Atom Index only work for single molecule now...')
        end
    end
    
    %% Draw atoms
    if Plot_Atoms
        plot3(Center_LF(:,1)    ,Center_LF(:,2)    ,Center_LF(:,3)   ,'LineStyle','none','Marker','d','MarkerFaceColor','w')
        plot3(Carbon_Pos(:,1)   ,Carbon_Pos(:,2)   ,Carbon_Pos(:,3)  ,'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerSize',10)
        plot3(O_Double_Pos(:,1) ,O_Double_Pos(:,2) ,O_Double_Pos(:,3),'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',10)
        plot3(O_Single_Pos(:,1) ,O_Single_Pos(:,2) ,O_Single_Pos(:,3),'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',10)
    end
    
    %% Draw molecular and Lab frame
    if Plot_Axis
        hAx = findobj(hF,'type','axes');
        Lab_Frame = [1,0,0;
                     0,1,0;
                     0,0,1 ];
                 
        Box_Cornor = [hAx.XLim(1),hAx.YLim(1),hAx.ZLim(1)];
        Box_length = [hAx.XLim(2),hAx.YLim(2),hAx.ZLim(2)] - [hAx.XLim(1),hAx.YLim(1),hAx.ZLim(1)];
        COA = Box_Cornor + (0.1).*Box_length;
        
        PlotRotMolFrame(hAx,Lab_Frame,R_MF_LF,COA)
    end
    
    %% Draw lattice, for MBA now
    if Plot_Lattice
        Vec_1  = Structure.Vec_1;
        Vec_2  = Structure.Vec_2;
        L_Vec1 = norm(Vec_1);
        L_Vec2 = norm(Vec_2);
        
        if eq(N_Vec1*N_Vec2,1)
            disp('No grid lines for single molecule...')
        else
            XYZ_Tmp = reshape(XYZ_LF,[],N_Vec1*N_Vec2,3);
            XYZ_S   = squeeze(XYZ_Tmp(13,:,:));

            Conn_V1_Max = Connectivity(XYZ_S,'BondLength',L_Vec1+0.1);
            Conn_V1_Min = Connectivity(XYZ_S,'BondLength',L_Vec1-0.1);
            Conn_V1 = and(Conn_V1_Max,~Conn_V1_Min);

            Conn_V2_Max = Connectivity(XYZ_S,'BondLength',L_Vec2+0.1);
            Conn_V2_Min = Connectivity(XYZ_S,'BondLength',L_Vec2-0.1);
            Conn_V2 = and(Conn_V2_Max,~Conn_V2_Min);

            gplot3(or(Conn_V1,Conn_V2),XYZ_S,'Color',[0,1,0],'LineWidth',2)
        end


    end
    hold off

%% Figure options
hAx = findobj(hF,'type','axes');
axis tight;
rotate3d on
view([54,12])
box on ;

ExtScal = 0.2;
hAx.XLim = hAx.XLim + ExtScal*sum(hAx.XLim.*[-1,1])*[-1,1];
hAx.YLim = hAx.YLim + ExtScal*sum(hAx.YLim.*[-1,1])*[-1,1];
hAx.ZLim = hAx.ZLim + ExtScal*sum(hAx.ZLim.*[-1,1])*[-1,1];
hAx.DataAspectRatio = [1,1,1];

xlabel('X')
ylabel('Y')
