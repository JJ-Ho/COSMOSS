function Output = Fig_Response_Scalar(hModel, GUI_Inputs, Structure, SpecData, GUI_Data_hMain) 
% Orientation dependence of SFG/2DSFG signal in scalar form

%% Debug
% N_Grid      = 30;
% ScaleFactor = 1;
% Plot3D      = 0;
% PlotCT      = 1;
% PlotSum     = 0;
% SpecType    = 2;
% 
% GUI_Data_hMain.A_IR    = 90;
% GUI_Data_hMain.A_Vis1D = 90;
% GUI_Data_hMain.A_Sig1D = 90;
% 
% GUI_Data_hMain.P_IR    = 0;
% GUI_Data_hMain.P_Vis1D = 0;
% GUI_Data_hMain.P_Sig1D = 0;
% 
% SpecData = OneDSFG_Main(Structure,GUI_Data_hMain);
% % SpecData = TwoDSFG_Main(Structure,GUI_Data_hMain);
% 
% % EigVec_Ind = [2,3,7,8];
% EigVec_Ind = 3;

%% define constants
N_Grid      = GUI_Inputs.Sig_NGrid;
ScaleFactor = GUI_Inputs.Sig_Scale;
Plot3D      = GUI_Inputs.Sig_Plot3D;
PlotCT      = GUI_Inputs.Sig_PlotCT;
PlotSum     = GUI_Inputs.Sig_PlotSum;

SpecType    = GUI_Inputs.SpecType;
EigVec_Ind  = GUI_Inputs.Mu_Alpha_Ind; % get a full mode list from Mu_Alpha_Ind

%% Generate ExJ(psi,theta,Ai...,Pi...) and Signal
% The relative orientation of laser beams are determined by the incidnet 
% angles(Ai) and the polarization angle (Pi) defined on the incident beam 
% path. The molecular response to this laser setup configuration is a 
% function of the molecular orientation relative to the fixed experiment
% frame. Instead of viewing the response as function of molecule rotatein
% the exp. frame, we can fix the molecule and rotate the whole laser setup
% instead. These two picture are mathematically equivenlent but cause less
% computation if we apply the rotation on ExJ on each laer beam follow by 
% tensor product each of ExJ. This is effectly a [1x2]*[2x3]*[3x3]
% operation and follow by tensor product of however many [1x3] vectors. It
% is way faster than generating [3^Nx3^N] rotational matrix to rotate the
% higher order (N =3 for SFG, =5 for 2DSFG) molecular response. 
% If we further assume that the molecule system is azimutal symetric then 
% we can eliminate the phi angle out of three Euler angle and thus we can
% map the response (R) as function of two Euler angle so we can visulize it
% in (R,psi,theta) 3D coordinate. 

[EJR2,EJR1,Phi,Theta] = EJR_Scalar(GUI_Data_hMain,N_Grid);

switch SpecType
    case 2 % 1DSFG
        A_All = SpecData.Alpha.M_Ex_01;
        M_All    = SpecData.Mu.M_Ex_01;
        
        Alpha = A_All(EigVec_Ind,:)';
        Mu    =    M_All(EigVec_Ind,:)';
        Sig   = (EJR2*Alpha).*(EJR1*Mu);
        
    case 5 % 2DSFG
        % decode the unified index to type + index
        PathType_Total = fieldnames(SpecData.Index);
        TypeVec  = [];
        IndexVec = [];
        for P = 1:length(PathType_Total)
            N_Path   = size(SpecData.Index.(PathType_Total{P}),1);
            TypeVec  = [TypeVec;ones(N_Path,1).*P];
            IndexVec = [IndexVec;(1:N_Path)'];
        end
        PathType = TypeVec(EigVec_Ind);
        TypeInd = IndexVec(EigVec_Ind);
        
        F = SpecData.Index.(PathType_Total{PathType})(TypeInd,:);
        
        
        A1_All = SpecData.Alpha.M_Ex_01;
        A2_All = SpecData.Alpha.M_Ex_12;
        M1_All = SpecData.Mu.M_Ex_01;
        M2_All = SpecData.Mu.M_Ex_12;
        
        switch PathType
            case 1
                a = F(4);
                b = F(2);

                A4 = A1_All(b,:)';
                M3 = M1_All(b,:)';
                M2 = M1_All(a,:)';
                M1 = M1_All(a,:)';
            case 2
                a = F(4);
                b = F(3);

                A4 = A1_All(b,:)';
                M3 = M1_All(a,:)';
                M2 = M1_All(b,:)';
                M1 = M1_All(a,:)';
            case 3
                a = F(4);
                b = F(3);
                x = F(2);
                
                A4 = -squeeze(A2_All(a,x,:));
                M3 =  squeeze(M2_All(b,x,:));
                M2 = M1_All(b,:)';
                M1 = M1_All(a,:)';
            case 4
                a = F(4);
                b = F(2);

                A4 = A1_All(b,:)';
                M3 = M1_All(b,:)';
                M2 = M1_All(a,:)';
                M1 = M1_All(a,:)';
            case 5
                a = F(4);
                b = F(3);

                A4 = A1_All(a,:)';
                M3 = M1_All(b,:)';
                M2 = M1_All(b,:)';
                M1 = M1_All(a,:)';
            case 6
                a = F(4);
                b = F(3);
                x = F(2);
                
                A4 = -squeeze(A2_All(b,x,:));
                M3 =  squeeze(M2_All(a,x,:));
                M2 = M1_All(b,:)';
                M1 = M1_All(a,:)';
        end

        Sig   = (EJR2*A4).*(EJR1*M3).*(EJR1*M2).*(EJR1*M1);
                
    otherwise
        disp('Spectral type has not supported yet...')
        Output = 'Null';
        return
end

Theta_D = Theta./pi.*180;
Phi_D   =   Phi./pi.*180;

%% Deal with Grid and scale
% Response Grid
Rho = reshape(Sig,N_Grid,2*N_Grid,[]);

% scale
N_Modes = length(EigVec_Ind);
Rho = Rho.*ScaleFactor;
Max_Rho = max(abs(reshape(Rho,[],N_Modes)));

%% Sum up multimodes response if needed
if PlotSum
    Rho = sum(Rho,3);
    Max_Rho = sum(Max_Rho);
    N_Modes = 1;
end

%% Make 3D hyper-ellipsoid
hF = struct;
if Plot3D
    
    % Plot molecule
    [hFunc_Model,~,~] = StructureModel(Structure.StructModel);
    hF_3D = hFunc_Model('PlotMolecule',hModel,'',guidata(hModel));
    hAx_3D = findobj(hF_3D,'type','axes');
    
    % Deal with orientation vector 
    PlotRotV = 1;

    % Center of Exciton mode
    EigVecM      = SpecData.H.Sort_Ex_V(2:end,2:end); % get ride of ground state
    EigVecM2     = EigVecM.^2;
    Center_Ex_MF = EigVecM2*(Structure.LocCenter);
    
    if PlotSum
        Center = sum(Center_Ex_MF(EigVec_Ind,:),1)./length(EigVec_Ind);
    else
        Center = Center_Ex_MF(EigVec_Ind,:);
        %Center = [0,0,0]; % tmp work-around for tensor form 2DSFG 
    end

    % shift hyperellipsoid to mode center
    Az = Phi;
    Elv = pi/2-Theta;
    [X,Y,Z] = sph2cart(Az,Elv,abs(Rho));
    X = X + Center(1);
    Y = Y + Center(2);
    Z = Z + Center(3);

    hold on
        if PlotRotV

            RotV = [0;0;1];
            RotV = RotV.*Max_Rho*1.1;

            quiver3(hAx_3D,...
                    Center(1),Center(2),Center(3),...
                    RotV(1),RotV(2),RotV(3),0,...
                    'LineWidth',2,...
                    'MaxHeadSize',0.6,...
                    'Color',[0,0,1]);
        end

        hSurf = surf(hAx_3D,X,Y,Z,sign(Rho)); % colormapping sign only
    hold off

    % colorbar
    colormap('cool')
    caxis([-1,1])
    
    % Adjust the surface
    Transparency = 0.5;
    hSurf.EdgeColor = 'interp';
    hSurf.FaceAlpha = Transparency;
    hSurf.EdgeAlpha = Transparency;

    % Figure title inherent the molecular plot
    Fig_Title = hAx_3D.Title.String;
    Mode_Ind_Str  = sprintf('#%d',EigVec_Ind);
    Mode_Freq_Str = sprintf(', @%6.2f cm^{-1}' ,SpecData.H.Sort_Ex_Freq(EigVec_Ind+1));
    Scaling_Str   = sprintf(', Scale= %2.1f',ScaleFactor);

    Fig_Title{length(Fig_Title)+1} = [Mode_Ind_Str, Mode_Freq_Str, Scaling_Str];
    hAx_3D.Title.String = Fig_Title;
    
    hF.hF_3D = hF_3D;
end

%% Contour plot
if PlotCT
    for j = 1:N_Modes
        hF_C = figure; 
        contourf(Phi_D,Theta_D,Rho(:,:,j),20)

        % figure adjustment
        hAx_C = findobj(hF_C,'type','axes');
        hAx_C.FontSize = 16;
        hAx_C.XGrid = 'on';
        hAx_C.YGrid = 'on';
        hAx_C.XMinorGrid = 'on';
        hAx_C.YMinorGrid = 'on';
        hAx_C.XTick = (0:60:360)';
        hAx_C.YTick = (0:30:180)';
        xlabel(hAx_C, '\phi (Degree)');
        ylabel(hAx_C, '\theta (Degree)');
        colorbar
        colormap('jet')
        caxis([-Max_Rho(j),Max_Rho(j)])
        
        % Figure title inherent the molecular plot
        if PlotSum
            M_Ind = EigVec_Ind;
            Mode_Ind_Str  = sprintf('#%d ',M_Ind);
            Fig_Title = ['Contour map of mode ', Mode_Ind_Str];

        else
            switch SpecType
                case 2
                    M_Ind = EigVec_Ind(j);
                    Mode_Ind_Str  = sprintf('#%d',M_Ind);
                    Mode_Freq_Str = sprintf(', @%6.2f cm^{-1}' ,SpecData.H.Sort_Ex_Freq(M_Ind+1));
                    Fig_Title = ['Contour map of mode ', Mode_Ind_Str, Mode_Freq_Str];
            
                case 4
                    PathType_Str = PathType_Total{PathType};
                    PathInd_Str = sprintf('%s,%d,%d,%d,%d',PathType_Str,F(1),F(2),F(3),F(4));
                    Fig_Title = ['Contour map of mode ', PathInd_Str];
            end
        end
        
        hAx_C.Title.String = Fig_Title;

        hF.hF_C = hF_C;
    end
end

%% Output
Output.hF    = hF;
Output.Phi   = Phi;
Output.Theta = Theta;
Output.Rho   = Rho;
