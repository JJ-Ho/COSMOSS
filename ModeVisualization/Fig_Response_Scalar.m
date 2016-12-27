% function Output = Fig_Response_Scalar(hModel, GUI_Inputs, Structure, SpecData, GUI_Data_hMain) 
% Plot hyper ellipsoid so that the radius = (ExJ)x(LxRxbeta)
%% Debug
N_Grid      = 30;
ScaleFactor = 1;
Plot3D      = 0;
PlotCT      = 1;
PlotSum     = 0;
SpecType    = 2;

GUI_Data_hMain.A_IR    = 90;
GUI_Data_hMain.A_Vis1D = 90;
GUI_Data_hMain.A_Sig1D = 90;

GUI_Data_hMain.P_IR    = 0;
GUI_Data_hMain.P_Vis1D = 0;
GUI_Data_hMain.P_Sig1D = 0;

SpecData = OneDSFG_Main(Structure,GUI_Data_hMain);
% SpecData = TwoDSFG_Main(Structure,GUI_Data_hMain);

% EigVec_Ind = [2,3,7,8];
EigVec_Ind = 3;

%% define constants
% N_Grid      = GUI_Inputs.Sig_NGrid;
% ScaleFactor = GUI_Inputs.Sig_Scale;
% Plot3D      = GUI_Inputs.Sig_Plot3D;
% PlotCT      = GUI_Inputs.Sig_PlotCT;
% PlotSum     = GUI_Inputs.Sig_PlotSum;

% SpecType    = GUI_Inputs.SpecType;
% EigVec_Ind  = GUI_Inputs.Mu_Alpha_Ind; % get a full mode list from Mu_Alpha_Ind

%% Generate ExJ(psi,theta,Ai...,Pi...)
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


switch SpecType
    case 2 % 1DSFG
        [M2,M1,Phi,Theta] = EJR_Scalar(GUI_Data_hMain,N_Grid);
        
        Alpha_All = SpecData.Alpha.M_Ex_01;
        Mu_All    = SpecData.Mu.M_Ex_01;
        
        Alpha = Alpha_All(EigVec_Ind,:)';
        Mu    =    Mu_All(EigVec_Ind,:)';
        Rho   = (M2*Alpha).*(M1*Mu);
    case 4 % 2DSFG
        disp('Spectral Type no supported yet...')
        %[M,Phi,Theta] = EJRR_2DSFG(GUI_Data_hMain,N_Grid);
        %Response = SpecData.MolFrame.R1; %% Need to implement which pathway to use
    otherwise
        disp('Spectral Type no supported yet...')
        return
end

Theta_D = Theta./pi.*180;
Phi_D   =   Phi./pi.*180;

%% Deal with response L x <R> x beta
% Response 
Rho = reshape(Rho,N_Grid,2*N_Grid,[]);

% scale
N_Modes = length(EigVec_Ind);
Rho = Rho.*ScaleFactor;
Max_Rho = max(abs(reshape(Rho,[],N_Modes)));

%% Sum up multimodes response
if PlotSum
    Rho = sum(Rho,3);
    Max_Rho = sum(Max_Rho);
    N_Modes = 1;
end

%% Make 3D figure
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
    Center_Ex_MF = EigVecM2*(Structure.center);
    
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
            M_Ind = EigVec_Ind(j);
            Mode_Ind_Str  = sprintf('#%d',M_Ind);
            Mode_Freq_Str = sprintf(', @%6.2f cm^{-1}' ,SpecData.H.Sort_Ex_Freq(M_Ind+1));
            Fig_Title = ['Contour map of mode ', Mode_Ind_Str, Mode_Freq_Str];
        end
        
        hAx_C.Title.String = Fig_Title;

        hF.hF_C = hF_C;
    end
end

%% Make contour ratio plot
% if PlotCT_R
% %     reverse = 1./(1:0.1:5);
% %     v = [5:-0.1:2,reverse];
%     for k = 2:N_Modes
%         hF_CR = figure; 
%         contourf(Phi_D,Theta_D,Rho_R(:,:,k),100)
% 
%         % figure adjustment
%         hAx_CR = findobj(hF_CR,'type','axes');
%         hAx_CR.FontSize = 16;
%         hAx_CR.XGrid = 'on';
%         hAx_CR.YGrid = 'on';
%         hAx_CR.XMinorGrid = 'on';
%         hAx_CR.YMinorGrid = 'on';
%         hAx_CR.XTick = (0:60:360)';
%         hAx_CR.YTick = (-90:30:90)';
%         xlabel(hAx_CR, '\phi (Degree)');
%         ylabel(hAx_CR, '\theta (Degree)');
%         
%         % Set colorbar
%         colorbar
%         CMAP = SelectColormap(7);
%         colormap(CMAP)   
%         
%         Max_Rho_R = max(abs(reshape(Rho_R(:,:,k),[],1)));
%         caxis([-Max_Rho_R,Max_Rho_R])
%         %caxis([0.2,5])
% 
%         
%         % Figure title inherent the molecular plot
%         Mode_Ind_Str  = sprintf('#%d to %d',EigVec_Ind(k),EigVec_Ind(1));
%         Mode_Freq_Str = sprintf(', @%6.2f cm^{-1}' ,OneDSFG.H.Sort_Ex_Freq(EigVec_Ind(k)+1));
% 
%         Fig_Title = ['Ratio map of mode ', Mode_Ind_Str, Mode_Freq_Str];
%         hAx_CR.Title.String = Fig_Title;
% 
%         hF.hF_CR = hF_CR;
%     end
% end

%% Output
if ~or(Plot3D,PlotCT)
    hF = 'no figure made';
end

Output.hF    = hF;
Output.Phi   = Phi;
Output.Theta = Theta;
Output.Rho   = Rho;
