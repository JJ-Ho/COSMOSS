function Output = Fig_Response(hModel, GUI_Inputs, Structure, OneDSFG, GUI_Data_hMain) 
% Plot hyper ellipsoid so that the radius = (ExJ)x(LxRxbeta)

%% define constants
N_Grid      = GUI_Inputs.Sig_NGrid;
ScaleFactor = GUI_Inputs.Sig_Scale;
Plot3D      = GUI_Inputs.Sig_Plot3D;
PlotCT      = GUI_Inputs.Sig_PlotCT;
PlotSum     = GUI_Inputs.Sig_PlotSum;

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

% [M,Phi,Theta] = EJ(Polar,Ang,N_Grid);
[M,Phi,Theta] = EJRR(GUI_Data_hMain,N_Grid);

Theta_D = Theta./pi.*180;
Phi_D   =   Phi./pi.*180;

%% Deal with response L x <R> x beta
% selecte mode
EigVec_Ind = GUI_Inputs.Mu_Alpha_Ind; % get a full mode list from Mu_Alpha_Ind

% Response 
Response = OneDSFG.MolFrame;
Rho = M*Response(:,EigVec_Ind);
Rho = reshape(Rho,N_Grid,2*N_Grid,[]);

% scale
N_Modes = length(EigVec_Ind);
Rho = Rho.*ScaleFactor;
Max_Rho = max(abs(reshape(Rho,[],N_Modes)));

% % deal with multi modes
% if gt(N_Modes,1)
%     disp('Multiple mode comparison, ')
%     disp('Contour plot shows the ratio of all the modes raletive to the first mode')
%     Plot3D = 0;
%     
%     Rho_R = zeros(size(Rho) - [0,0,1]);
% 
%     for i = 2:N_Modes
%         %Rho_R(:,:,i) = Rho(:,:,i)./Rho(:,:,1);
%         %Rho_R(:,:,i) = (Rho(:,:,i)./Max_Rho(i)) ./ (Rho(:,:,1)./Max_Rho(1));
%         %Rho_R(:,:,i) = (Rho(:,:,i)./Max_Rho(i)) - (Rho(:,:,1)./Max_Rho(1));
%         %Rho_R(:,:,i) = 1./(Rho(:,:,i)./Max_Rho(i)) - (Rho(:,:,1)./Max_Rho(1));
%         Rho_R(:,:,i) = (Rho(:,:,i)./Max_Rho(i)) - (Rho(:,:,1)./Max_Rho(1));
%         %Rho_R(:,:,i) = (Rho(:,:,i)+Max_Rho(i)) ./ (Rho(:,:,1)+Max_Rho(1));
%         %Rho_R(:,:,i) = (Rho(:,:,i)+1) ./ (Rho(:,:,1)+1);
%         %Rho_R(:,:,i) = (Rho(:,:,i)+Max_Rho(i)) ./ (Rho(:,:,1)+Max_Rho(1));
%         %Rho_R(:,:,i) = (abs(Rho(:,:,i))+1) ./ (abs(Rho(:,:,1))+1);
%         %Rho_R(:,:,i) = (abs(Rho(:,:,1))+1) ./ (abs(Rho(:,:,i))+1);
%     end
% end

%% Sum up multimodes response
if PlotSum
    Rho = sum(Rho,3);
    Max_Rho = sum(Max_Rho);
    N_Modes = 1;
end

%% Make 3D figure
if Plot3D
    
    % Plot molecule
    [hFunc_Model,~,~] = StructureModel(Structure.StructModel);
    hF_3D = hFunc_Model('PlotMolecule',hModel,'',guidata(hModel));
    hAx_3D = findobj(hF_3D,'type','axes');
    
    % Deal with orientation vector 
    PlotRotV = 1;

    % Center of Exciton mode
    EigVecM      = OneDSFG.H.Sort_Ex_V(2:end,2:end); % get ride of ground state
    EigVecM2     = EigVecM.^2;
    Center_Ex_MF = EigVecM2*(Structure.center);
    Center       = Center_Ex_MF(EigVec_Ind,:);

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
                    'Color',[255,128,0]./256);
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
    Mode_Freq_Str = sprintf(', @%6.2f cm^{-1}' ,OneDSFG.H.Sort_Ex_Freq(EigVec_Ind+1));
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
        Mode_Ind_Str  = sprintf('#%d',EigVec_Ind(j));
        Mode_Freq_Str = sprintf(', @%6.2f cm^{-1}' ,OneDSFG.H.Sort_Ex_Freq(EigVec_Ind(j)+1));

        Fig_Title = ['Contour map of mode ', Mode_Ind_Str, Mode_Freq_Str];
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
