function hF = Fig_Response(hModel, GUI_Inputs, Structure, OneDSFG, GUI_Data_hMain) 
% Plot hyper ellipsoid so that the radius = (ExJ)x(LxRxbeta)

%% define constants
N_Grid      = GUI_Inputs.Sig_NGrid;
ScaleFactor = GUI_Inputs.Sig_Scale;
Plot3D      = GUI_Inputs.Sig_Plot3D;
PlotCT      = GUI_Inputs.Sig_PlotCT;

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

%% Deal with response L x <R> x beta
% selecte mode
EigVec_Ind = GUI_Inputs.EigVec_Ind;

% Center of Exciton mode
EigVecM      = OneDSFG.H.Sort_Ex_V(2:end,2:end); % get ride of ground state
EigVecM2     = EigVecM.^2;
Center_Ex_MF = EigVecM2*(Structure.center);
Center       = Center_Ex_MF(EigVec_Ind,:);

% Response 
Response = OneDSFG.MolFrame;
% Response = OneDSFG.LabFrame;
Rho = M*Response(:,EigVec_Ind);
Rho = reshape(Rho,N_Grid,N_Grid);

% scale
Rho = Rho.*ScaleFactor;

[X,Y,Z] = sph2cart(Phi,Theta,abs(Rho));
X = X + Center(1);
Y = Y + Center(2);
Z = Z + Center(3);

Max_Rho = max(abs(Rho(:)));

%% Make 3D figure
if Plot3D
    
    % Plot molecule
    [hFunc_Model,~,~] = StructureModel(Structure.StructModel);
    hF_3D = hFunc_Model('PlotMolecule',hModel,'',guidata(hModel));
    hAx_3D = findobj(hF_3D,'type','axes');
    
    % Deal with orientation vector 
    PlotRotV = 1;

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

    % Figure setting
    % Inherent the molecular plot title
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
    hF_C = figure; 
    contour(Phi./pi.*180,Theta./pi.*180,Rho,20)
    hAx_C = findobj(hF_C,'type','axes');
    hAx_C.FontSize = 16;
    hAx_C.XGrid = 'on';
    hAx_C.YGrid = 'on';
    hAx_C.XMinorGrid = 'on';
    hAx_C.YMinorGrid = 'on';
    hAx_C.XTick = (0:60:360)';
    hAx_C.YTick = (-90:30:90)';
    xlabel(hAx_C, '\phi (Degree)');
    ylabel(hAx_C, '\theta (Degree)');
    colorbar
    colormap('jet')
    caxis([-Max_Rho,Max_Rho])
    
    % Figure setting
    % Inherent the molecular plot title
    Mode_Ind_Str  = sprintf('#%d',EigVec_Ind);
    Mode_Freq_Str = sprintf(', @%6.2f cm^{-1}' ,OneDSFG.H.Sort_Ex_Freq(EigVec_Ind+1));

    Fig_Title = ['Contour map of mode ', Mode_Ind_Str, Mode_Freq_Str];
    hAx_C.Title.String = Fig_Title;

    hF.hF_C = hF_C;
end

%% Output
if ~or(Plot3D,PlotCT)
    hF = 'no figure made';
end
