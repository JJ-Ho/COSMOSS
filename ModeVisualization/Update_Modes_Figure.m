function Output = Update_Modes_Figure(hF,GUI_Inputs, Structure, OneDSFG)
%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultAvg_Phi    = 0;
defaultAvg_Theta  = 0;
defaultAvg_Psi    = 0;
defaultPlot_Loc       = 0;
defaultLoc_Ind        = [];
defaultPlot_Ex        = 0;
defaultEx_Ind         = [];
defaultPlot_TDV       = 1;
defaultScale_TDV      = 1;
defaultPlot_Raman     = 1;
defaultScale_Raman    = 1;
defaultNormalize      = 1;
defaultPlot_EigenVec  = 0;
defaultEigneVec_Ind   = [];


% add options
addOptional(INPUT,'Avg_Phi'   ,defaultAvg_Phi);
addOptional(INPUT,'Avg_Theta' ,defaultAvg_Theta);
addOptional(INPUT,'Avg_Psi'   ,defaultAvg_Psi);
addOptional(INPUT,'Plot_Loc'      , defaultPlot_Loc      );
addOptional(INPUT,'Loc_Ind'       , defaultLoc_Ind       );
addOptional(INPUT,'Plot_Ex'       , defaultPlot_Ex       );
addOptional(INPUT,'Ex_Ind'        , defaultEx_Ind        );
addOptional(INPUT,'Plot_TDV'      , defaultPlot_TDV      );
addOptional(INPUT,'Scale_TDV'     , defaultScale_TDV     );
addOptional(INPUT,'Plot_Raman'    , defaultPlot_Raman    );
addOptional(INPUT,'Scale_Raman'   , defaultScale_Raman   );
addOptional(INPUT,'Normalize'     , defaultNormalize     );
addOptional(INPUT,'Plot_EigenVec' , defaultPlot_EigenVec );
addOptional(INPUT,'EigneVec_Ind'  , defaultEigneVec_Ind  );

parse(INPUT,GUI_Inputs_C{:});

% Reassign Variable names
Avg_Phi    = INPUT.Results.Avg_Phi;
Avg_Theta  = INPUT.Results.Avg_Theta;
Avg_Psi    = INPUT.Results.Avg_Psi;
Plot_Loc       = INPUT.Results.Plot_Loc      ;
Loc_Ind        = INPUT.Results.Loc_Ind       ;
Plot_Ex        = INPUT.Results.Plot_Ex       ;
Ex_Ind         = INPUT.Results.Ex_Ind        ;
Plot_TDV       = INPUT.Results.Plot_TDV      ;
Scale_TDV      = INPUT.Results.Scale_TDV     ;
Plot_Raman     = INPUT.Results.Plot_Raman    ;
Scale_Raman    = INPUT.Results.Scale_Raman   ;
Normalize      = INPUT.Results.Normalize     ;
Plot_EigenVec  = INPUT.Results.Plot_EigenVec ;
EigneVec_Ind   = INPUT.Results.EigneVec_Ind  ;

%% Generate distinguisable colors for modes
if Plot_EigenVec
    Loc_Ind = 1:Structure.Num_Modes;
end

N_Loc_Mode = length(Loc_Ind);
N_Ex_Mode  = length(Ex_Ind );
Total_N_Plot_Modes = N_Loc_Mode + N_Ex_Mode;

if exist('distinguishable_colors','file')
    UnWanted = [0,0,0;1,1,1;1,0,0;0,1,0;0,0,1];
    Mode_colors = distinguishable_colors(Total_N_Plot_Modes,UnWanted);
else
    Mode_colors = bsxfun(@times,ones(Total_N_Plot_Modes,3),[255,128,0]./256);
end
Loc_Mode_colors = Mode_colors(           1:        N_Loc_Mode,:);
Ex_Mode_colors  = Mode_colors(N_Loc_Mode+1:Total_N_Plot_Modes,:);

%% Rotate from molecule frame to lab frame

% Orientation = Orientation/180*pi; % turn to radius unit
Avg_Phi_R   =   Avg_Phi/180*pi;
Avg_Psi_R   =   Avg_Psi/180*pi;
Avg_Theta_R = Avg_Theta/180*pi;
R_MF_LF     = R1_ZYZ_0(Avg_Phi_R,Avg_Psi_R,Avg_Theta_R);

%% molecular frame
% center
Center_Loc_MF = Structure.center(Loc_Ind,:);

EigVecM      = OneDSFG.H.Sort_Ex_V(2:end,2:end); % get ride of ground state
EigVecM2     = EigVecM.^2;
Center_Ex_MF = EigVecM2*(Structure.center);
Center_Ex_MF = Center_Ex_MF(Ex_Ind,:);

% Transition dipole
Mu_Loc_MF     = squeeze(OneDSFG.Mu.Trans_Loc(1,Loc_Ind+1,:)); % shift by 1 to avoid ground state
Mu_Ex_MF      = squeeze(OneDSFG.Mu. Trans_Ex(1, Ex_Ind+1,:)); % shift by 1 to avoid ground state

% Raman Tensor
Alpha_Loc_MF  = squeeze(OneDSFG.Alpha.Trans_Loc(1,Loc_Ind+1,:));
Alpha_Ex_MF   = squeeze(OneDSFG.Alpha. Trans_Ex(1, Ex_Ind+1,:));

%% lab frame
% permute the matix dimension for spectial case
if eq(N_Loc_Mode,1)
    Mu_Loc_MF = Mu_Loc_MF';
end

if eq(N_Ex_Mode,1)
    Mu_Ex_MF = Mu_Ex_MF';
end

% center
Center_Loc_LF = (R_MF_LF*Center_Loc_MF')';
Center_Ex_LF  = (R_MF_LF*Center_Ex_MF')';

% transition dipole
Mu_Loc_LF = (R_MF_LF*Mu_Loc_MF')';
Mu_Ex_LF  = (R_MF_LF*Mu_Ex_MF')';

% Raman tensor
AlphaM_Loc_MF = reshape(Alpha_Loc_MF,N_Loc_Mode,3,3);
AlphaM_Loc_LF = zeros(size(AlphaM_Loc_MF));
for Loc_i = 1:N_Loc_Mode
    AlphaM_Loc_LF(Loc_i,:,:)  = R_MF_LF*squeeze(AlphaM_Loc_MF(Loc_i,:,:))/(R_MF_LF);
end
Alpha_Loc_LF = reshape(AlphaM_Loc_LF,N_Loc_Mode,9);

AlphaM_Ex_MF = reshape(Alpha_Ex_MF,N_Ex_Mode,3,3);
AlphaM_Ex_LF = zeros(size(AlphaM_Ex_MF));
for Ex_i = 1:N_Ex_Mode
    AlphaM_Ex_LF(Ex_i,:,:)  = R_MF_LF*squeeze(AlphaM_Ex_MF(Ex_i,:,:))/(R_MF_LF);
end
Alpha_Ex_LF = reshape(AlphaM_Ex_LF,N_Ex_Mode,9);

%% Retreive Axes from input figure with molecule plotted
hAx = findobj(hF,'type','axes');
hold on

    %% Plot local modes
    if Plot_Loc
        
        if Plot_EigenVec
            Mu_Loc_LF    = bsxfun(@times,   Mu_Loc_LF,EigVecM(:,EigneVec_Ind));
            Alpha_Loc_LF = bsxfun(@times,Alpha_Loc_LF,EigVecM(:,EigneVec_Ind));
        end
        
        Plot_Mu_Alpha(hAx,...
                      N_Loc_Mode,...
                      Center_Loc_LF,...
                      Mu_Loc_LF,...
                      Alpha_Loc_LF,...
                      Loc_Mode_colors,...
                      Plot_TDV,...
                      Scale_TDV,...
                      Plot_Raman,...
                      Scale_Raman,...
                      Normalize)
    end

    %% Plot Exciton modes
    if Plot_Ex
        
        Plot_Mu_Alpha(hAx,...
                      N_Ex_Mode,...
                      Center_Ex_LF,...
                      Mu_Ex_LF,...
                      Alpha_Ex_LF,...
                      Ex_Mode_colors,...
                      Plot_TDV,...
                      Scale_TDV,...
                      Plot_Raman,...
                      Scale_Raman,...
                      Normalize)
    end

    %% Plot Mixing coefficients
    if Plot_EigenVec
       Mix_Coeft  = EigVecM(:,EigneVec_Ind);
       [X0,Y0,Z0] = sphere(50);
       for k = 1:Structure.Num_Modes
           
           R_Scaling = 3;
           RR   = abs(Mix_Coeft(k)) .* R_Scaling;
           Sign = sign(Mix_Coeft(k));
           switch Sign
               case 1
                   F_Color = [1,0,0];
               case -1
                   F_Color = [0,0,1];
           end
           
           surf(hAx,...
                RR.*X0 + Center_Loc_LF(k,1),...
                RR.*Y0 + Center_Loc_LF(k,2),...
                RR.*Z0 + Center_Loc_LF(k,3),...
                'FaceColor',F_Color,...
                'FaceAlpha',0.5,...
                'LineStyle','none')
       end
    end
    
hold off

%% Figure setting
% Inherent the molecular plot title
Fig_Title = hAx.Title.String;

Mode_Ind_Str  = sprintf('#%d', Ex_Ind);
Mode_Freq_Str = sprintf(', @%6.2f cm^{-1}' ,OneDSFG.H.Sort_Ex_Freq(Ex_Ind +1));
Scaling_Str   = sprintf(', S\\mu= %2.1f, S\\alpha= %2.1f',Scale_TDV,Scale_Raman);

Fig_Title{length(Fig_Title)+1} = [Mode_Ind_Str, Mode_Freq_Str, Scaling_Str];
hAx.Title.String = Fig_Title;

lightangle(0,90)
camlight

view([-20,16])

%% Output
Output.Loc_Ind = Loc_Ind;
Output.Ex_Ind  = Ex_Ind;


function Plot_Mu_Alpha(hAx,N_Plot_Mode,Center,Mu,Alpha,Mode_colors,Plot_TDV,Scale_TDV,Plot_Raman,Scale_Raman,Normalize)
% Plot Transition dipoles
if Plot_TDV
    if Normalize
        % normalize to unit vector for direction comparison
        Mu_Loc_Int = sqrt(sum(Mu.^2,2));
        Mu = bsxfun(@rdivide,Mu,Mu_Loc_Int);
    end
    Mu_Loc_S = Scale_TDV .* Mu; % Scale TDV vector in plot 
    for j = 1: N_Plot_Mode
        quiver3(hAx,...
                Center(j,1),Center(j,2),Center(j,3),...
                Mu_Loc_S(j,1),Mu_Loc_S(j,2),Mu_Loc_S(j,3),0,...
                'LineWidth',2,...
                'Color',Mode_colors(j,:));
    end
end

% plot Raman tensors
if Plot_Raman
    N_mesh   = 20;
    if Normalize
        % normalize to unit vector for direction comparison
        Alpha_Norm = sqrt(sum(Alpha(:,:).^2,2)); % Norm defined in Silby's paper: JCP 1992, 97, 5607?5615.
        Alpha = bsxfun(@rdivide,Alpha,Alpha_Norm);
    end

    for i = 1: N_Plot_Mode
        RamanM = reshape(Alpha(i,:),3,3);
        plot_Raman(hAx,RamanM,Center(i,:),Scale_Raman,N_mesh,Mode_colors(i,:))
    end
end

