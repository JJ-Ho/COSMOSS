function Output = Update_Modes_Figure(hF,GUI_Inputs, Structure, OneDSFG)
%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
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


%% Calculate the exciton center
Center_Loc = Structure.center;

EigVecM   = OneDSFG.H.Sort_Ex_V(2:end,2:end); % get ride of ground state
EigVecM2  = EigVecM.^2; 
Center_Ex = EigVecM2 * Center_Loc;

%% Retreive Axes from input figure with molecule plotted 
hAx = findobj(hF,'type','axes');
hold on

    %% Plot local modes
    if Plot_Loc
        
        N_Plot_Mode = length(Loc_Ind);
        Center_Loc  = Center_Loc(Loc_Ind,:);
        Mu_Loc      = squeeze(OneDSFG.Mu.   Trans_Loc(1,Loc_Ind+1,:)); % shift by 1 to avoid ground state
        Alpha_Loc   = squeeze(OneDSFG.Alpha.Trans_Loc(1,Loc_Ind+1,:));

        if Plot_EigenVec
            Mu_Loc    = bsxfun(@times,   Mu_Loc,EigVecM(:,EigneVec_Ind));
            Alpha_Loc = bsxfun(@times,Alpha_Loc,EigVecM(:,EigneVec_Ind));
        end
        
        Plot_Mu_Alpha(hAx,...
                      N_Plot_Mode,...
                      Center_Loc,...
                      Mu_Loc,...
                      Alpha_Loc,...
                      Loc_Mode_colors,...
                      Plot_TDV,...
                      Scale_TDV,...
                      Plot_Raman,...
                      Scale_Raman,...
                      Normalize)
    end

    %% Plot Exciton modes
    if Plot_Ex
        N_Plot_Mode = length(Ex_Ind);
        Center_Ex   = Center_Ex(Ex_Ind,:);
        Mu_Ex       = squeeze(OneDSFG.Mu.   Trans_Ex(1,Ex_Ind+1,:)); % shift by 1 to avoid ground state
        Alpha_Ex    = squeeze(OneDSFG.Alpha.Trans_Ex(1,Ex_Ind+1,:));
        
        Plot_Mu_Alpha(hAx,...
                      N_Plot_Mode,...
                      Center_Ex,...
                      Mu_Ex,...
                      Alpha_Ex,...
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
                RR.*X0 + Center_Loc(k,1),...
                RR.*Y0 + Center_Loc(k,2),...
                RR.*Z0 + Center_Loc(k,3),...
                'FaceColor',F_Color,...
                'FaceAlpha',0.5,...
                'LineStyle','none')
       end
    end
    
hold off

%% Figure setting
Fig_Title = ['Loc #: ', num2str(Loc_Ind ), '; Ex #: ', num2str(Ex_Ind ) ];
hAx.Title.String = Fig_Title;
lightangle(-45,30)

%% Output
Output.Loc_Ind = Loc_Ind;
Output.Ex_Ind  = Ex_Ind;


function Plot_Mu_Alpha(hAx,N_Plot_Mode,Center,Mu,Alpha,Mode_colors,Plot_TDV,Scale_TDV,Plot_Raman,Scale_Raman,Normalize)
% permute the matix dimension for spectial case
if N_Plot_Mode == 1
    Mu    = Mu';
    Alpha = Alpha';
end

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
        Alpha_Loc_Tr = sum(abs(Alpha(:,[1,5,9])),2);
        Alpha = bsxfun(@rdivide,Alpha,Alpha_Loc_Tr);
    end

    for i = 1: N_Plot_Mode
        RamanM_Loc = reshape(Alpha(i,:),3,3);
        plot_Raman(hAx,RamanM_Loc,Center(i,:),Scale_Raman,N_mesh,Mode_colors(i,:))
    end
end

