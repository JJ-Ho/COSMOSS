function Output = Update_Modes_Figure(GUI_Inputs, Structure, SpecData)
%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultMode_Type       = 'Exciton Modes';
defaultMode_Ind        = [];
defaultTDV_Plot        = 1;
defaultTDV_Scale       = 1;
defaultRaman_Plot      = 1;
defaultRaman_Scale     = 1;
defaultRaman_Type      = 1;
defaultTDV_Normalize   = 1;
defaultRaman_Normalize = 1;
defaultPlot_EigVec     = 0;
defaultEigVec_Conv     = 1;
defaultEigVec_Scale    = 1;
defaultSpecType        = 1; % FTIR

% add options
addOptional(INPUT,'Mode_Type'      , defaultMode_Type      );
addOptional(INPUT,'Mode_Ind'       , defaultMode_Ind       );
addOptional(INPUT,'TDV_Plot'       , defaultTDV_Plot       );
addOptional(INPUT,'TDV_Scale'      , defaultTDV_Scale      );
addOptional(INPUT,'Raman_Plot'     , defaultRaman_Plot     );
addOptional(INPUT,'Raman_Scale'    , defaultRaman_Scale    );
addOptional(INPUT,'Raman_Type'     , defaultRaman_Type     );
addOptional(INPUT,'TDV_Normalize'  , defaultTDV_Normalize  );
addOptional(INPUT,'Raman_Normalize', defaultRaman_Normalize);
addOptional(INPUT,'Plot_EigVec'    , defaultPlot_EigVec    );
addOptional(INPUT,'EigVec_Conv'    , defaultEigVec_Conv    );
addOptional(INPUT,'EigVec_Scale'   , defaultEigVec_Scale   );
addOptional(INPUT,'SpecType'       , defaultSpecType       );

parse(INPUT,GUI_Inputs_C{:});

% Reassign Variable names
Mode_Type       = INPUT.Results.Mode_Type      ;
Mode_Ind        = INPUT.Results.Mode_Ind       ;
TDV_Plot        = INPUT.Results.TDV_Plot       ;
TDV_Scale       = INPUT.Results.TDV_Scale      ;
Raman_Plot      = INPUT.Results.Raman_Plot     ;
Raman_Scale     = INPUT.Results.Raman_Scale    ;
Raman_Type      = INPUT.Results.Raman_Type     ;
TDV_Normalize   = INPUT.Results.TDV_Normalize  ;
Raman_Normalize = INPUT.Results.Raman_Normalize;
Plot_EigVec     = INPUT.Results.Plot_EigVec    ;
EigVec_Conv     = INPUT.Results.EigVec_Conv    ;
EigVec_Scale    = INPUT.Results.EigVec_Scale   ;
SpecType        = INPUT.Results.SpecType       ;

%% Draw Molecule
hF = Structure.Draw;

%% Define useful parameters
N_Mode_Total  = Structure.Nmodes;
N_Mode_Select = length(Mode_Ind);

%% molecular frame (need to rename to ab frame)
% center
Center_Loc_MF = Structure.LocCenter;

% Transition dipole
Mu_Loc_MF = SpecData.Mu.M_Lo_01;
Mu_Ex_MF  = SpecData.Mu.M_Ex_01;

% Raman Tensor
switch SpecType
    case 'FTIR'
        Alpha_Loc_MF = zeros(size(Mu_Loc_MF,1),9); % size [N x 9]
        Alpha_Ex_MF  = zeros(size(Mu_Ex_MF,1),9);  % size [N x 9]
    case 'SFG'
        Alpha_Loc_MF = SpecData.Alpha.M_Lo_01;
        Alpha_Ex_MF  = SpecData.Alpha.M_Ex_01;
end

%% Plot TDV/Raman
hAx = findobj(hF,'type','axes');
hold on

% Plot TDV/Raman for modes
if or(TDV_Plot,Raman_Plot)
    Mode_Colors = bsxfun(@times,ones(N_Mode_Select,3),[255,128,0]./256);
    switch Mode_Type
        case 'Local modes'
            Center = Center_Loc_MF(Mode_Ind,:);
            Mu     =     Mu_Loc_MF(Mode_Ind,:);
            Alpha  =  Alpha_Loc_MF(Mode_Ind,:);            
        case 'Exciton modes'
            Center = repmat(Structure.CoM,N_Mode_Select,1);
            Mu     =     Mu_Ex_MF(Mode_Ind,:);
            Alpha  =  Alpha_Ex_MF(Mode_Ind,:);
    end
    
    if TDV_Plot
        if TDV_Normalize % normalize to unit vector for direction comparison
            Mu_Int = sqrt(sum(Mu.^2,2));
            Mu = bsxfun(@rdivide,Mu,Mu_Int);
        end
        Mu_S = TDV_Scale .* Mu; % Scale TDV vector in plot  
        
        Plot_Mu(hAx,Center,Mu_S,Mode_Colors)
    end
    
    if Raman_Plot
        if Raman_Normalize % normalize to unit vector for direction comparison
            Alpha_Norm = sqrt(sum(Alpha(:,:).^2,2)); % Norm defined in Silby's paper: JCP 1992, 97, 5607?5615.
            Alpha = bsxfun(@rdivide,Alpha,Alpha_Norm);
        end
        Alpha_S = Raman_Scale .* Alpha;
          
        Plot_Alpha(hAx,Center,Alpha_S,Mode_Colors,Raman_Type)
    end
end

%% Plot Mixing coefficients
% Eigen vectors
EigVecM = SpecData.H.Sort_Ex_V1;
if Plot_EigVec
    if gt(length(Mode_Ind),1)
        disp('Eigenvevtor visulization only take 1 mode at a time')
        Mode_Ind = Mode_Ind(1);
        disp(['Drawing mode #',num2str(Mode_Ind),' only...'])
    end
    Mode_Colors = bsxfun(@times,ones(N_Mode_Total,3),[255,128,0]./256);
    Mix_Coeft   = EigVecM(:,Mode_Ind).* EigVec_Scale;
    switch EigVec_Conv
        case 'Sphere' % just the coefficient
           [X0,Y0,Z0] = sphere(50);
           for k = 1:Structure.Nmodes
               
               RR   = abs(Mix_Coeft(k));
               Sign = sign(Mix_Coeft(k));
               switch Sign
                   case 1
                       F_Color = [1,0,0];
                   case -1
                       F_Color = [0,0,1];
               end

               surf(hAx,...
                    RR.*X0 + Center_Loc_MF(k,1),...
                    RR.*Y0 + Center_Loc_MF(k,2),...
                    RR.*Z0 + Center_Loc_MF(k,3),...
                    'FaceColor',F_Color,...
                    'FaceAlpha',0.5,...
                    'LineStyle','none')
           end
           
        case 'TDV' % with TDV
            Mu_Loc_MF  = bsxfun(@times,Mu_Loc_MF,Mix_Coeft);
            Center     = Center_Loc_MF;
            Mu         = Mu_Loc_MF;
            Mu_S = TDV_Scale .* Mu; % Scale TDV vector in plot 
            
            Plot_Mu(hAx,Center,Mu_S,Mode_Colors)

        case 'Raman' % with Raman tensor
            Alpha_Loc_MF = bsxfun(@times,Alpha_Loc_MF,Mix_Coeft);
            Center       = Center_Loc_MF;
            Alpha        = Alpha_Loc_MF;
            Alpha_S = Raman_Scale .* Alpha;
            
            Plot_Alpha(hAx,Center,Alpha_S,Mode_Colors,Raman_Type)
    end
end

hold off

%% Figure setting
% Inherent the molecular plot title
Fig_Title = hAx.Title.String;

Mode_Ind_Str  = sprintf('#%d',Mode_Ind);
Mode_Freq_Str = sprintf(', @%6.2f cm^{-1}' ,SpecData.H.Sort_Ex_Freq(Mode_Ind +1));
Scaling_Str   = sprintf(', S\\mu= %2.1f, S\\alpha= %2.1f',TDV_Scale,Raman_Scale);

Fig_Title{length(Fig_Title)+1} = [Mode_Ind_Str, Mode_Freq_Str, Scaling_Str];
hAx.Title.String = Fig_Title;

% lightangle(0,90)
% camlight

view([-20,16])

% deal with box size
axis tight;
ExtScal = 0.2;
hAx.XLim = hAx.XLim + ExtScal*sum(hAx.XLim.*[-1,1])*[-1,1];
hAx.YLim = hAx.YLim + ExtScal*sum(hAx.YLim.*[-1,1])*[-1,1];
hAx.ZLim = hAx.ZLim + ExtScal*sum(hAx.ZLim.*[-1,1])*[-1,1];
hAx.DataAspectRatio = [1,1,1];

%% Output
Output.Mu_Alpha_Ind = Mode_Ind;

function Plot_Mu(hAx,Center,Mu_S,Mode_colors)
N_Plot_Mode = size(Center,1);
for j = 1: N_Plot_Mode
    quiver3(hAx,...
            Center(j,1),Center(j,2),Center(j,3),...
            Mu_S(j,1),Mu_S(j,2),Mu_S(j,3),0,...
            'LineWidth',2,...
            'MaxHeadSize',0.6,...
            'Color',Mode_colors(j,:));
end

function Plot_Alpha(hAx,Center,Alpha,Mode_colors,Raman_Type)
N_Plot_Mode = size(Center,1);
N_mesh = 20;
for i = 1: N_Plot_Mode
    RamanM = reshape(Alpha(i,:),3,3);
    plot_Raman(hAx,RamanM,Center(i,:),N_mesh,Mode_colors(i,:),Raman_Type)
end
