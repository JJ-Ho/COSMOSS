function Output = Update_Modes_Figure(hF,GUI_Inputs, Structure, SpecData)
%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultMu_Alpha_Plot  = 1;
defaultMu_Alpha_Type  = 1;
defaultMu_Alpha_Ind   = [];
defaultTDV_Plot       = 1;
defaultTDV_Scale      = 1;
defaultRaman_Plot     = 1;
defaultRaman_Scale    = 1;
defaultRaman_Type     = 1;
defaultNormalize      = 1;
defaultPlot_EigVec    = 0;
defaultEigVec_Ind     = [];
defaultEigVec_Conv    = 1;
defaultEigVec_Scale   = 1;
defaultSpecType       = 1; % FTIR

% add options
addOptional(INPUT,'Mu_Alpha_Plot' , defaultMu_Alpha_Plot );
addOptional(INPUT,'Mu_Alpha_Type' , defaultMu_Alpha_Type );
addOptional(INPUT,'Mu_Alpha_Ind'  , defaultMu_Alpha_Ind  );
addOptional(INPUT,'TDV_Plot'      , defaultTDV_Plot      );
addOptional(INPUT,'TDV_Scale'     , defaultTDV_Scale     );
addOptional(INPUT,'Raman_Plot'    , defaultRaman_Plot    );
addOptional(INPUT,'Raman_Scale'   , defaultRaman_Scale   );
addOptional(INPUT,'Raman_Type'    , defaultRaman_Type    );
addOptional(INPUT,'Normalize'     , defaultNormalize     );
addOptional(INPUT,'Plot_EigVec'   , defaultPlot_EigVec   );
addOptional(INPUT,'EigVec_Ind'    , defaultEigVec_Ind    );
addOptional(INPUT,'EigVec_Conv'   , defaultEigVec_Conv   );
addOptional(INPUT,'EigVec_Scale'  , defaultEigVec_Scale  );
addOptional(INPUT,'SpecType'      , defaultSpecType      );

parse(INPUT,GUI_Inputs_C{:});

% Reassign Variable names
Mu_Alpha_Plot  = INPUT.Results.Mu_Alpha_Plot ;
Mu_Alpha_Type  = INPUT.Results.Mu_Alpha_Type ;
Mu_Alpha_Ind   = INPUT.Results.Mu_Alpha_Ind  ;
TDV_Plot       = INPUT.Results.TDV_Plot      ;
TDV_Scale      = INPUT.Results.TDV_Scale     ;
Raman_Plot     = INPUT.Results.Raman_Plot    ;
Raman_Scale    = INPUT.Results.Raman_Scale   ;
Raman_Type     = INPUT.Results.Raman_Type    ;
Normalize      = INPUT.Results.Normalize     ;
Plot_EigVec    = INPUT.Results.Plot_EigVec   ;
EigVec_Ind     = INPUT.Results.EigVec_Ind    ;
EigVec_Conv    = INPUT.Results.EigVec_Conv   ;
EigVec_Scale   = INPUT.Results.EigVec_Scale  ;
SpecType       = INPUT.Results.SpecType      ;

%% Define useful parameters
N_Mode_Total  = Structure.Nmodes;
%Mode_Index    = (1:N_Mode_Total)+1;
N_Mode_Select = length(Mu_Alpha_Ind);

%% molecular frame (need to rename to ab frame)
% center
Center_Loc_MF = Structure.LocCenter;

% Transition dipole
Mu_Loc_MF = SpecData.Mu.M_Lo_01;
Mu_Ex_MF  = SpecData.Mu.M_Ex_01;

% Raman Tensor
if or(eq(SpecType,1),eq(SpecType,4)) % check if FTIR or 2DIR
    Alpha_Loc_MF = zeros(size(Mu_Loc_MF,1),9); % size [N x 9]
    Alpha_Ex_MF  = zeros(size(Mu_Ex_MF,1),9);  % size [N x 9]
else
    Alpha_Loc_MF = SpecData.Alpha.M_Lo_01;
    Alpha_Ex_MF  = SpecData.Alpha.M_Ex_01;
end


% Eigen vectors
EigVecM      = SpecData.H.Sort_Ex_V1;
EigVecM2     = EigVecM.^2;
Center_Ex_MF = EigVecM2*(Structure.LocCenter);

%% Retreive Axes from input figure with molecule plotted
hAx = findobj(hF,'type','axes');
hold on

% Plot TDV/Raman for modes
if Mu_Alpha_Plot
    Mode_Colors = bsxfun(@times,ones(N_Mode_Select,3),[255,128,0]./256);
    switch Mu_Alpha_Type
        case 1 % local mode
            Center = Center_Loc_MF(Mu_Alpha_Ind,:);
            Mu     =     Mu_Loc_MF(Mu_Alpha_Ind,:);
            Alpha  =  Alpha_Loc_MF(Mu_Alpha_Ind,:);            
        case 2 % Exciton mode
            Center = Center_Ex_MF(Mu_Alpha_Ind,:);
            Mu     =     Mu_Ex_MF(Mu_Alpha_Ind,:);
            Alpha  =  Alpha_Ex_MF(Mu_Alpha_Ind,:);
    end
    
    if TDV_Plot
       Plot_Mu(hAx,N_Mode_Select,Center,Mu,Mode_Colors,TDV_Scale,Normalize)
    end
    
    if Raman_Plot
       Plot_Alpha(hAx,N_Mode_Select,Center,Alpha,Mode_Colors,Raman_Scale,Raman_Type,Normalize)
    end
end

% Plot Mixing coefficients
if Plot_EigVec
    Mode_Colors = bsxfun(@times,ones(N_Mode_Total,3),[255,128,0]./256);
    Mix_Coeft   = EigVecM(:,EigVec_Ind).* EigVec_Scale;
    switch EigVec_Conv
        case 1 % just the coefficient
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
           
        case 2 % with TDV
            Mu_Loc_MF  = bsxfun(@times,Mu_Loc_MF,Mix_Coeft);
            Center     = Center_Loc_MF;
            Mu         = Mu_Loc_MF;
            
            Plot_Mu(hAx,N_Mode_Total,Center,Mu,Mode_Colors,TDV_Scale,Normalize)

        case 3 % with Raman tensor
            Alpha_Loc_MF = bsxfun(@times,Alpha_Loc_MF,Mix_Coeft);
            Center       = Center_Loc_MF;
            Alpha        = Alpha_Loc_MF;
            
            Plot_Alpha(hAx,N_Mode_Total,Center,Alpha,Mode_Colors,Raman_Scale,Raman_Type,Normalize)
    end
end

hold off

%% Figure setting
% Inherent the molecular plot title
Fig_Title = hAx.Title.String;

Mode_Ind_Str  = sprintf('#%d',Mu_Alpha_Ind);
Mode_Freq_Str = sprintf(', @%6.2f cm^{-1}' ,SpecData.H.Sort_Ex_Freq(Mu_Alpha_Ind +1));
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
Output.Mu_Alpha_Ind = Mu_Alpha_Ind;

function Plot_Mu(hAx,N_Plot_Mode,Center,Mu,Mode_colors,TDV_Scale,Normalize)
if Normalize
    % normalize to unit vector for direction comparison
    Mu_Int = sqrt(sum(Mu.^2,2));
    Mu = bsxfun(@rdivide,Mu,Mu_Int);
end
Mu_Loc_S = TDV_Scale .* Mu; % Scale TDV vector in plot 

for j = 1: N_Plot_Mode
    quiver3(hAx,...
            Center(j,1),Center(j,2),Center(j,3),...
            Mu_Loc_S(j,1),Mu_Loc_S(j,2),Mu_Loc_S(j,3),0,...
            'LineWidth',2,...
            'MaxHeadSize',0.6,...
            'Color',Mode_colors(j,:));
end

function Plot_Alpha(hAx,N_Plot_Mode,Center,Alpha,Mode_colors,Raman_Scale,Raman_Type,Normalize)
% plot Raman tensors

N_mesh   = 20;
if Normalize
    % normalize to unit vector for direction comparison
    Alpha_Norm = sqrt(sum(Alpha(:,:).^2,2)); % Norm defined in Silby's paper: JCP 1992, 97, 5607?5615.
    Alpha = bsxfun(@rdivide,Alpha,Alpha_Norm);
end

for i = 1: N_Plot_Mode
    RamanM = reshape(Alpha(i,:),3,3);
    plot_Raman(hAx,RamanM,Center(i,:),Raman_Scale,N_mesh,Mode_colors(i,:),Raman_Type)
end
