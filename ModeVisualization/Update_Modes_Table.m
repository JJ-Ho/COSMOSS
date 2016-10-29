function Output = Update_Modes_Table(Structure, MainGUI_Inputs)
% Run OneDSFG to get the corresponding mu and alpha of exciton modes
% Retrieve Label index and Coupling model from COSMOSS GUI if any, if
% running Plot_Exciton stand alone for debug testing, give a field 'debug'
% to use the default values in OneDSFG_Main.m

OneDSFG = OneDSFG_Main(Structure,MainGUI_Inputs);

Ex_Freq     = OneDSFG.H.Sort_Ex_Freq(2:end);
Num_Ex_Mode = length(Ex_Freq);
Ex_Ind      = (1:Num_Ex_Mode)';
Ex_Mu       = squeeze(OneDSFG.Mu.Trans_Ex(1,2:end,:));
% permute the matix dimension for spectial case
if eq(Num_Ex_Mode,1)
    Ex_Mu     = Ex_Mu';
end

Ex_Mu_x     = Ex_Mu(:,1);
Ex_Mu_y     = Ex_Mu(:,2);
Ex_Mu_z     = Ex_Mu(:,3);
Ex_Mu_Int   = sqrt(sum(Ex_Mu.^2,2));

Ex_Alpha      = squeeze(OneDSFG.Alpha.Trans_Ex(1,2:end,:));
% permute the matix dimension for spectial case
if eq(Num_Ex_Mode,1)
    Ex_Alpha     = Ex_Alpha';
end

Ex_Alpha_zz   = Ex_Alpha(:,9);
Ex_Alpha_Norm = sqrt(sum(Ex_Alpha(:,:).^2,2)); % Norm defined in Silby's paper: JCP 1992, 97, 5607?5615.

% Diagonalze Raman Tensor so I can look at their priciple values
Ex_AlphaM = reshape(Ex_Alpha,Num_Ex_Mode,3,3);
EigenV_Alpha = zeros(Num_Ex_Mode,3);
for i = 1: Num_Ex_Mode
    [~,D] = eig(squeeze(Ex_AlphaM(i,:,:)));
    EigenV_Alpha(i,:) = diag(D)';
end

Norm_1D    =  Ex_Mu_Int    .*Ex_Alpha_Norm;
Norm_2D    = (Ex_Mu_Int.^3).*Ex_Alpha_Norm;

% diaplay mode properties
ModeList = [ Ex_Ind,...
             Ex_Freq,...
             Norm_1D,...
             Norm_2D,...
             Ex_Mu_Int,...
             Ex_Mu_x,...
             Ex_Mu_y,...
             Ex_Mu_z,...
             Ex_Alpha_Norm,...
             EigenV_Alpha(:,1),...
             EigenV_Alpha(:,2),...
             EigenV_Alpha(:,3),...
             ];

%% Output
Output.ModeList  = ModeList;
Output.OneDSFG   = OneDSFG;
         