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
Ex_Mu_Z     = Ex_Mu(:,3);
Ex_Mu_Int   = sqrt(sum(Ex_Mu.^2,2));

Ex_Alpha    = squeeze(OneDSFG.Alpha.Trans_Ex(1,2:end,:));
Ex_Alpha_ZZ = Ex_Alpha(:,9);
Ex_Alpha_Tr = sum(abs(Ex_Alpha(:,[1,5,9])),2); % take trace of abosolute value!

Sig_Z_1D     =  Ex_Mu_Z    .*Ex_Alpha_ZZ;
Sig_Z_2D     = (Ex_Mu_Z.^3).*Ex_Alpha_ZZ;

Norm_1D    =  Ex_Mu_Int    .*Ex_Alpha_Tr;
Norm_2D    = (Ex_Mu_Int.^3).*Ex_Alpha_Tr;

% diaplay mode properties
ModeList = [ Ex_Ind,...
             Ex_Freq,...
             Norm_1D,...
             Norm_2D,...
             Ex_Mu_Int,...
             Ex_Alpha_Tr,...
             Sig_Z_1D,...
             Sig_Z_2D,...
             Ex_Mu_Z,...
             Ex_Alpha_ZZ,...
             ];

%% Output
Output.ModeList  = ModeList;
Output.OneDSFG   = OneDSFG;
         