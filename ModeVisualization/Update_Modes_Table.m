function Output = Update_Modes_Table(SpecData)
% This function update the table content of selected spectrum type

SpecType = SpecData.SpecType;
switch SpecType
    case 'FTIR'
        Ex_Freq     = SpecData.H.Sort_Ex_Freq(2:end);
        Num_Ex_Mode = length(Ex_Freq);
        Ex_Ind      = (1:Num_Ex_Mode)';
        Ex_Mu       = squeeze(SpecData.Mu.Trans_Ex(1,2:end,:));
        % permute the matix dimension if only one mode
        if eq(Num_Ex_Mode,1)
            Ex_Mu     = Ex_Mu';
        end
        
        Ex_Mu_x     = Ex_Mu(:,1);
        Ex_Mu_y     = Ex_Mu(:,2);
        Ex_Mu_z     = Ex_Mu(:,3);
        Ex_Mu_Int   = sqrt(sum(Ex_Mu.^2,2));

        ModeList = [ Ex_Ind,...
                     Ex_Freq,...
                     Ex_Mu_Int,...
                     Ex_Mu_x,...
                     Ex_Mu_y,...
                     Ex_Mu_z,...
                     ];
                 
    case 'SFG'
        Ex_Freq     = SpecData.H.Sort_Ex_Freq(2:end);
        Num_Ex_Mode = length(Ex_Freq);
        Ex_Ind      = (1:length(Ex_Freq))';
        Ex_Mu       = squeeze(SpecData.Mu.Trans_Ex(1,2:end,:));
        % permute the matix dimension if only one mode
        if eq(Num_Ex_Mode,1)
            Ex_Mu     = Ex_Mu';
        end

        Ex_Mu_x     = Ex_Mu(:,1);
        Ex_Mu_y     = Ex_Mu(:,2);
        Ex_Mu_z     = Ex_Mu(:,3);
        Ex_Mu_Int   = sqrt(sum(Ex_Mu.^2,2));

        Ex_Alpha      = squeeze(SpecData.Alpha.Trans_Ex(1,2:end,:));
        % permute the matix dimension for spectial case
        if eq(Num_Ex_Mode,1)
            Ex_Alpha     = Ex_Alpha';
        end

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
end


%% Output
Output.ModeList  = ModeList;
Output.OneDSFG   = SpecData;
         