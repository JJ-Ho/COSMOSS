function Output = Update_Modes_Table(SpecType,Structure,COSMOSS_Inputs)
% This function update the table content of selected spectrum type

%% Defining column formats using ColumnFormat class
C_Index  = ColumnFormat('Index'  ,'short',40);
C_Freq   = ColumnFormat('Freq.'  ,'bank' ,50);
C_P_Num  = ColumnFormat('P. %'   ,'bank' ,50);
C_Mu_V   = ColumnFormat({'Mu x','Mu y','Mu z'},{'short','short','short'},{60,60,60});
C_Mu_Int = ColumnFormat('Mu Int' ,'bank' ,50);
C_A_V    = ColumnFormat({'A p1','A p2','A p3'},{'short','short','short'},{60,60,60});
C_A_Int  = ColumnFormat('Norm[A]','bank' ,50);
C_Norm1D = ColumnFormat('Norm 1D','bank' ,60);
C_Norm2D = ColumnFormat('Norm 2D','bank' ,60);

%% Calculate Spectral data
switch SpecType
    case 'FTIR'
        SD = FTIR_Main(Structure,COSMOSS_Inputs);
    case 'SFG'
        SD = OneDSFG_Main(Structure,COSMOSS_Inputs);
end

%% Common part
% Mode Frequency
Ex_Freq     = SD.H.Sort_Ex_Freq(2:end);
Num_Ex_Mode = length(Ex_Freq);
C_Freq.Data = Ex_Freq;

% mode index
Ex_Ind       = (1:Num_Ex_Mode)';
C_Index.Data = Ex_Ind;

% Participation number, percentage of local mode involve 
EigV         = SD.H.Sort_Ex_V(2:end,2:end);
P_Num        = 1./sum(EigV.^4,1)'./Num_Ex_Mode.*100;
C_P_Num.Data = P_Num;

% Transition dipole
Ex_Mu = squeeze(SD.Mu.Trans_Ex(1,2:end,:));
    % permute the matix dimension if only one mode
    if eq(Num_Ex_Mode,1)
        Ex_Mu = Ex_Mu';
    end
C_Mu_V.Data = Ex_Mu;

Ex_Mu_Int = sqrt(sum(Ex_Mu.^2,2));
C_Mu_Int.Data = Ex_Mu_Int;

%% Spectrum specific part and create corresponding TableColumn Class
switch SpecType
    case 'FTIR'
        ModeList = merge2cell(C_Index,...
                              C_Freq,...
                              C_P_Num,...
                              C_Mu_Int,...
                              C_Mu_V);
                 
    case 'SFG'
        Ex_Alpha = squeeze(SD.Alpha.Trans_Ex(1,2:end,:));
            % permute the matix dimension for spectial case
            if eq(Num_Ex_Mode,1)
                Ex_Alpha     = Ex_Alpha';
            end

        Ex_Alpha_Norm = sqrt(sum(Ex_Alpha(:,:).^2,2)); % Norm defined in Silby's paper: JCP 1992, 97, 5607?5615.
        C_A_Int.Data = Ex_Alpha_Norm;
        
        % Diagonalze Raman Tensor so I can look at their priciple values
        Ex_AlphaM = reshape(Ex_Alpha,Num_Ex_Mode,3,3);
        EigenV_Alpha = zeros(Num_Ex_Mode,3);
        for i = 1: Num_Ex_Mode
            [~,D] = eig(squeeze(Ex_AlphaM(i,:,:)));
            EigenV_Alpha(i,:) = diag(D)';
        end
        C_A_V.Data = EigenV_Alpha;
        
        
        Norm1D =  Ex_Mu_Int    .*Ex_Alpha_Norm;
        Norm2D = (Ex_Mu_Int.^3).*Ex_Alpha_Norm;
        C_Norm1D.Data = Norm1D;
        C_Norm2D.Data = Norm2D;
        
        % diaplay mode properties
        ModeList = merge2cell(C_Index,...
                              C_Freq,...
                              C_P_Num,...
                              C_Norm1D,...
                              C_Norm2D,...
                              C_Mu_Int,...
                              C_Mu_V,...
                              C_A_Int,...
                              C_A_V);
end


%% Output
Output.ModeList  = ModeList;
Output.SpecData  = SD;
         