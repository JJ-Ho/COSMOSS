function Output = Update_Modes_Table(SpecType,Structure,COSMOSS_Inputs)
% This function update the table content of selected spectrum type

%% Defining column formats using ColumnFormat class
C_Index  = ColumnFormat('Index'  ,'short',40);
C_Freq   = ColumnFormat('Freq.'  ,'bank' ,50);
C_P_Num  = ColumnFormat('P. %'   ,'bank' ,40);
C_Mu_V   = ColumnFormat({'Mu x','Mu y','Mu z'},{'bank','bank','bank'},{40,40,40});
C_Mu_Int = ColumnFormat('Mu Int' ,'bank' ,50);
C_A_V    = ColumnFormat({'A p1','A p2','A p3'},{'bank','bank','bank'},{40,40,40});
C_A_Int  = ColumnFormat('Norm[A]','bank' ,60);
C_Norm1D = ColumnFormat('Norm 1D','bank' ,60);
C_Norm2D = ColumnFormat('Norm 2D','bank' ,60);

C_PumpF  = ColumnFormat('Pump F.' ,'bank' ,55);
C_ProbF  = ColumnFormat('Probe F.','bank' ,55);
C_2D_Int = ColumnFormat('2D Int.' ,'bank' ,70);

%% Calculate Spectral data
switch SpecType
    case 'FTIR'
        SD = FTIR_Main(Structure,COSMOSS_Inputs);
    case 'SFG'
        SD = OneDSFG_Main(Structure,COSMOSS_Inputs);
    case 'TwoDIR'
        [~,SD] = TwoDIR_Main(Structure,COSMOSS_Inputs);
end

%% 1D spectra Common part 
if or(strcmp(SpecType,'FTIR'),strcmp(SpecType,'SFG'))
    % Mode Frequency
    Pump_F     = SD.H.Sort_Ex_F1;
    Num_Ex_Mode = SD.H.Num_Modes;
    C_Freq.Data = Pump_F;

    % mode index
    Ex_Ind       = (1:Num_Ex_Mode)';
    C_Index.Data = Ex_Ind;

    % Participation number, percentage of local mode involve 
    EigV         = SD.H.Sort_Ex_V1;
    P_Num        = 1./sum(EigV.^4,1)'./Num_Ex_Mode.*100;
    C_P_Num.Data = P_Num;

    % Transition dipole
    Ex_Mu = SD.Mu.M_Ex_01;
        % permute the matix dimension if only one mode
        if eq(Num_Ex_Mode,1)
            Ex_Mu = Ex_Mu';
        end
    C_Mu_V.Data = Ex_Mu;

    Ex_Mu_Int = sqrt(sum(Ex_Mu.^2,2));
    C_Mu_Int.Data = Ex_Mu_Int;
end

%% 2D spectrum Common part
if or(strcmp(SpecType,'TwoDIR'),strcmp(SpecType,'TwoDSFG'))
    % Pump Frequency
    Pump_F = abs(SD.Freq.R1(:,1));
    N_Path = length(Pump_F);
    C_PumpF.Data = Pump_F;
    
    % Probe Frequency
    Prob_F = abs(SD.Freq.R1(:,3));
    C_ProbF.Data = Prob_F;

    % mode index
    Ex_Ind       = (1:N_Path)';
    C_Index.Data = Ex_Ind;
    
    % signal intensity
    Int = SD.EJRBeta.R1';
    C_2D_Int.Data = Int;
end

%% 1D Spectrum specific part and create corresponding TableColumn Class
switch SpecType
    case 'FTIR'
        ModeList = merge2cell(C_Index,...
                              C_Freq,...
                              C_P_Num,...
                              C_Mu_Int,...
                              C_Mu_V);
                 
    case 'SFG'
        Ex_Alpha = SD.Alpha.M_Ex_01;
            % permute the matix dimension for spectial case
            if eq(Num_Ex_Mode,1)
                Ex_Alpha = Ex_Alpha';
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
                              
    case 'TwoDIR'
        ModeList = merge2cell(C_Index,...
                              C_PumpF,...
                              C_ProbF,...
                              C_2D_Int);
end


%% Output
Output.ModeList  = ModeList;
Output.SpecData  = SD;
         