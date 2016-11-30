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

C_PumpF  = ColumnFormat('Pump' ,'bank' ,55);
C_ProbF  = ColumnFormat('Probe','bank' ,55);
C_Int    = ColumnFormat('Intensity' ,'bank' ,70);

C_Path_Label = ColumnFormat('Path' ,'short' ,30);

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
    PathName = {'R1','R2','R3','NR1','NR2','NR3'};
    Pump_F = [];
    Prob_F = [];
    Ex_Ind = [];
    Int    = [];
    PathName_Str = [];
    
    for L = 1:length(PathName)
        % Pump Frequency
        New_Pump_F = abs(SD.Freq.(PathName{L})(:,1));
        Pump_F = [Pump_F ; New_Pump_F ];
        N_Path = length(New_Pump_F);

        % Probe Frequency
        New_Prob_F = abs(SD.Freq.(PathName{L})(:,3));
        Prob_F = [Prob_F ; New_Prob_F ];    

        % mode index
        New_Ind = (1:N_Path)';
        Ex_Ind = [Ex_Ind ; New_Ind ];    

        % signal intensity
        New_Int = SD.EJRBeta.(PathName{L})';
        Int = [Int ; New_Int ];
        
        % Pathway label (Number version)
        Tmp = ones(N_Path,1).*L;
        PathName_Str = [PathName_Str ; Tmp];
    end
    
    % Apply intensity cutoff
    CutOff_R = 0.01;
    CutOff_I = Int < max(abs(Int)*CutOff_R);
    
    Int(CutOff_I) = [];
    Pump_F(CutOff_I) = [];
    Prob_F(CutOff_I) = [];
    Ex_Ind(CutOff_I) = [];
    PathName_Str(CutOff_I) = [];
    
    C_PumpF.Data = Pump_F;
    C_ProbF.Data = Prob_F;
    C_Index.Data = Ex_Ind;
    C_Int.Data = Int;
    C_Path_Label.Data = PathName_Str;
end

%% Spectrum specific part and create corresponding TableColumn Class
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
                              C_Int,...
                              C_Path_Label);
                              
end


%% Output
Output.ModeList  = ModeList;
Output.SpecData  = SD;
         