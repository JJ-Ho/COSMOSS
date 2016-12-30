function Output = Update_Modes_Table(SpecType,Structure,COSMOSS_Inputs)
% This function update the table content of selected spectrum type

%% Defining column formats using ColumnFormat class
C_Index  = ColumnFormat('Index'  ,'short',40);
C_Freq   = ColumnFormat('Freq.'  ,'bank' ,50);
C_P_Num  = ColumnFormat('P. %'   ,'bank' ,40);
C_Mu_V   = ColumnFormat({'Mu x%','Mu y%','Mu z%'},{'bank','bank','bank'},{40,40,40});
C_A_V    = ColumnFormat({'A p1','A p2','A p3'},{'bank','bank','bank'},{40,40,40});
C_Mu_Int = ColumnFormat('N[Mu]' ,'bank' ,40);
C_A_Int  = ColumnFormat('N[A]'  ,'bank' ,40);
C_Norm1D = ColumnFormat('N[1D]' ,'bank' ,40);
C_Sig    = ColumnFormat('Signal','bank' ,40);

C_PumpF  = ColumnFormat('Pump' ,'bank' ,55);
C_ProbF  = ColumnFormat('Probe','bank' ,55);
C_Int    = ColumnFormat('Intensity' ,'bank' ,70);

C_PathName = ColumnFormat('Path' ,'short' ,30);
C_PathInd  = ColumnFormat({'4th','3rd','2nd','1st'},{'short','short','short'},{40,40,40,40});

%% Calculate Spectral data
switch SpecType
    case 'FTIR'
        SD = FTIR_Main(Structure,COSMOSS_Inputs);
    case 'SFG'
        SD = OneDSFG_Main(Structure,COSMOSS_Inputs);
    case 'TwoDIR'
        [~,SD] = TwoDIR_Main(Structure,COSMOSS_Inputs);
    case 'TwoDSFG'
        [~,SD] = TwoDSFG_Main_Sparse(Structure,COSMOSS_Inputs);
end

%% 1D spectra Common part 
PathType = 'None';
if or(strcmp(SpecType,'FTIR'),strcmp(SpecType,'SFG'))
    % Mode Frequency
    Pump_F = SD.H.Sort_Ex_F1;
    C_Freq = ImportSortInd(C_Freq,Pump_F);

    % mode index
    Num_Ex_Mode = SD.H.Num_Modes;
    Ex_Ind      = (1:Num_Ex_Mode)';
    C_Index = ImportSortInd(C_Index,Ex_Ind);
    
    % Participation number, percentage of local mode involve 
    EigV         = SD.H.Sort_Ex_V1;
    P_Num        = 1./sum(EigV.^4,1)'./Num_Ex_Mode.*100;
    C_P_Num = ImportSortInd(C_P_Num,P_Num);

    % Transition dipole vector direction
    Ex_Mu = SD.Mu.M_Ex_01; % [Nx3]
    Ex_Mu_N = bsxfun(@rdivide,Ex_Mu,sqrt(sum(Ex_Mu.^2,2))); % normalized vector components
    C_Mu_V = ImportSortInd(C_Mu_V,Ex_Mu_N.^2); % (x^2,y^2,z^2)./(x^2+y^2+z^2) so when the sum of the three components = 1

    % Transition dipole Intensity
    Ex_Mu_Int = sqrt(sum(Ex_Mu.^2,2));
    C_Mu_Int = ImportSortInd(C_Mu_Int,Ex_Mu_Int);
end

%% 2D spectrum Common part
if or(strcmp(SpecType,'TwoDIR'),strcmp(SpecType,'TwoDSFG'))
    PathType = {'R1','R2','R3','NR1','NR2','NR3'}; % [Improve] add GUI input for this 
    Pump_F       = [];
    Prob_F       = [];
    Int          = [];
    PathName_N   = [];
    PathInd_N    = [];
    
    % unify path index from type+index to accumulated index
    % this will have size difference problem when PathType no longer
    % including all. Need to fix the cutOff Index size issue later
    PathType_Total = fieldnames(SD.Index);
    N_Path_Total = 0;
    for P = 1:length(PathType_Total)
        N_Path_Total = N_Path_Total + size(SD.Index.(PathType_Total{P}),1);
    end
    Index = (1:N_Path_Total)';
    
    for L = 1:length(PathType)
        % Pump Frequency
        New_Pump_F = abs(SD.Freq.(PathType{L})(:,1));
        Pump_F = [Pump_F ; New_Pump_F ];
        N_Path = length(New_Pump_F);

        % Probe Frequency
        New_Prob_F = abs(SD.Freq.(PathType{L})(:,3));
        Prob_F = [Prob_F ; New_Prob_F ];    

        % signal intensity
        New_Int = SD.Int.(PathType{L});
        Int = [Int ; New_Int(:) ];
        
        % Pathway Name (Number version)
        Tmp = ones(N_Path,1).*L;
        PathName_N = [PathName_N ; Tmp];
        
        % Pathways index (Number version)
        New_PathInd_N = SD.Index.(PathType{L});
        PathInd_N = [PathInd_N ; New_PathInd_N];     
    end
    
    % Apply intensity cutoff
    CutOff_R = 1E-2; % [Improve] add GUI input for this 
    CutOff_I = Int < max(abs(Int)*CutOff_R);
    
    Int(CutOff_I)           = [];
    Pump_F(CutOff_I)        = [];
    Prob_F(CutOff_I)        = [];
    Index(CutOff_I)         = [];
    PathName_N(CutOff_I)    = [];
    PathInd_N(CutOff_I,:)   = [];
    
    % save data to formated column
    C_PumpF     = ImportSortInd(C_PumpF    ,Pump_F);
    C_ProbF     = ImportSortInd(C_ProbF    ,Prob_F);
    C_Index     = ImportSortInd(C_Index    ,Index);
    C_Int       = ImportSortInd(C_Int      ,Int);
    C_PathName  = ImportSortInd(C_PathName ,PathName_N);
    C_PathInd   = ImportSortInd(C_PathInd  ,PathInd_N);
    
    % Replace face value of Data to proper stings
    % Pathway Name (String version)
    FormattedPathNameStr = strcat('<html><pre style="color:black ; font-weight: bold;">', ...
                                   strcat({' '},PathType(PathName_N)'), ...
                                   '</pre></html>');   
    C_PathName.Data = FormattedPathNameStr;

    % Pathways index (String version)
    FormattedPathInd = TableFromatPathways(PathName_N,PathInd_N);        
    C_PathInd.Data   = FormattedPathInd;
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
        Ex_Alpha = SD.Alpha.M_Ex_01; % [Nx9]
        
        % deal with norm(alpha)
        Ex_Alpha_Norm = sqrt(sum(Ex_Alpha(:,:).^2,2)); % Norm defined in Silby's paper: JCP 1992, 97, 5607?5615.
        C_A_Int = ImportSortInd(C_A_Int,Ex_Alpha_Norm);
        
        % deal with norm(1DSFG)
        Norm1D =  Ex_Mu_Int.*Ex_Alpha_Norm;
        C_Norm1D = ImportSortInd(C_Norm1D,Norm1D);
        
        % Diagonalze Raman Tensor so I can look at their principle values
        % Also project transition dipoles onto the principle axis of Raman tensors 
        Ex_AlphaM = reshape(Ex_Alpha,Num_Ex_Mode,3,3);
        EigenV_Alpha = zeros(Num_Ex_Mode,3);
        Ex_Mu_N_P_Raman = zeros(Num_Ex_Mode,3);
        for i = 1: Num_Ex_Mode
            [V,D] = eig(squeeze(Ex_AlphaM(i,:,:)));
            EigenV_Alpha(i,:) = diag(D)';
            EigenV_Alpha(i,:) = EigenV_Alpha(i,:)./norm(EigenV_Alpha(i,:)); % normalize
            
            Ex_Mu_N_P_Raman(i,:) = Ex_Mu_N(i,:) * V;
        end
        C_A_V = ImportSortInd(C_A_V,EigenV_Alpha.^2); % (x^2,y^2,z^2)./(x^2+y^2+z^2) so when the sum of the three components = 1
        
        % update the Principle axis projected transtion dipole components
        C_Mu_V.Name = {'Mu p1','Mu p2','Mu p3'};
        C_Mu_V = ImportSortInd(C_Mu_V,Ex_Mu_N_P_Raman.^2); % (x^2,y^2,z^2)./(x^2+y^2+z^2) so when the sum of the three components = 1
        
        % signal size
        Signal = SD.EJLabFrame';
        C_Sig = ImportSortInd(C_Sig,Signal);
        
        % diaplay mode properties
        ModeList = merge2cell(C_Index,...
                              C_Freq,...
                              C_P_Num,...
                              C_Sig,...
                              C_Norm1D,...
                              C_Mu_Int,...
                              C_Mu_V,...
                              C_A_Int,...
                              C_A_V);
                              
    case 'TwoDIR'
        ModeList = merge2cell(C_Index,...
                              C_PumpF,...
                              C_ProbF,...
                              C_Int,...
                              C_PathName,...
                              C_PathInd);
                          
    case 'TwoDSFG'
        ModeList = merge2cell(C_Index,...
                              C_PumpF,...
                              C_ProbF,...
                              C_Int,...
                              C_PathName,...
                              C_PathInd);
                              
end

%% Output
Output.ModeList  = ModeList;
Output.SpecData  = SD;
Output.PathType  = PathType;