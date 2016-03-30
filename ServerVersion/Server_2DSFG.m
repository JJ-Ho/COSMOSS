function [SpectraGrid,Response] = Server_2DSFG(Structure,GUI_Inputs)
%% Calculate TwoD response

Sample_Num = GUI_Inputs.Sample_Num;
TSTART_whole = tic;
if eq(GUI_Inputs.Sampling,1)
    % Pre-allocate
    GridSize     = length(GUI_Inputs.FreqRange);

    Rephasing    = zeros(GridSize,GridSize,Sample_Num);
    NonRephasing = zeros(GridSize,GridSize,Sample_Num);
    SpecR1       = zeros(GridSize,GridSize,Sample_Num);
    SpecR2       = zeros(GridSize,GridSize,Sample_Num);
    SpecR3       = zeros(GridSize,GridSize,Sample_Num);
    SpecNR1      = zeros(GridSize,GridSize,Sample_Num);
    SpecNR2      = zeros(GridSize,GridSize,Sample_Num);
    SpecNR3      = zeros(GridSize,GridSize,Sample_Num);
    
    Num_Modes = Structure.Num_Modes;
    Freq_Orig = Structure.freq;
    
    StandardDiv = GUI_Inputs.FWHM./(2*sqrt(2*log(2)));
    P_FlucCorr  = GUI_Inputs.P_FlucCorr/100; % turn percentage to number within 0~1
    
    % Add diagonal disorder
    for j = 1:Sample_Num
        Structure_Fluc(j) = Structure;
        
        Correlation_Dice = rand;

        if Correlation_Dice < P_FlucCorr
            Fluctuation = StandardDiv'.*(randn(1,1).*ones(Num_Modes,1));
        else 
            Fluctuation = StandardDiv'.*randn(Num_Modes,1); 
        end
        Structure_Fluc(j).freq = Freq_Orig + Fluctuation;
    end
    
    
    TSTART = zeros(GUI_Inputs.Sample_Num,1,'uint64');
    TIME   = zeros(GUI_Inputs.Sample_Num,1);
    
    parfor i = 1:GUI_Inputs.Sample_Num
        
        TSTART(i) = tic;
                
        [Tmp_SG,Tmp_Res] = TwoDSFG_Main(Structure_Fluc(i),GUI_Inputs);
        
        Rephasing(:,:,i)    = Tmp_SG.Rephasing   ;
        NonRephasing(:,:,i) = Tmp_SG.NonRephasing;
        SpecR1(:,:,i)       = Tmp_SG.SpecAccuR1  ;
        SpecR2(:,:,i)       = Tmp_SG.SpecAccuR2  ;
        SpecR3(:,:,i)       = Tmp_SG.SpecAccuR3  ;
        SpecNR1(:,:,i)      = Tmp_SG.SpecAccuNR1 ;
        SpecNR2(:,:,i)      = Tmp_SG.SpecAccuNR2 ;
        SpecNR3(:,:,i)      = Tmp_SG.SpecAccuNR3 ;   
        
        TIME(i) = toc(TSTART(i));
        disp(['Run ' num2str(i) ' finished within '  num2str(TIME(i)) '...'])
        
    end

        SpectraGrid.Rephasing    = sum(Rephasing,3)   ;
        SpectraGrid.NonRephasing = sum(NonRephasing,3);
        SpectraGrid.SpecAccuR1   = sum(SpecR1,3)      ;
        SpectraGrid.SpecAccuR2   = sum(SpecR2,3)      ;
        SpectraGrid.SpecAccuR3   = sum(SpecR3,3)      ;
        SpectraGrid.SpecAccuNR1  = sum(SpecNR1,3)     ;
        SpectraGrid.SpecAccuNR2  = sum(SpecNR2,3)     ;
        SpectraGrid.SpecAccuNR3  = sum(SpecNR3,3)     ;
        
        Response = Tmp_Res;
    
else
    [SpectraGrid,Response] = TwoDSFG_Main(Structure,GUI_Inputs);
end

Total_TIME = toc(TSTART_whole);
disp(['Total time: ' num2str(Total_TIME)])
