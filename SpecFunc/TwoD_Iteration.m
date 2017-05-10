function [SpectraGrid,Response] = TwoD_Iteration(h2DFunc,GUI_data,GUI_Inputs,hGUIs)
%% Read GUI
DynamicUpdate = GUI_Inputs.DynamicUpdate;
FreqRange     = GUI_Inputs.FreqRange;
Structure     = GUI_data.Structure;

%% Pre-calculation settings
% Create figure object for dynamics figure update
if DynamicUpdate
    hF = figure;
end

%% Calculate TwoD response
if eq(GUI_Inputs.Sampling,1)
    % Pre-allocate
    FreqRange = FreqRange(1):FreqRange(end)+100; % add 100 cm-1 range to prevent fluctuation out of range
    GUI_Inputs.FreqRange = FreqRange; % pass the extended Frequency Range to the TwoD main function
    GridSize  = FreqRange(end); 
    
    Rephasing    = sparse(GridSize,GridSize);
    NonRephasing = sparse(GridSize,GridSize);
    SpecAccuR1   = sparse(GridSize,GridSize);
    SpecAccuR2   = sparse(GridSize,GridSize);
    SpecAccuR3   = sparse(GridSize,GridSize);
    SpecAccuNR1  = sparse(GridSize,GridSize);
    SpecAccuNR2  = sparse(GridSize,GridSize);
    SpecAccuNR3  = sparse(GridSize,GridSize);

    Num_Modes = Structure.Nmodes;
    Freq_Orig = Structure.LocFreq;

    StandardDiv = GUI_Inputs.FWHM./(2*sqrt(2*log(2)));
    P_FlucCorr  = GUI_Inputs.P_FlucCorr/100; % turn percentage to number within 0~1

    TSTART = zeros(GUI_Inputs.Sample_Num,1,'uint64');
    TIME   = zeros(GUI_Inputs.Sample_Num,1);

    for i = 1:GUI_Inputs.Sample_Num

        DynamicUpdate = hGUIs.DynamicUpdate.Value; % directly access the GUI elment so can get the most recnt values
        UpdateStatus  = hGUIs.UpdateStatus.Value;
        if and(~eq(i,1), and(eq(DynamicUpdate,1),~eq(UpdateStatus,1)))
            break
        end

        TSTART(i) = tic;

        % Add diagonal disorder
        Correlation_Dice = rand;
        if Correlation_Dice < P_FlucCorr
            Fluctuation = StandardDiv'.*(randn(1,1).*ones(Num_Modes,1));
        else 
            Fluctuation = StandardDiv'.*randn(Num_Modes,1); 
        end
        Structure.LocFreq = Freq_Orig + Fluctuation;

        % run main function
        [Tmp_SG,Tmp_Res] = h2DFunc(Structure,GUI_Inputs);
        
        % Accumulate result
        try
            Rephasing    = Rephasing    + Tmp_SG.Rephasing   ;
            NonRephasing = NonRephasing + Tmp_SG.NonRephasing;
            SpecAccuR1   = SpecAccuR1   + Tmp_SG.SpecAccuR1  ;
            SpecAccuR2   = SpecAccuR2   + Tmp_SG.SpecAccuR2  ;
            SpecAccuR3   = SpecAccuR3   + Tmp_SG.SpecAccuR3  ;
            SpecAccuNR1  = SpecAccuNR1  + Tmp_SG.SpecAccuNR1 ;
            SpecAccuNR2  = SpecAccuNR2  + Tmp_SG.SpecAccuNR2 ;
            SpecAccuNR3  = SpecAccuNR3  + Tmp_SG.SpecAccuNR3 ;   
        catch
            disp(['Frequency fluctuation out of range: ', num2str(GridSize),', dop this run...'])
            continue
        end

        SpectraGrid.Rephasing    = Rephasing    ;
        SpectraGrid.NonRephasing = NonRephasing ;
        SpectraGrid.SpecAccuR1   = SpecAccuR1   ;
        SpectraGrid.SpecAccuR2   = SpecAccuR2   ;
        SpectraGrid.SpecAccuR3   = SpecAccuR3   ;
        SpectraGrid.SpecAccuNR1  = SpecAccuNR1  ;
        SpectraGrid.SpecAccuNR2  = SpecAccuNR2  ;
        SpectraGrid.SpecAccuNR3  = SpecAccuNR3  ;
        Response = Tmp_Res;

        % Dynamic update of figure % update every 10 run
        while and(~eq(DynamicUpdate,0),eq(mod(i,10),0))
            CVL = Conv2D(SpectraGrid,GUI_Inputs);
            CVL.FilesName = [Structure.FilesName,' ',num2str(i),'-th run...']; % pass filesname for figure title
            SpecType = Response.SpecType;
            Plot2D(hF,CVL,GUI_Inputs,SpecType);
            drawnow
            DynamicUpdate = 0;
        end
        
        TIME(i) = toc(TSTART(i));
        disp(['Run ' num2str(i) ' finished within '  num2str(TIME(i)) '...'])
    end
    Total_TIME = sum(TIME);
    disp(['Total time: ' num2str(Total_TIME)])

else
    [SpectraGrid,Response] = h2DFunc(Structure,GUI_Inputs);
end
