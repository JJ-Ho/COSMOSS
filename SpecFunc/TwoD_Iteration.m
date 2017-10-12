function [SpectraGrid,Response] = TwoD_Iteration(h2DFunc,app)
%% Read GUI
I = Parse_COSMOSS(app);
Sampling      = I.Sampling;
DynamicUpdate = I.DynamicUpdate;
Sample_Num    = I.Sample_Num;
FreqRange     = I.FreqRange;
Structure     = app.Structure;

%% Pre-calculation settings
% Create figure object for dynamics figure update
if DynamicUpdate && Sampling
    hF = figure;
    hAx = axes('Parent',hF);
end

%% Calculate TwoD response
if Sampling
    % Pre-allocate
    FreqRange = FreqRange(1):FreqRange(end)+100; % add 100 cm-1 range to prevent fluctuation out of range
    I.FreqRange = FreqRange; % pass the extended Frequency Range to the TwoD main function
    GridSize  = FreqRange(end); 
    
    R1   = sparse(GridSize,GridSize);
    R2   = sparse(GridSize,GridSize);
    R3   = sparse(GridSize,GridSize);
    NR1  = sparse(GridSize,GridSize);
    NR2  = sparse(GridSize,GridSize);
    NR3  = sparse(GridSize,GridSize);

    TSTART = zeros(Sample_Num,1,'uint64');
    TIME   = zeros(Sample_Num,1);

    for i = 1:Sample_Num
        TSTART(i) = tic;
        
        % Read GUI input directly from obj handle
        DynamicUpdate = app.CheckBox_DynamicFigUpdate.Value;
        UpdateStatus  = app.CheckBox_Continue.Value;
        if and(~eq(i,1), and(eq(DynamicUpdate,1),~eq(UpdateStatus,1)))
            break
        end

        % run main function
        [Tmp_SG,Tmp_Res] = h2DFunc(Structure,I);
        
        % Accumulate result
        try
            R1   = R1   + Tmp_SG.R1  ;
            R2   = R2   + Tmp_SG.R2  ;
            R3   = R3   + Tmp_SG.R3  ;
            NR1  = NR1  + Tmp_SG.NR1 ;
            NR2  = NR2  + Tmp_SG.NR2 ;
            NR3  = NR3  + Tmp_SG.NR3 ;   
        catch
            disp(['Frequency fluctuation out of range: ', num2str(GridSize),', dop this run...'])
            continue
        end

        SpectraGrid.R1   = R1   ;
        SpectraGrid.R2   = R2   ;
        SpectraGrid.R3   = R3   ;
        SpectraGrid.NR1  = NR1  ;
        SpectraGrid.NR2  = NR2  ;
        SpectraGrid.NR3  = NR3  ;
        Response = Tmp_Res;

        % Dynamic update of figure % update every 10 run
        while and(~eq(DynamicUpdate,0),eq(mod(i,10),0))
            CVL = Conv2D(SpectraGrid,I);
            CVL.FilesName = [Structure.FilesName,' ',num2str(i),'-th run...']; % pass filesname for figure title
            SpecType = Response.SpecType;
            Plot2D(hAx,CVL,I,SpecType);
            drawnow
            DynamicUpdate = 0;
        end
        
        TIME(i) = toc(TSTART(i));
        disp(['Run ' num2str(i) ' finished within '  num2str(TIME(i)) '...'])
    end
    Total_TIME = sum(TIME);
    disp(['Total time: ' num2str(Total_TIME)])

else
    [SpectraGrid,Response] = h2DFunc(Structure,I);
end
