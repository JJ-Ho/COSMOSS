function output = iteration(app,specFunc,simulationType)
%% Deal with common parameters
I = app.Parse_GUI;
S = app.Structure;

Sampling   = I.Sampling;
Sample_Num = I.Sample_Num;
existFig   = I.existFig;
hFig       = I.hFig;
FilesName  = S.FilesName;

%% deal with figure handles
if existFig
    hAx = findobj(evalin('base',hFig),'Type','Axes');
else
    hF  = figure;
    hAx = axes('Parent',hF);
end

%% iteration or single run
if eq(Sampling,1)   
    
    % prepare axes for histogram
    hF_histo = figure;
    hAx_hiso = axes('Parent',hF_histo);
    
    % Pre-allocate data ouputs
    Freq = nan(S.Nmodes,Sample_Num); % this is for analysing local mode frequency distribution
    switch simulationType
        case 'oneD'
            GridSize   = length(I.FreqRange_1D);
            Response1D = zeros(GridSize,1);
        case 'twoD'
            % run the 2D simulation once
            [~,Res1]   = specFunc(S,I);
            GridSize   = Res1.SparseMax+ 100; 
            I.F_Max_2D = GridSize; % force all the iterations to take this value as it's grid size
            R1   = sparse(GridSize,GridSize);
            R2   = sparse(GridSize,GridSize);
            R3   = sparse(GridSize,GridSize);
            NR1  = sparse(GridSize,GridSize);
            NR2  = sparse(GridSize,GridSize);
            NR3  = sparse(GridSize,GridSize);
    end
    
    % prepare the timer
    TSTART = zeros(Sample_Num,1,'uint64');
    TIME   = zeros(Sample_Num,1);
    
    % iteration loops
    for i = 1:Sample_Num
        TSTART(i) = tic;
                
        % Calculate and accumulate 1D/2D respsonses
        switch simulationType
            case 'oneD' 
                Tmp_1D = specFunc(S,I);
                Response1D = Response1D + Tmp_1D.Response1D; % recursive accumulation of signal, note freq is binned and sorted, so direct addition work
                
                Tmp_1D.Response1D = Response1D;
                Freq(:,i) = Tmp_1D.H.dLocFreq;
            case 'twoD'
                [Tmp_SG,Tmp_Res] = specFunc(S,I);
                try % Accumulate result
                    R1   = R1   + Tmp_SG.R1  ;
                    R2   = R2   + Tmp_SG.R2  ;
                    R3   = R3   + Tmp_SG.R3  ;
                    NR1  = NR1  + Tmp_SG.NR1 ;
                    NR2  = NR2  + Tmp_SG.NR2 ;
                    NR3  = NR3  + Tmp_SG.NR3 ;   
                catch
                    disp([num2str(Tmp_Res.SparseMax),' is larger than the set grid size: ',num2str(GridSize),' drop this run...'])
                    continue
                end

                SpectraGrid.R1   = R1   ;
                SpectraGrid.R2   = R2   ;
                SpectraGrid.R3   = R3   ;
                SpectraGrid.NR1  = NR1  ;
                SpectraGrid.NR2  = NR2  ;
                SpectraGrid.NR3  = NR3  ;
                Response = Tmp_Res;
                Freq(:,i) = Tmp_Res.H.dLocFreq;
        end
        
        
        % Dynamic update of figure and export output
        DynamicUpdate = app.CheckBox_DynamicFigUpdate.Value;
        UpdateStatus  = app.CheckBox_Continue.Value;
        if and(~eq(i,1), and(eq(DynamicUpdate,1),~eq(UpdateStatus,1)))
            break
        end
        cla(hAx)
        output.Freq = Freq;
        
        while eq(DynamicUpdate,1)
            % draw histogram of sampling frequency
            histogram(hAx_hiso,Freq(:))
            Title_String = [FilesName,' ',num2str(i),'-th run...'];
            title(hAx_hiso,Title_String,'FontSize',16);
            
            switch simulationType            
                case 'oneD' 
                    Tmp_1D.FilesName = [FilesName,' local mode sampling ',num2str(i),'-th run...'];
                    Plot1D(hAx,Tmp_1D,I);
                    drawnow
                case 'twoD' % Dynamic update of figure, update every 10 run 
                    if eq(mod(i,3),0)
                        CVL = Conv2D(SpectraGrid,I);
                        CVL.FilesName = [FilesName,' ',num2str(i),'-th run...'];
                        Plot2D(hAx,CVL,I,Response.SpecType);
                        drawnow
                    end
            end
            DynamicUpdate = 0;
        end

        TIME(i) = toc(TSTART(i));
        disp(['Run ' num2str(i) ' finished within '  num2str(TIME(i)) '...'])
    end
    
    % draw the final spectra and export output
    switch simulationType
        case 'oneD'
            Plot1D(hAx,Tmp_1D,I);
            output.Response = Tmp_1D;
        case 'twoD'
            CVL = Conv2D(SpectraGrid,I);
            CVL.FilesName = [FilesName,' ',num2str(i),'-th run...'];
            Plot2D(hAx,CVL,I,Response.SpecType);
            output.SpectraGrid = SpectraGrid;
            output.Response    = Response;
            output.CVL         = CVL;
    end
    
    Total_TIME = sum(TIME);
    disp(['Total time: ' num2str(Total_TIME)])
        
else
    % if not sampling, run single simulation
    output.Freq = S.LocFreq;
    switch simulationType
        case 'oneD'
            Response = specFunc(S,I);
            Plot1D(hAx,Response,I);
            
            output.Response = Response;
            
        case 'twoD'
            [SpectraGrid,Response] = specFunc(S,I);
            CVL = Conv2D(SpectraGrid,I);
            CVL.FilesName = FilesName;
            Plot2D(hAx,CVL,I,Response.SpecType);
            
            output.SpectraGrid = SpectraGrid;
            output.Response    = Response;
            output.CVL         = CVL;
    end
end