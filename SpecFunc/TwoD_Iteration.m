function [SpectraGrid,Response,CVL] = TwoD_Iteration(h2DFunc,app)
%% Read GUI
I = app.Parse_GUI;
S = app.Structure;

Sampling   = I.Sampling;
Sample_Num = I.Sample_Num;
existFig   = I.existFig;
hFig       = I.hFig;
FilesName  = S.FilesName;

if existFig
    hAx = findobj(evalin('base',hFig),'Type','Axes');
else
    hF  = figure;
    hAx = axes('Parent',hF);
end

%% Calculate TwoD response
if Sampling
    % run the 2D simulation once
    [SG1,~] = h2DFunc(S,I);
    
    % Pre-allocate
    GridSize   = size(SG1.R1,1) + 100; 
    I.F_Max_2D = GridSize; % force all the iterations to take this value as it's grid size
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
        [Tmp_SG,Tmp_Res] = h2DFunc(S,I);
        
        % Accumulate result
        try
            R1   = R1   + Tmp_SG.R1  ;
            R2   = R2   + Tmp_SG.R2  ;
            R3   = R3   + Tmp_SG.R3  ;
            NR1  = NR1  + Tmp_SG.NR1 ;
            NR2  = NR2  + Tmp_SG.NR2 ;
            NR3  = NR3  + Tmp_SG.NR3 ;   
        catch
            disp([num2str(Tmp_Res.SparseMax),' is larger than the set grid size: ',num2str(GridSize),' dop this run...'])
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
            CVL.FilesName = [FilesName,' ',num2str(i),'-th run...']; % pass filesname for figure title

            cla(hAx)
            Plot2D(hAx,CVL,I,Response.SpecType);
            drawnow
            DynamicUpdate = 0;
        end
        
        TIME(i) = toc(TSTART(i));
        disp(['Run ' num2str(i) ' finished within '  num2str(TIME(i)) '...'])
    end
    Total_TIME = sum(TIME);
    disp(['Total time: ' num2str(Total_TIME)])

else
    [SpectraGrid,Response] = h2DFunc(S,I);
    CVL = Conv2D(SpectraGrid,I);
    CVL.FilesName = FilesName;
    Plot2D(hAx,CVL,I,Response.SpecType);
end
