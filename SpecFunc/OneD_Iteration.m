function OneD = OneD_Iteration(h1DFunc,app)
I = app.Parse_GUI;
S = app.Structure;

Sampling   = I.Sampling;
FreqRange  = I.FreqRange_1D;
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




if eq(Sampling,1)   
    % Pre-allocate
    GridSize   = length(FreqRange);
    Response1D = zeros(GridSize,1);
    
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
        
%         % deal with MD snapshots(Need to modify the reading mechanism later)
%         if isfield(S,'NStucture') && gt(S.NStucture,1)
%             S.LocMu = squeeze(S.LocMu(:,:,i));
%             S.LocCenter = S.LocCenter(:,:,i);
%         end
        
        % Calculate 1D respsonse
        OneD = h1DFunc(S,I);
        
        % recursive accumulation of signal
        Response1D = Response1D + OneD.Response1D; % note freq is binned and sported, so direct addition work
        OneD.Response1D = Response1D;
        
        TIME(i) = toc(TSTART(i));
        disp(['Run ' num2str(i) ' finished within '  num2str(TIME(i)) '...'])
        
        while ~eq(DynamicUpdate,0)
            OneD.FilesName = [FilesName,' ',num2str(i),'-th run...']; % pass filesname for figure title
            cla(hAx)
            Plot1D(hAx,OneD,I);
            drawnow
            DynamicUpdate = 0;
        end
    end
    
        Plot1D(hAx,OneD,I);
        Total_TIME = sum(TIME);
        disp(['Total time: ' num2str(Total_TIME)])
        
else
    OneD = h1DFunc(S,I);
    Plot1D(hAx,OneD,I);
end