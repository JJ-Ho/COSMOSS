function OneD = OneD_Iteration(h1DFunc,Structure,GUI_Inputs,hGUIs)
%% Create figure object for dynamics figure update
DynamicUpdate = GUI_Inputs.DynamicUpdate;
Sampling      = GUI_Inputs.Sampling;

if DynamicUpdate && Sampling
    hF = figure;
    hAx = axes('Parent',hF);
end

%% Calculation of 1D response 
if eq(GUI_Inputs.Sampling,1)
    % Pre-allocate
    GridSize   = length(GUI_Inputs.FreqRange);
    Response1D = zeros(GridSize,1);
    
    TSTART = zeros(GUI_Inputs.Sample_Num,1,'uint64');
    TIME   = zeros(GUI_Inputs.Sample_Num,1);
    
    for i = 1:GUI_Inputs.Sample_Num
        TSTART(i) = tic;
        
        % Read GUI input directly from obj handle
        DynamicUpdate = hGUIs.DynamicUpdate.Value;
        UpdateStatus  = hGUIs.UpdateStatus.Value;
        if and(~eq(i,1), and(eq(DynamicUpdate,1),~eq(UpdateStatus,1)))
            break
        end
        
%         % deal with MD snapshots(Need to modify the reading mechanism later)
%         if isfield(S,'NStucture') && gt(S.NStucture,1)
%             S.LocMu = squeeze(S.LocMu(:,:,i));
%             S.LocCenter = S.LocCenter(:,:,i);
%         end
        
        % Calculate 1D respsonse
        OneD = h1DFunc(Structure,GUI_Inputs);
        
        % recursive accumulation of signal
        Response1D = Response1D + OneD.Response1D; % note freq is binned and sported, so direct addition work
        OneD.Response1D = Response1D;
        
        TIME(i) = toc(TSTART(i));
        disp(['Run ' num2str(i) ' finished within '  num2str(TIME(i)) '...'])
        
        while ~eq(DynamicUpdate,0)
            OneD.FilesName = [Structure.FilesName,' ',num2str(i),'-th run...']; % pass filesname for figure title
            Plot1D(hAx,OneD,GUI_Inputs);
            drawnow
            DynamicUpdate = 0;
        end
    end
    
        Total_TIME = sum(TIME);
        disp(['Total time: ' num2str(Total_TIME)])
        
else
    OneD = h1DFunc(Structure,GUI_Inputs);
end