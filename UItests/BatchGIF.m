%% Inputs
hCOSMOSS = Data_COSMOSS.hCOSMOSS;
hModel   = Data_COSMOSS.hModel;

ScanP = 0:10:350;

%- Figure options ---------------------------------------------------------
hF = figure;
hF.Position = [100,100,800,800];
hAx_Mode     = subplot(2,2,1);
hAx_1D       = subplot(2,2,2);
hAx_Molecule = subplot(2,2,3);
hAx_2D       = subplot(2,2,4);


YLim_hAx_1D_FTIR = [-1,1];
YLim_hAx_1D_SFG  = [-1,1];
% Fig_Save = 0;
saveGIF  = 0;
BaseFileName = 'APB_R5S3_';
% PathName = pwd;
SavePath = '~/Desktop';

%% Get handles of all GUI elements
hGUIs_COSMOSS = guihandles(hCOSMOSS);
hGUIs_Model   = guihandles(hModel);

%% Update structure
for i = 1:length(ScanP)
    cla(hAx_Mode)
    cla(hAx_1D)
    cla(hAx_Molecule)
    cla(hAx_2D)
    
    eventdata.Source = 'External';
    eventdata.hF  = hF;
    
    % Structure
    hObject = hGUIs_Model.Theta_D;
    hObject.String = num2str(ScanP(i));
    eventdata.hAx = hAx_Molecule;
    Model_Betasheet_AmideI('UpdateStructure',hObject,eventdata,guidata(hObject))
    Data = guidata(hObject);
    Data.Structure.Draw(hAx_Molecule);
    
    % FTIR/SFG
    hObject = hGUIs_COSMOSS.hMain;
    eventdata.hAx = hAx_1D;
    yyaxis(hAx_1D,'left')
    hAx_1D.YLim = YLim_hAx_1D_FTIR;
    COSMOSS('FTIR_Callback',hObject,eventdata,guidata(hObject))
    yyaxis(hAx_1D,'right')
    hAx_1D.YLim = YLim_hAx_1D_SFG;
    COSMOSS('SFG_Callback',hObject,eventdata,guidata(hObject))
    
    % 2DSFG
    hObject = hGUIs_COSMOSS.hMain;
    eventdata.hAx = hAx_2D;
    COSMOSS('TwoDSFG_Callback',hObject,eventdata,guidata(hObject))
    
    drawnow
    Frame_all(i) = getframe(hF);
end
close all

%% Play movie and save GIF
hF = figure;
hF.Position = [100,100,800,800];
hAx_Mode     = subplot(2,2,1);
hAx_1D       = subplot(2,2,2);
hAx_Molecule = subplot(2,2,3);
hAx_2D       = subplot(2,2,4);
movie(hF,Frame_all,10)

%% 
if saveGIF
    %% Save Frame
    SaveName = [SavePath,'/',BaseFileName];
    FrameName = [SaveName,'Frame_0-360'];
    save(FrameName,'Frame_all') 

    %% save gif
    load(FrameName)
    filename = [SaveName,'0-360.gif'];
    DT = 0.5;

    for j = 1:length(ScanP)

        frame = Frame_all(j);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);

        if j == 1;
          imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',DT);
        else
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',DT);
        end

    end
    disp('GIF saved...')
end
 No newline at end of file
