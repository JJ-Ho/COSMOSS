%% Setup canning parameters
distanceArray = 3:0.1:7;

GUI_Inputs = struct;
GUI_Inputs.Theta_D1 = -30;
GUI_Inputs.Theta_D2 = 30;
GUI_Inputs.FreqRange = 1650:1800;

% figure setting
hF  = figure;
hF.Position = [100,100,800,300];
hAx_FTIR     = subplot(1,2,1,'Parent',hF);
hAx_Molecule = subplot(1,2,2,'Parent',hF);

%% Strueture generation, simulation, and make figure
for i = 1:length(distanceArray)
    
    % clear the axises before start
    cla(hAx_FTIR)
    cla(hAx_Molecule)
    
    % Prepare initial structure
    GUI_Inputs.Displacement = [distanceArray(i),0,0];
    SD_TCO = ConstructTCO(GUI_Inputs);

    % Run the Spectral simulation
    OneD = FTIR_Main(SD_TCO,GUI_Inputs);

    % make figure
    Plot1D(hAx_FTIR,OneD,GUI_Inputs);
    SD_Draw(SD_TCO,hAx_Molecule);
    hAx_Molecule.XLim = [-6,6];
    
    drawnow
    Frame_all(i) = getframe(hF);
    
end 

%% Save GIF
saveGIF  = 0;
BaseFileName = 'TCO_X3_X10';
SavePath = '~/Desktop';
SaveName = [SavePath,'/',BaseFileName,'.gif'];
DT = 0.1;

if saveGIF
    for j = 1:length(distanceArray)

        frame = Frame_all(j);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);

        if j == 1
          imwrite(imind,cm,SaveName,'gif','Loopcount',inf,'DelayTime',DT);
        else
          imwrite(imind,cm,SaveName,'gif','WriteMode','append','DelayTime',DT);
        end

    end
    disp([BaseFileName, 'GIF saved...'])
end

%% Play movie frames
hF  = figure;
hF.Position = [100,100,800,300];
hAx_FTIR     = subplot(1,2,1,'Parent',hF);
hAx_Molecule = subplot(1,2,2,'Parent',hF);
movie(hF,Frame_all,3)