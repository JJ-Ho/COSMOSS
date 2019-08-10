function OverlayExp(hF_Sim)
% Load 2D Experiment 
Default_Path = '~/Desktop';

[FilesName,PathName,~] = uigetfile({'*.fig','2D SFG fig'; ...
                                    '*,*','All Files'},...
                                    'Select inputs',Default_Path);
disp([FilesName,' loaded...'])
hF_Exp = openfig([PathName,FilesName]);
hContour_Exp = findobj(hF_Exp,'Type','Contour');
hContour_Exp = hContour_Exp(1);

hAx = findobj(hF_Sim,'Type','Axes');
copyobj(hContour_Exp,hAx)
close(hF_Exp)