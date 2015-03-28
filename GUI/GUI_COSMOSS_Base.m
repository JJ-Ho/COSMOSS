function hMainFig = GUI_COSMOSS_Base(Singleton)
% this function create base figure. To utilize the original GUI building
% mechanism of Matlab and avoid Matlab default way to export GUI element 
% handles, this function only create the base layer with no GUI elements.
% Aftter calling "gui_mainfcn" we will call "GUI_COSMOSS" to build GUI
% elements.

%% Create base figure
hMainFig = figure;
hMainFig.Position = [103.8333 61.667 600.0000 600.0000];
hMainFig.Name = 'COSMOSS';
hMainFig.ToolBar = 'none';
hMainFig.MenuBar = 'none';
hMainFig.NumberTitle = 'off';
hMainFig.IntegerHandle = 'off';
hMainFig.Tag = 'Main';
gui_Options.syscolorfig = 1;
setappdata(hMainFig,'GUIDEOptions',gui_Options);
hMainFig.HandleVisibility= 'Callback';

disp('Using GUI Layout Toolbox!')