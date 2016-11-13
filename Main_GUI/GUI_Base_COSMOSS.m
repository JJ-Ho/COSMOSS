function hCOSMOSS = GUI_COSMOSS_Base(Singleton)
% this function create base figure. To utilize the original GUI building
% mechanism of Matlab and avoid Matlab default way to export GUI element 
% handles, this function only create the base layer with no GUI elements.
% Aftter calling "gui_mainfcn" we will call "GUI_COSMOSS" to build GUI
% elements.

%% Create base figure
hCOSMOSS = figure;

hCOSMOSS.Position         = [103.8333 61.667 600.0000 600.0000];
hCOSMOSS.Name             = 'COSMOSS';
hCOSMOSS.ToolBar          = 'none';
hCOSMOSS.MenuBar          = 'none';
hCOSMOSS.NumberTitle      = 'off';
hCOSMOSS.IntegerHandle    = 'off';
hCOSMOSS.Tag              = 'Main'; % tag to distinguish type of GUI
hCOSMOSS.HandleVisibility = 'Callback';

gui_Options.syscolorfig = 1;
setappdata(hCOSMOSS,'GUIDEOptions',gui_Options);

disp('Creating COSMOSS GUI Using GUI Layout Toolbox!')