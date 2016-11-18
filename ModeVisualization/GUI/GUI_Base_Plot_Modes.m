function hPlot_Modes = GUI_Base_Plot_Modes(Singleton)
% this function create base figure. To utilize the original GUI building
% mechanism of Matlab and avoid Matlab default way to export GUI element 
% handles, this function only create the base layer with no GUI elements.
% Aftter calling "gui_mainfcn" we will call "GUI_COSMOSS" to build GUI
% elements.

%% Create base figure
hPlot_Modes = figure;

hPlot_Modes.Units            = 'Pixels';
hPlot_Modes.Position         = [2 53 752 537];
hPlot_Modes.Name             = 'Plot_Modes';
hPlot_Modes.ToolBar          = 'none';
hPlot_Modes.MenuBar          = 'none';
hPlot_Modes.NumberTitle      = 'off';
hPlot_Modes.IntegerHandle    = 'off';
hPlot_Modes.Tag              = 'hPlot_Modes'; % tag to distinguish type of GUI
hPlot_Modes.HandleVisibility = 'Callback';

gui_Options.syscolorfig = 1;
setappdata(hPlot_Modes,'GUIDEOptions',gui_Options);

disp('Creating Plot_Modes GUI Using GUI Layout Toolbox!')
disp('...')