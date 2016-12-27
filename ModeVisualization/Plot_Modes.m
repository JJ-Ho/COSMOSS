function varargout = Plot_Modes(varargin)
% PLOT_MODES MATLAB code for Plot_Modes.fig
%      PLOT_MODES, by itself, creates a new PLOT_MODES or raises the existing
%      singleton*.
%
%      H = PLOT_MODES returns the handle to a new PLOT_MODES or the handle to
%      the existing singleton*.
%
%      PLOT_MODES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOT_MODES.M with the given input arguments.
%
%      PLOT_MODES('Property','Value',...) creates a new PLOT_MODES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Plot_Modes_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Plot_Modes_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Plot_Modes

% Last Modified by GUIDE v2.5 01-Nov-2016 23:52:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Plot_Modes_OpeningFcn, ...
                   'gui_OutputFcn',  @Plot_Modes_OutputFcn, ...
                   'gui_LayoutFcn',  @GUI_Base_Plot_Modes, ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function hPlot_Modes = GUI_Base_Plot_Modes(Singleton)
% this function create base figure. To utilize the original GUI building
% mechanism of Matlab and avoid Matlab default way to export GUI element 
% handles, this function only create the base layer with no GUI elements.
% Aftter calling "gui_mainfcn" we will call "GUI_COSMOSS" to build GUI
% elements.

% Create base figure
hPlot_Modes = figure;

hPlot_Modes.Units            = 'Pixels';
hPlot_Modes.Position         = [2 53 610 600];
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

function Plot_Modes_OpeningFcn(hPlot_Modes, eventdata, GUI_data, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Plot_Modes (see VARARGIN)
if nargin > 3    
    % Get Structural modeling GUI's handles
    if ishandle(varargin{1}) 
       hModel = varargin{1};
    end
else
    disp('Plot_Modes cannot run without COSMOSS... ')  
    return
end

% Call createInterface to create GUI elements
hGUIs = GUI_Plot_Modes(hPlot_Modes);

% Change Names on Plot_Modes GUI to identify which Structural model is
% plotting
Model_Name = hModel.Name;
Title_Str  = hGUIs.hPlot_Modes.Name;
hGUIs.hPlot_Modes.Name = [Title_Str, ': ', Model_Name];

% Prep necessary data to be saved in GUI_data
GUI_data.hPlot_Modes = hPlot_Modes;
GUI_data.hGUIs       = hGUIs; % export GUI handles to handles
GUI_data.hModel      = hModel;
GUI_data.hCOSMOSS    = hModel.UserData; % get hCOSMOSS from the UserData of model fig file
guidata(hPlot_Modes, GUI_data);

% update exciton info in Plot_Exciton
Update_Modes(hPlot_Modes, eventdata, GUI_data)

function varargout = Plot_Modes_OutputFcn(hPlot_Modes, eventdata, GUI_data) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = hPlot_Modes;

function Export_Handle_Callback(hPlot_Modes, eventdata, GUI_data)
% export handles back to work space
assignin('base', 'Data_Plot_Modes', GUI_data)
disp('Updated GUI Data_Plot_Modes exported!')
%^ GUI Setup ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




function Update_Modes(hObject, eventdata, GUI_data)
%% Gather GUI inputs
GUI_Data_hMain  = guidata(GUI_data.hCOSMOSS);
COSMOSS_Inputs  = ParseGUI_Main(GUI_Data_hMain.hGUIs);

GUI_Data_hModel = guidata(GUI_data.hModel);
Structure       = GUI_Data_hModel.Structure;

GUI_Inputs      = ParseGUI_Modes(GUI_data.hGUIs);
SpecType        = GUI_Inputs.SpecType;

%% Update GUI
% update table contents
switch SpecType
    case 1 %'FTIR'
        T = Update_Modes_Table('FTIR',Structure,COSMOSS_Inputs);
        GUI_data.hGUIs.Raman_Plot.Value = 0;
        GUI_data.hGUIs.Raman_Plot.Enable = 'off';
    case 2 %'SFG'
        T = Update_Modes_Table('SFG',Structure,COSMOSS_Inputs);
        GUI_data.hGUIs.Raman_Plot.Value = 1;
        GUI_data.hGUIs.Raman_Plot.Enable = 'on';
    case 3 %'2DIR'
        T = Update_Modes_Table('TwoDIR',Structure,COSMOSS_Inputs);
        GUI_data.hGUIs.Raman_Plot.Value = 0;
        GUI_data.hGUIs.Raman_Plot.Enable = 'off';
    case 4 %'2DSFG'
        T = Update_Modes_Table('TwoDSFG',Structure,COSMOSS_Inputs);
        GUI_data.hGUIs.Raman_Plot.Value = 1;
        GUI_data.hGUIs.Raman_Plot.Enable = 'on';
end

GUI_data.hGUIs.ModeList.ColumnName   = T.ModeList.Name;
GUI_data.hGUIs.ModeList.ColumnFormat = T.ModeList.Format;
GUI_data.hGUIs.ModeList.ColumnWidth  = T.ModeList.Width;
GUI_data.hGUIs.ModeList.Data         = T.ModeList.Data;
GUI_data.hGUIs.SortInd.String        = T.ModeList.Name;

% call sorting to sort table with the same GUI setting
uitable_Sort(hObject, eventdata, GUI_data)

%% Update handles structure
GUI_data.SpecData  = T.SpecData;
GUI_data.Structure = Structure;

guidata(hObject, GUI_data);

function Update_TDV_Raman(hObject, eventdata, GUI_data)
%% Gather necessary inputs
GUI_Inputs = ParseGUI_Modes(GUI_data.hGUIs);
Structure  = GUI_data.Structure;
SpecData   = GUI_data.SpecData;
hModel     = GUI_data.hModel;

%% Draw molecule by calling the PlotMolecule function in each model
[hFunc_Model,~,~] = StructureModel(Structure.StructModel);
hF = hFunc_Model('PlotMolecule',hModel,eventdata,guidata(hModel));

%% Call Update Figure function
Fig_Output = Update_Modes_Figure(hF, GUI_Inputs, Structure, SpecData);

%% update handles
GUI_data.hF           = hF;
GUI_data.Mu_Alpha_Ind = Fig_Output.Mu_Alpha_Ind;
guidata(hObject,GUI_data)

function Update_Response(hObject, eventdata, GUI_data)
%% Gather necessary inputs
GUI_Inputs      = ParseGUI_Modes(GUI_data.hGUIs);
Structure       = GUI_data.Structure;
SpecData        = GUI_data.SpecData;
hModel          = GUI_data.hModel;
GUI_Data_Main   = guidata(GUI_data.hCOSMOSS);
GUI_Inputs_Main = ParseGUI_Main(GUI_Data_Main.hGUIs);

%% Call Update RespF function
Response = Fig_Response(hModel, GUI_Inputs, Structure, SpecData, GUI_Inputs_Main);

%% update handles
GUI_data.Response = Response;
guidata(hObject,GUI_data)

function uitable_CellSelection(hObject, eventdata, GUI_data)
hGUIs = GUI_data.hGUIs;

TableData   = hGUIs.ModeList.Data;
CurrentCell = eventdata.Indices;

if isempty(CurrentCell)
    disp('No mode selected, please select at least one mode...')
    return
else
    CurrentRowInd = CurrentCell(:,1)';
    Mode_Ind_Str_Full  = num2str(cell2mat(TableData(CurrentRowInd   ,1)'));
    Mode_Ind_Str_1st   = num2str(cell2mat(TableData(CurrentRowInd(1),1)'));

    % Update the Mode index on GUI
    hGUIs.Mu_Alpha_Ind.String = Mode_Ind_Str_Full;
    hGUIs.EigVec_Ind.String   = Mode_Ind_Str_1st; % only take the first index of slection, since mixxing coefficient only take one mode
end
%% update handles
GUI_data.Mode_Ind_Str = Mode_Ind_Str_Full;
guidata(hObject,GUI_data)

function uitable_Sort(hObject, eventdata, GUI_data)
TableData  = GUI_data.hGUIs.ModeList.Data;
SortColumn = GUI_data.hGUIs.SortInd.Value;

SC = abs(cell2mat(TableData(:,SortColumn))); % only sort the absolute value

if eq(SortColumn,1) 
    % for sorting index
    [~,SortInd] = sort(SC,'ascend');
else
    [~,SortInd] = sort(SC,'descend');
end

SortedData = TableData(SortInd,:);

% Update Table on GUI
set(GUI_data.hGUIs.ModeList,'Data', SortedData);


