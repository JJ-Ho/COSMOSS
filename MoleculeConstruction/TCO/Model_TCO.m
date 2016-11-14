%^ GUI Setup ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function varargout = Model_TCO(varargin)
% Model_TCO MATLAB code for Model_TCO.fig
%      Model_TCO, by itself, creates a new Model_TCO or raises the existing
%      singleton*.
%
%      H = Model_TCO returns the handle to a new Model_TCO or the handle to
%      the existing singleton*.
%
%      Model_TCO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Model_TCO.M with the given input arguments.
%
%      Model_TCO('Property','Value',...) creates a new Model_TCO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Model_TCO_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Model_TCO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Model_TCO

% Last Modified by GUIDE v2.5 01-Oct-2014 16:16:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Model_TCO_OpeningFcn, ...
                   'gui_OutputFcn',  @Model_TCO_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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

function Model_TCO_OpeningFcn(hModel_TCO, eventdata, GUI_data, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Model_TCO (see VARARGIN)
if nargin > 3
    switch varargin{1}
        case 'COSMOSS'
            hCOSMOSS = varargin{2};
            Data_COSMOSS = guidata(hCOSMOSS);
            
            %PRE ASSIGN VALUES TO SUBSTITUTE MAIN GUI VALUES
            hGUIs_COSMOSS  = Data_COSMOSS.hGUIs;
            set(hGUIs_COSMOSS.X_Min  ,'String','1650')
            set(hGUIs_COSMOSS.X_Max  ,'String','1750')
            
            GUI_data.hCOSMOSS = hCOSMOSS;
            
            disp('Running Model_TCO directly from COSMOSS...')
        case 'Comb2'
            hModel_Comb2 = varargin{2};
            Comb2_Order  = varargin{3};
            
            GUI_data.hModel_Comb2 = hModel_Comb2;
            GUI_data.Comb2_Order  = Comb2_Order;
            
            % Add comb2 order # to GUI title, if necessary
            TitleStr = hModel_TCO.Name;
            if ~strcmp(TitleStr(1),'#')
                hModel_TCO.Name = ['#',int2str(Comb2_Order),', ',TitleStr];
            end
            
            disp('Running Model_TCO as a sub GUI of Comb2...')
    end
else
    disp('Running Model_TCO in stand alone mode...')
end

% Call create Interface to create GUI elements
hGUIs = GUI_TCO(hModel_TCO);

% Prep necessary data to be saved in GUI_data
GUI_data.hModel_TCO = hModel_TCO;
GUI_data.hGUIs      = hGUIs;

% Update handles structure
guidata(hModel_TCO, GUI_data);

% Generate the default structure
UpdateStructure(hModel_TCO, eventdata, GUI_data)

function varargout = Model_TCO_OutputFcn(hModel_TCO, eventdata, GUI_data) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = hModel_TCO;

function Export_Handle_Callback(hModel_TCO, eventdata, GUI_data)
% export handles back to work space
assignin('base', 'Data_TCO', GUI_data)
disp('Updated GUI Data_TCO exported!')
%^ GUI Setup ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




function UpdateStructure(hObject, eventdata, GUI_data)
%% Read GUI variables
hGUIs = GUI_data.hGUIs;
GUI_Inputs = ParseGUI_TCO(hGUIs);

%% Construct molecule
Structure = GetAcid(GUI_Inputs);

% Export into Structure so it can be passsed around different GUIs
Structure.StructModel = 1;                

%% Export result to Main guidata
GUI_data.Structure = Structure;

% include FieldName of GUI Inputs
[~,~,~,hGUIParser] = StructureModel(Structure.StructModel);
[~,GUI_FieldName] = hGUIParser(hGUIs);
GUI_data.GUI_FieldName = GUI_FieldName;

guidata(hObject,GUI_data)

% update to other GUIs
Export2GUIs(GUI_data)

disp('Structure file generated!')

function hF = PlotMolecule(hObject, eventdata, GUI_data)
%% Read GUI variables
hGUIs = GUI_data.hGUIs;
GUI_Inputs = ParseGUI_TCO(hGUIs);

hF = PlotXYZfiles_Acid(GUI_data.Structure,GUI_Inputs);

function PlotModes(hObject, eventdata, GUI_data)
Plot_Modes(GUI_data.hModel);

