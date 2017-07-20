function varargout = Model_Betasheet_AmideI(varargin)
% MODEL_BETASHEET_AMIDEI MATLAB code for Model_Betasheet_AmideI.fig
%      MODEL_BETASHEET_AMIDEI, by itself, creates a new MODEL_BETASHEET_AMIDEI or raises the existing
%      singleton*.
%
%      H = MODEL_BETASHEET_AMIDEI returns the handle to a new MODEL_BETASHEET_AMIDEI or the handle to
%      the existing singleton*.
%
%      MODEL_BETASHEET_AMIDEI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODEL_BETASHEET_AMIDEI.M with the given input arguments.
%
%      MODEL_BETASHEET_AMIDEI('Property','Value',...) creates a new MODEL_BETASHEET_AMIDEI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Model_Betasheet_AmideI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Model_Betasheet_AmideI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Model_Betasheet_AmideI

% Last Modified by GUIDE v2.5 20-Feb-2016 16:34:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Model_Betasheet_AmideI_OpeningFcn, ...
                   'gui_OutputFcn',  @Model_Betasheet_AmideI_OutputFcn, ...
                   'gui_LayoutFcn',  @GUI_Base_Betasheet_AmideI, ...
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

function hModel = GUI_Base_Betasheet_AmideI(Singleton)
% this function create base figure. To utilize the original GUI building
% mechanism of Matlab and avoid Matlab default way to export GUI element 
% handles, this function only create the base layer with no GUI elements.
% Aftter calling "gui_mainfcn" we will call "GUI_COSMOSS" to build GUI
% elements.

% Create base figure
hModel = figure;

hModel.Units            = 'Pixels';
hModel.Position         = [2 406 250 600];
hModel.Name             = 'Model_Betasheet_AmideI';
hModel.ToolBar          = 'none';
hModel.MenuBar          = 'none';
hModel.NumberTitle      = 'off';
hModel.IntegerHandle    = 'off';
hModel.Tag              = 'hModel'; % tag to distinguish type of GUI
hModel.HandleVisibility = 'Callback';

gui_Options.syscolorfig = 1;
setappdata(hModel,'GUIDEOptions',gui_Options);


function Model_Betasheet_AmideI_OpeningFcn(hModel_Betasheet_AmideI, eventdata, GUI_data, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Model_Betasheet_AmideI (see VARARGIN)
if nargin > 3
    switch varargin{1}
        case 'COSMOSS'
            hCOSMOSS = varargin{2};
            Data_COSMOSS = guidata(hCOSMOSS);
            
            %PRE ASSIGN VALUES TO SUBSTITUTE MAIN GUI VALUES
            hGUIs_COSMOSS  = Data_COSMOSS.hGUIs;
            set(hGUIs_COSMOSS.Beta_NN,'String','0.8')
            set(hGUIs_COSMOSS.F_Min  ,'String','1550')
            set(hGUIs_COSMOSS.F_Max  ,'String','1750')
            
            GUI_data.hCOSMOSS = hCOSMOSS;
            
            disp('Running Model_Betasheet_AmideI directly from COSMOSS...')
        case 'Comb2'
            hModel_Comb2 = varargin{2};
            Comb2_Order  = varargin{3};
            
            GUI_data.hModel_Comb2 = hModel_Comb2;
            GUI_data.Comb2_Order  = Comb2_Order;
            
            % Add comb2 order # to GUI title, if necessary
            TitleStr = hModel_Betasheet_AmideI.Name;
            if ~strcmp(TitleStr(1),'#')
                hModel_Betasheet_AmideI.Name = ['#',int2str(Comb2_Order),', ',TitleStr];
            end
            
            disp('Running Model_Betasheet_AmideI as a sub GUI of Comb2...')
    end
else
    disp('Running Model_Betasheet_AmideI in stand alone mode...')
end

% Call createInterface to create GUI elements
hGUIs = GUI_Betasheet_AmideI(hModel_Betasheet_AmideI);

% Prep necessary data to be saved in GUI_data
GUI_data.hModel_Betasheet_AmideI = hModel_Betasheet_AmideI;
GUI_data.hGUIs                   = hGUIs;

% Update handles structure
guidata(hModel_Betasheet_AmideI, GUI_data);

% Generate default structure from GUI inputs
UpdateStructure(hModel_Betasheet_AmideI, eventdata, GUI_data)

function varargout = Model_Betasheet_AmideI_OutputFcn(hModel_Betasheet_AmideI, eventdata, GUI_data) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = hModel_Betasheet_AmideI;

function Export_Handle_Callback(hObject, eventdata, GUI_data)
% export handles back to work space
assignin('base', 'Data_Betasheet_AmideI', GUI_data)
disp('Updated GUI Data_Betasheet_AmideI exported!')
%^ GUI Setup ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




function UpdateStructure(hObject, eventdata, GUI_data)
%% Read GUI variables
hGUIs  = GUI_data.hGUIs;
GUI_Inputs = ParseGUI_Betasheet(hGUIs);

%% Construct molecule
BB        = ConstuctBetaSheet(GUI_Inputs);
Structure = GetAmideI(...
                      BB.XYZ,...
                      BB.AtomName,...
                      BB.FilesName,...
                      GUI_Inputs);
                  
% Export into Structure so it can be passsed around different GUIs
Structure.StructModel = 4;
%% Export extra info into Structure
Structure.Extra.N_Residue         = BB.N_Residue;
Structure.Extra.N_Strand          = BB.N_Strand;
Structure.Extra.N_Mode_per_Starnd = BB.N_Residue-1;

% C terminus Index
Structure.Extra.Ind_H = BB.Ind_H;
Structure.Extra.Ind_O = BB.Ind_O;

% Betasheet orientation info export
Structure.Extra.TransV = BB.TransV;
Structure.Extra.TwistV = BB.TwistV;
Structure.Extra.RotV   = [GUI_Inputs.Phi_D,GUI_Inputs.Psi_D,GUI_Inputs.Theta_D];

% Include the whole BB info for debug
Structure.Extra.BB = BB;

% export necessary handle and functions
Structure.hPlotFunc = @Plot_Betasheet_AmideI;
Structure.hParseGUIFunc = @ParseGUI_Betasheet;
Structure.hGUIs = hGUIs;

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
% Read GUI variables
hGUIs  = GUI_data.hGUIs;
GUI_Inputs = ParseGUI_Betasheet(hGUIs);

if isstruct(eventdata)
    if strcmp(eventdata.Source,'External')
        GUI_Inputs.External.hF  = eventdata.hF;
        GUI_Inputs.External.hAx = eventdata.hAx;
    end
end

hAx = 'New';
hF = Plot_Betasheet_AmideI(hAx,GUI_data.Structure,GUI_Inputs);

function PlotModes(hObject, eventdata, GUI_data)
Plot_Modes(GUI_data.hModel_Betasheet_AmideI);

