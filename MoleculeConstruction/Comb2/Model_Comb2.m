function varargout = Model_Comb2(varargin)
% MODEL_COMB2 MATLAB code for Model_Comb2.fig
%      MODEL_COMB2, by itself, creates a new MODEL_COMB2 or raises the existing
%      singleton*.
%
%      H = MODEL_COMB2 returns the handle to a new MODEL_COMB2 or the handle to
%      the existing singleton*.
%
%      MODEL_COMB2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODEL_COMB2.M with the given input arguments.
%
%      MODEL_COMB2('Property','Value',...) creates a new MODEL_COMB2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Model_Comb2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Model_Comb2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Model_Comb2

% Last Modified by GUIDE v2.5 09-Dec-2015 11:11:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Model_Comb2_OpeningFcn, ...
                   'gui_OutputFcn',  @Model_Comb2_OutputFcn, ...
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


% --- Executes just before Model_Comb2 is made visible.
function Model_Comb2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Model_Comb2 (see VARARGIN)

% Choose default command line output for Model_Comb2
handles.output = hObject;

% Call createInterface to create GUI elements
GUI_Struc = GUI_Comb2(hObject);
handles.GUI_Struc = GUI_Struc; % export GUI handles to handles

% Get Main function's handles
if nargin > 3    
    if ishandle(varargin{1}) 
       hMain = varargin{1};
       Data_Main = guidata(hMain);

       handles.hMain = hMain;
       handles.Data_Main = Data_Main;
    else
        disp('Running in stand alone mode.')
    end
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Model_Comb2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Model_Comb2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function UpdateStructure(hObject, eventdata, handles)
GUI_Struc = handles.GUI_Struc;

%% Construct molecule
GUI_Inputs = ParseGUI_Comb2(GUI_Struc);
% export handles back to work space
assignin('base', 'GUI', GUI_Inputs)