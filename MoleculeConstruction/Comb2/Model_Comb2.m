function varargout = Model_Comb2(varargin)
% HMODEL MATLAB code for hModel.fig
%      HMODEL, by itself, creates a new HMODEL or raises the existing
%      singleton*.
%
%      H = HMODEL returns the handle to a new HMODEL or the handle to
%      the existing singleton*.
%
%      HMODEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HMODEL.M with the given input arguments.
%
%      HMODEL('Property','Value',...) creates a new HMODEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Model_Comb2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Model_Comb2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help hModel

% Last Modified by GUIDE v2.5 28-Feb-2016 23:01:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
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

% --- Executes just before hModel is made visible.
function Model_Comb2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hModel (see VARARGIN)

% Choose default command line output for hModel
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
    end
else
    disp('Running Model_Comb2 in stand alone mode.')    
end
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Model_Comb2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;

function Struc1(hObject, eventdata, handles)

if isfield(handles,'LoadStructModel')
    StructModel = handles.LoadStructModel;
else
    GUI_Struc   = handles.GUI_Struc;
    StructModel = get(GUI_Struc.StructBox1,'Value');
end 

[fhStructure1,~,hPlotFunc1] = StructureModel(StructModel);

guidata_Struc1 = feval(fhStructure1,'Comb2_Mode');

% pass handles of comb2 to sub-GUI so when click update stracture in
% sub-GUI, it knows where to push data to.
guidata_Struc1.hComb2      = handles.hModel;
guidata_Struc1.Comb2_Order = 1;
guidata_Struc1.Structure   = handles.Structure;
guidata(guidata_Struc1.hModel,guidata_Struc1)

% Update the GUI inputs if Load Structure 
if isfield(handles,'LoadStructModel')
    hGUI = guidata_Struc1.GUI_Struc;
    GUI_Tag  = fieldnames(handles.GUI_Inputs1);
    GUI_Type = struct2cell(handles.GUI_FieldName1);
    GUI_Para = struct2cell(handles.GUI_Inputs1);

    %ParseGUI_Comb2(hGUI,GUI_Tag,GUI_Para);
    for i = 1:length(GUI_Tag)
        if strcmp(GUI_Type{i},'String')
            hGUI.(GUI_Tag{i}).(GUI_Type{i}) = num2str(GUI_Para{i});
        else
            hGUI.(GUI_Tag{i}).(GUI_Type{i}) = GUI_Para{i};
        end
    end
end

% update guidata in comb2 and export 
handles.hStruc1    = guidata_Struc1.hModel;
handles.hPlotFunc1 = hPlotFunc1;
guidata(hObject,handles)

% update Structure in sub-GUI
hSubGUI = guidata_Struc1.hModel;
fhStructure1('UpdateStructure',hSubGUI,eventdata,guidata(hSubGUI));

function Struc2(hObject, eventdata, handles)

if isfield(handles,'LoadStructModel')
    StructModel = handles.LoadStructModel;
else
    GUI_Struc   = handles.GUI_Struc;
    StructModel = get(GUI_Struc.StructBox2,'Value');
end 

[fhStructure2,~,hPlotFunc2] = StructureModel(StructModel);

guidata_Struc2 = feval(fhStructure2,'Comb2_Mode');

% pass handles of comb2 to sub-GUI so when click update stracture in
% sub-GUI, it knows where to push data to.
guidata_Struc2.hComb2      = handles.hModel;
guidata_Struc2.Comb2_Order = 2;
guidata_Struc2.Structure   = handles.Structure;
guidata(guidata_Struc2.hModel,guidata_Struc2)

% Update the GUI inputs if Load Structure 
if isfield(handles,'LoadStructModel')
    hGUI = guidata_Struc2.GUI_Struc;
    GUI_Tag  = fieldnames(handles.GUI_Inputs2);
    GUI_Type = struct2cell(handles.GUI_FieldName2);
    GUI_Para = struct2cell(handles.GUI_Inputs2);

    %ParseGUI_Comb2(hGUI,GUI_Tag,GUI_Para);
    for i = 1:length(GUI_Tag)
        if strcmp(GUI_Type{i},'String')
            hGUI.(GUI_Tag{i}).(GUI_Type{i}) = num2str(GUI_Para{i});
        else
            hGUI.(GUI_Tag{i}).(GUI_Type{i}) = GUI_Para{i};
        end
    end
end

handles.hStruc2    = guidata_Struc2.hModel;
handles.hPlotFunc2 = hPlotFunc2;
guidata(hObject,handles)

% update Structure in sub-GUI
hSubGUI = guidata_Struc2.hModel;
fhStructure2('UpdateStructure',hSubGUI,eventdata,guidata(hSubGUI));

function LoadStructure(hObject, eventdata, handles)
%% load previously saved comb2 output
PWD = pwd;
PDB_Path = [PWD, '/StructureFiles/Comb2/'];

[FilesName,PathName,~] = uigetfile({'*.mat','Comb2 output'; ...
                                    '*,*','All Files'},...
                                    'Select inputs',PDB_Path);
L = load([PathName FilesName]);

Structure   = L.Structure;
GUI_Inputs0 = L.GUI_Inputs0;
GUI_Inputs1 = L.GUI_Inputs1;
GUI_Inputs2 = L.GUI_Inputs2;

GUI_FieldName0 = L.GUI_FieldName0;
GUI_FieldName1 = L.GUI_FieldName1;
GUI_FieldName2 = L.GUI_FieldName2;

%% Update GUI front end in Comb2
hGUI = handles.GUI_Struc;
GUI_Tag  = fieldnames(GUI_Inputs0);
GUI_Type = struct2cell(GUI_FieldName0);
GUI_Para = struct2cell(GUI_Inputs0);

%ParseGUI_Comb2(hGUI,GUI_Tag,GUI_Para);
for i = 1:length(GUI_Tag)
    if strcmp(GUI_Type{i},'String')
        hGUI.(GUI_Tag{i}).(GUI_Type{i}) = num2str(GUI_Para{i});
    else
        hGUI.(GUI_Tag{i}).(GUI_Type{i}) = GUI_Para{i};
    end
end

%% Call sub-model GUIs 
% construct the 1st model
handles.LoadStructModel = Structure.StrucData1.StructModel;
handles.GUI_Inputs1     = GUI_Inputs1;
handles.GUI_FieldName1  = GUI_FieldName1;
handles.Structure       = Structure.StrucData1;
Struc1(hObject, eventdata, handles)

% construct the 2nd model
handles.LoadStructModel = Structure.StrucData2.StructModel;
handles.GUI_Inputs2     = GUI_Inputs2;
handles.GUI_FieldName2  = GUI_FieldName2;
handles.Structure       = Structure.StrucData2;
Struc2(hObject, eventdata, handles)

%% export to guidata
handles.Structure  = Structure;
guidata(hObject,handles)

% update Structure in sub-GUI 
% for some reason UpdateStructure does not see handles.hStruc1 paased
% after excuted the UpdateStructure in sub-GUIs
% UpdateStructure(hObject, eventdata, guidata(hObject))

function SaveStructure(hObject, eventdata, handles)
%% Collecting outputs
% GUI inputs from both sub GUIs
hStruc1 = handles.hStruc1;
hStruc2 = handles.hStruc2;
GUI_Data1 = guidata(hStruc1);
GUI_Data2 = guidata(hStruc2);

GUI_Inputs0 =   handles.GUI_Inputs;
GUI_Inputs1 = GUI_Data1.GUI_Inputs;
GUI_Inputs2 = GUI_Data2.GUI_Inputs;

Output.Structure   = handles.Structure;
Output.GUI_Inputs0 = GUI_Inputs0;
Output.GUI_Inputs1 = GUI_Inputs1;
Output.GUI_Inputs2 = GUI_Inputs2;

Output.GUI_FieldName0 = handles.GUI_FieldName;
Output.GUI_FieldName1 = handles.GUI_FieldName1;
Output.GUI_FieldName2 = handles.GUI_FieldName2;

%% Determine path and save
PWD = pwd;
PDB_Path = [PWD, '/StructureFiles/Comb2/'];

[FilesName,PathName,~] = uiputfile({'*.mat','Comb2 output'; ...
                                    '*,*','All Files'},...
                                    'Select inputs',PDB_Path);

save([PathName FilesName],'-struct','Output')
                                                                
function UpdateStructure(hObject, eventdata, handles)
GUI_Struc = handles.GUI_Struc;

hStruc1        = handles.hStruc1;
hStruc2        = handles.hStruc2;
StrucGUI_Data1 = guidata(hStruc1);
StrucGUI_Data2 = guidata(hStruc2);
StrucData1     = StrucGUI_Data1.Structure;
StrucData2     = StrucGUI_Data2.Structure;

%% Retreive GUI inputs
GUI_Inputs = ParseGUI_Comb2(GUI_Struc);

Conc_Scaling = GUI_Inputs.Conc_Scaling;
Trans_X      = GUI_Inputs.Trans_X;
Trans_Y      = GUI_Inputs.Trans_Y;
Trans_Z      = GUI_Inputs.Trans_Z;
Rot_Phi      = GUI_Inputs.Rot_Phi/180*pi;
Rot_Psi      = GUI_Inputs.Rot_Psi/180*pi;
Rot_Theta    = GUI_Inputs.Rot_Theta/180*pi;

TransV = [Trans_X,Trans_Y,Trans_Z];
% RM = R1_ZYZ_0(Rot_Phi,Rot_Psi,Rot_Theta);
RM = Rx(Rot_Phi)*Ry(Rot_Psi)*Rz(Rot_Theta);

%% Shift the center of mass of each structure to origin
Center1 = StrucData1.center;
% COM1 = sum(Center1,1)./size(Center1,1);
COM1 = [0,0,0];

Center2 = StrucData2.center;
% COM2 = sum(Center2,1)./size(Center2,1);
COM2 = [0,0,0];

%% Move the second molecule and merge
% Center 
Center1_0 = bsxfun(@minus,Center1,COM1);
Center2_0 = bsxfun(@minus,Center2,COM2);

XYZ2_R = (RM*Center2_0')';
Center2_T = bsxfun(@plus,XYZ2_R,TransV);

M_Center = [Center1_0; Center2_T];

% mu 
TDV1 = StrucData1.mu;
TDV2 = StrucData2.mu;

TDV2_R = Conc_Scaling .* (RM*TDV2')';

M_TDV = [TDV1;TDV2_R];

% alpha matrix
Raman_Matrix1 = StrucData1.alpha_matrix;
Raman_Matrix2 = StrucData2.alpha_matrix;
Num_Modes2    = StrucData2.Num_Modes;

Raman_Matrix2_R = zeros(size(Raman_Matrix2));
for i=1:Num_Modes2
    Raman_Matrix2_R(i,:,:) = Conc_Scaling .* RM*squeeze(Raman_Matrix2(i,:,:))*RM';
end

M_Raman_Matrix = cat(1,Raman_Matrix1,Raman_Matrix2_R);

% alpha vector
Raman1   = reshape(Raman_Matrix1  ,[],9);
Raman2_R = reshape(Raman_Matrix2_R,[],9);
M_Raman  = reshape(M_Raman_Matrix ,[],9);

% XYZ 
XYZ1 = StrucData1.XYZ;
XYZ2 = StrucData2.XYZ;

XYZ1_0 = bsxfun(@minus,XYZ1,COM1);
XYZ2_0 = bsxfun(@minus,XYZ2,COM2);

XYZ2_R = (RM*XYZ2_0')';
XYZ2_T = bsxfun(@plus,XYZ2_R,TransV);

M_XYZ = [XYZ1_0;XYZ2_T];

%% update Structure1 and 2

StrucData1.center       = Center1_0;
StrucData1.mu           = TDV1;
StrucData1.alpha_matrix = Raman_Matrix1;
StrucData1.alpha        = Raman1;
StrucData1.XYZ          = XYZ1_0;

StrucData2.center       = Center2_T;
StrucData2.mu           = TDV2_R;
StrucData2.alpha_matrix = Raman_Matrix2_R;
StrucData2.alpha        = Raman2_R;
StrucData2.XYZ          = XYZ2_T;

%% Merge each catagory that doesn't need to move
% freq
Freq1 = StrucData1.freq;
Freq2 = StrucData2.freq;
M_Freq = [Freq1;Freq2];

% anharmonicity
Anharm1 = StrucData1.anharm;
Anharm2 = StrucData2.anharm;
M_Anharm = [Anharm1;Anharm2];

% % Atom Serial Number
% EAtomSerNo  = StrucData1.AtomSerNo;
% Shift_Index = size(EXYZ,1);
% FAtomSerNo  = StrucData2.AtomSerNo + Shift_Index;
% M_AtomSerNo = [EAtomSerNo;FAtomSerNo];

%% Output
Structure.center       = M_Center;
Structure.freq         = M_Freq;
Structure.anharm       = M_Anharm;
Structure.mu           = M_TDV;
Structure.alpha        = M_Raman;
Structure.alpha_matrix = M_Raman_Matrix;
Structure.Num_Modes    = size(M_TDV,1);
Structure.XYZ          = M_XYZ;
Structure.FilesName    = 'Comb2';

% Export into Structure so it can be passsed around different GUIs
Structure.StrucData1  = StrucData1;
Structure.StrucData2  = StrucData2;
Structure.StructModel = 5;

%% export back to handles and Main GUI if any
handles.GUI_Inputs  = GUI_Inputs;
handles.Structure   = Structure;

% include FieldName of GUI Inputs
[~,~,~,hGUIParser] = StructureModel(Structure.StructModel);
[~,GUI_FieldName] = hGUIParser(GUI_Struc);
handles.GUI_FieldName = GUI_FieldName;

guidata(hObject,handles)

% update to other GUIs
Export2GUIs(handles)

disp('Structure file generated!')

function hF = PlotMolecule(hObject, eventdata, handles)
GUI_Struc = handles.GUI_Struc;
GUI_Inputs = ParseGUI_Comb2(GUI_Struc);

hF = PlotComb2(handles,GUI_Inputs);

UpdateStructure(hObject, eventdata, handles)

function PlotModes(hObject, eventdata, handles)
Plot_Modes(handles.hModel);

function Export_Handle_Callback(hObject, eventdata, handles)
% export handles back to work space
assignin('base', 'hModel_Comb2', handles)
disp('Updated handles exported!')
