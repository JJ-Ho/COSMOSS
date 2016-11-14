% ^ GUI Skeleton ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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

function Model_Comb2_OpeningFcn(hModel_Comb2, eventdata, GUI_data, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hModel (see VARARGIN)
if nargin > 3
    switch varargin{1}
        case 'COSMOSS'
            hCOSMOSS = varargin{2};
            Data_COSMOSS = guidata(hCOSMOSS);
            
            %PRE ASSIGN VALUES TO SUBSTITUTE MAIN GUI VALUES
            hGUIs_COSMOSS  = Data_COSMOSS.hGUIs;
            set(hGUIs_COSMOSS.X_Min  ,'String','1500')
            set(hGUIs_COSMOSS.X_Max  ,'String','1800')
            
            GUI_data.hCOSMOSS = hCOSMOSS;
            
            disp('Running Model_Comb2 durectly from COSMOSS...')
        case 'Comb2'
            hModel_Comb2 = varargin{2};
            Comb2_Order  = varargin{3};
            
            GUI_data.hModel_Comb2 = hModel_Comb2;
            GUI_data.Comb2_Order  = Comb2_Order;
            
            % Add comb2 order # to GUI title, if necessary
            TitleStr = hModel_Comb2.Name;
            if ~strcmp(TitleStr(1),'#')
                hModel_Comb2.Name = ['#',int2str(Comb2_Order),', ',TitleStr];
            end
            
            disp('Running Model_Comb2 as a sub GUI of Comb2...')
    end
else
    disp('Running Model_Comb2 in stand alone mode...')
end

% Call createInterface to create GUI elements
hGUIs = GUI_Comb2(hModel_Comb2);

% Prep necessary data to be export
GUI_data.hModel_Comb2 = hModel_Comb2;
GUI_data.hGUIs        = hGUIs; % export GUI handles

guidata(hModel_Comb2,GUI_data);

function varargout = Model_Comb2_OutputFcn(hModel_Comb2, eventdata, GUI_data) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = hModel_Comb2;

function Export_Handle_Callback(hObject, eventdata, GUI_data)
% export handles back to work space
assignin('base', 'Data_Comb2', GUI_data)
disp('Updated GUI Data_Comb2 exported!')
% ^ GUI Skeleton ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



function BuildModel(GUI_data, Comb2_Order)
%% Construct model GUI 
if isfield(GUI_data,'Load')
    StructModel = GUI_data.Load.StructModel;
else
    hGUIs       = GUI_data.hGUIs;
    switch Comb2_Order
        case 1
            StructModel = get(hGUIs.StructBox1,'Value');
        case 2
            StructModel = get(hGUIs.StructBox2,'Value');
    end
end 

[fhStructure,~,~] = StructureModel(StructModel);

hModel = feval(fhStructure,'Comb2',GUI_data.hModel_Comb2,Comb2_Order);

%% Update GUI_data_Model
% pass empty Structure for ???
GUI_data_Model              = guidata(hModel);
GUI_data_Model.Structure    = [];

% Update the GUI inputs if Load Structure 
if isfield(GUI_data,'Load')
    hGUIs_Model = GUI_data_Model.hGUIs;
    Para_S      = GUI_data.Load.GUI_Inputs;
    FieldName   = GUI_data.Load.GUI_FieldName;
    
    UpdateGUIs(hGUIs_Model,Para_S,FieldName)
    
    % push the loaded structure data to sub-GUI
    GUI_data_Model.GUI_FieldName = FieldName;
    GUI_data_Model.Structure     = GUI_data.Load.Structure;
end

guidata(hModel,GUI_data_Model)

%% Update GUI_Data in Comb2 and propergate to COSMOSS
% Save the modle handle to proper name
switch Comb2_Order
    case 1
        GUI_data.hStruc1 = hModel;
    case 2
        GUI_data.hStruc2 = hModel;
end

guidata(GUI_data.hModel_Comb2,GUI_data)

Export2GUIs(GUI_data_Model);

function LoadStructure(hObject, eventdata, GUI_data)
%% load previously saved comb2 output
PWD = pwd;
PDB_Path = [PWD, '/StructureFiles/Comb2/'];

[FilesName,PathName,~] = uigetfile({'*.mat','Comb2 output'; ...
                                    '*,*','All Files'},...
                                    'Select inputs',PDB_Path);
L = load([PathName FilesName]);

Structure   = L.Structure;
GUI_Inputs_0 = L.GUI_Inputs0;
GUI_Inputs_1 = L.GUI_Inputs1;
GUI_Inputs_2 = L.GUI_Inputs2;

GUI_FieldName_0 = L.GUI_FieldName0;
GUI_FieldName_1 = L.GUI_FieldName1;
GUI_FieldName_2 = L.GUI_FieldName2;

%% Update GUI front-end in Comb2
UpdateGUIs(GUI_data.hGUIs,GUI_Inputs_0,GUI_FieldName_0)

%% Call sub-model GUIs 
% construct the 1st model
GUI_data.Load.StructModel   = Structure.StrucData1.StructModel;
GUI_data.Load.GUI_Inputs    = GUI_Inputs_1;
GUI_data.Load.GUI_FieldName = GUI_FieldName_1;
GUI_data.Load.Structure     = Structure.StrucData1;
BuildModel(GUI_data, 1)

% Update GUI_data after calling Struct1
GUI_data = guidata(hObject);

% construct the 2nd model
GUI_data.Load.StructModel   = Structure.StrucData2.StructModel;
GUI_data.Load.GUI_Inputs    = GUI_Inputs_2;
GUI_data.Load.GUI_FieldName = GUI_FieldName_2;
GUI_data.Load.Structure     = Structure.StrucData2;
BuildModel(GUI_data, 2)

% Update GUI_data after calling Struct2  
GUI_data = guidata(hObject);

%% export to guidata
GUI_data.Structure  = Structure;
guidata(hObject,GUI_data)

function SaveStructure(hObject, eventdata, GUI_data)
%% Collecting outputs
hStruc1 = GUI_data.hStruc1;
hStruc2 = GUI_data.hStruc2;
GUI_Data1 = guidata(hStruc1);
GUI_Data2 = guidata(hStruc2);

GUI_Inputs0 =   GUI_data.GUI_Inputs;
GUI_Inputs1 = GUI_Data1.GUI_Inputs;
GUI_Inputs2 = GUI_Data2.GUI_Inputs;

Output.Structure   = GUI_data.Structure;
Output.GUI_Inputs0 = GUI_Inputs0;
Output.GUI_Inputs1 = GUI_Inputs1;
Output.GUI_Inputs2 = GUI_Inputs2;

Output.GUI_FieldName0 = GUI_data.GUI_FieldName;
Output.GUI_FieldName1 = GUI_data.GUI_FieldName1;
Output.GUI_FieldName2 = GUI_data.GUI_FieldName2;

%% Determine path and save
PWD = pwd;
PDB_Path = [PWD, '/StructureFiles/Comb2/'];

[FilesName,PathName,~] = uiputfile({'*.mat','Comb2 output'; ...
                                    '*,*','All Files'},...
                                    'Select inputs',PDB_Path);

save([PathName FilesName],'-struct','Output')
                                                                
function UpdateStructure(hObject, eventdata, GUI_data)
hGUIs = GUI_data.hGUIs;

hStruc1        = GUI_data.hStruc1;
hStruc2        = GUI_data.hStruc2;
StrucGUI_Data1 = guidata(hStruc1);
StrucGUI_Data2 = guidata(hStruc2);
StrucData1     = StrucGUI_Data1.Structure;
StrucData2     = StrucGUI_Data2.Structure;

%% Retreive GUI inputs
GUI_Inputs = ParseGUI_Comb2(hGUIs);

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
GUI_data.GUI_Inputs  = GUI_Inputs;
GUI_data.Structure   = Structure;

% include FieldName of GUI Inputs
[~,~,~,fhGUIParser] = StructureModel(Structure.StructModel);
[~,GUI_FieldName] = fhGUIParser(hGUIs);
GUI_data.GUI_FieldName = GUI_FieldName;

guidata(hObject,GUI_data)

% update to other GUIs
Export2GUIs(GUI_data)

disp('Structure file generated!')

function hF = PlotMolecule(hObject, eventdata, GUI_data)
%% Read GUI
hGUIs = GUI_data.hGUIs;
GUI_Inputs = ParseGUI_Comb2(hGUIs);

hF = PlotComb2(GUI_data,GUI_Inputs);

function PlotModes(hObject, eventdata, GUI_data)
Plot_Modes(GUI_data.hModel_Comb2);

