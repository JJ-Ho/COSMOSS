function GUI = GUI_PDB_AmideI(fig) 
%% Add base layout
MainLayout = uix.VBoxFlex(...
    'Parent',fig,...
    'Spacing',1,...
    'Padding',5);

%% Title
uicontrol('Parent',MainLayout,...
    'Style','text',...
    'String','Rotation',...
    'BackgroundColor',[0,1,1]);

%% oreintation inputs
Rot_HBox = uix.HBoxFlex(...
    'Parent',MainLayout,...
    'Spacing',3);

    uicontrol('Parent',Rot_HBox,...
        'Style','text',...
        'String','Phi');
    
    uicontrol('Parent',Rot_HBox,...
        'Style','edit',...
        'String','0',...
        'Tag','Phi',...
        'Callback',@(hObject,eventdata)Model_PDB_AmideI('UpdateStructure',hObject,eventdata,guidata(hObject)));

    uicontrol('Parent',Rot_HBox,...
        'Style','text',...
        'String','Psi');
    
    uicontrol('Parent',Rot_HBox,...
        'Style','edit',...
        'String','0',...
        'Tag','Psi',...
        'Callback',@(hObject,eventdata)Model_PDB_AmideI('UpdateStructure',hObject,eventdata,guidata(hObject)));

    uicontrol('Parent',Rot_HBox,...
        'Style','text',...
        'String','Theta');

    uicontrol('Parent',Rot_HBox,...
        'Style','edit',...
        'String','0',...
        'Tag','Theta',...
        'Callback',@(hObject,eventdata)Model_PDB_AmideI('UpdateStructure',hObject,eventdata,guidata(hObject)));

% %% Translation inputs
% 
% uicontrol('Parent',MainLayout,...
%     'Style','text',...
%     'String','Translation vector',...
%     'BackgroundColor',[0,1,1]);
% 
% uicontrol('Parent',MainLayout,...
%     'Style','edit',...
%     'String','5,0,0',...
%     'Tag','Trans',...
%     'Callback',@(hObject,eventdata)Model_PDB_AmideI('UpdateStructure',hObject,eventdata,guidata(hObject)));

%% Show PDB name
uicontrol('Parent',MainLayout,...
    'Style','text',...
    'String','');

uicontrol('Parent',MainLayout,...
    'Style','text',...
    'Tag','PDB_Name',...
    'String','PDB not loaded...');

%% Button
uicontrol( 'Parent', MainLayout, ...
'String', 'Load PDB',...
'Callback',@(hObject,eventdata)Model_PDB_AmideI('LoadStructure',hObject,eventdata,guidata(hObject)));

uicontrol( 'Parent', MainLayout, ...
'String', 'Generate',...
'Callback',@(hObject,eventdata)Model_PDB_AmideI('UpdateStructure',hObject,eventdata,guidata(hObject)));

uicontrol( 'Parent', MainLayout, ...
'String', 'Plot Molecule',...
'Callback',@(hObject,eventdata)Model_PDB_AmideI('PlotMolecule',hObject,eventdata,guidata(hObject)));

%% set sizes
set(MainLayout,'Height',[15,25,5,25,-1,-1,-1])

%% output handles
GUI = guihandles(fig);
    