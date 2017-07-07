function hFcomb2 = PlotComb2(GUI_data)
%% retreive data from handles
% hStruc1 = GUI_data.hStruc1;
% hStruc2 = GUI_data.hStruc2;
% S1      = GUI_data.S1;
% S2      = GUI_data.S2;
S1 = GUI_data.Structure.Children(1);
S2 = GUI_data.Structure.Children(2);

%% make figure
hPlotFunc1 = S1.hPlotFunc;
hPlotFunc2 = S2.hPlotFunc;

GUI_Data1 = feval(S1.hParseGUIFunc,S1.hGUIs);
GUI_Data2 = feval(S2.hParseGUIFunc,S2.hGUIs);

hF1 = feval(hPlotFunc1,S1,GUI_Data1);
hF2 = feval(hPlotFunc2,S2,GUI_Data2);

% [fhFunc_Model1,~,~] = StructureModel(S1.StructModel);
% [fhFunc_Model2,~,~] = StructureModel(S2.StructModel);

% retrieve GUI data of each selected model and update he corresponding
% structure.
% GUI_Data1 = guidata(hStruc1);
% GUI_Data2 = guidata(hStruc2);
% GUI_Data1.Structure = S1;
% GUI_Data2.Structure = S2;

% hF1 = fhFunc_Model1('PlotMolecule',hStruc1,'',GUI_Data1);
% hF2 = fhFunc_Model2('PlotMolecule',hStruc2,'',GUI_Data2);


hAx1 = findobj(hF1,'type','axes');
hAx2 = findobj(hF2,'type','axes');

hFcomb2 = figure;
hAx_Comb2 = axes;
copyobj(allchild(hAx1),hAx_Comb2);
copyobj(allchild(hAx2),hAx_Comb2);

close(hF1)
close(hF2)

%% figure options
axis image;
rotate3d on
grid on

xlabel('X')
ylabel('Y')
view([0,90])