function hFcomb2 = PlotComb2(GUI_data,GUI_Inputs)
%% retreive data from handles
hStruc1 = GUI_data.hStruc1;
hStruc2 = GUI_data.hStruc2;
S1      = GUI_data.S1;
S2      = GUI_data.S2;

%% Retreive GUI inputs
% % GUI_Inputs = ParseGUI_Comb2(GUI_Struc);
% 
% % translation part has been taking cared by UpdateStructure in Model_Comb2
% % Trans_X    = GUI_Inputs.Trans_X;
% % Trans_Y    = GUI_Inputs.Trans_Y;
% % Trans_Z    = GUI_Inputs.Trans_Z;
% % TransV = [Trans_X,Trans_Y,Trans_Z];
% 
% % Rotation and center of axes of the second molecule
% Rot_Phi    = GUI_Inputs.Rot_Phi/180*pi;
% Rot_Psi    = GUI_Inputs.Rot_Psi/180*pi;
% Rot_Theta  = GUI_Inputs.Rot_Theta/180*pi;
% RM = Rx(Rot_Phi)*Ry(Rot_Psi)*Rz(Rot_Theta);
% 
% Center2 = S2.center;
% COM2 = sum(Center2,1)./size(Center2,1);

%% make figure
[fhFunc_Model1,~,~] = StructureModel(S1.StructModel);
[fhFunc_Model2,~,~] = StructureModel(S2.StructModel);

% retrieve GUI data of each selected model and update he corresponding
% structure.
GUI_Data1 = guidata(hStruc1);
GUI_Data2 = guidata(hStruc2);
GUI_Data1.Structure = S1;
GUI_Data2.Structure = S2;

hF1 = fhFunc_Model1('PlotMolecule',hStruc1,'',GUI_Data1);
hF2 = fhFunc_Model2('PlotMolecule',hStruc2,'',GUI_Data2);

hAx1 = findobj(hF1,'type','axes');
hAx2 = findobj(hF2,'type','axes');

hFcomb2 = figure;
hAx_Comb2 = axes;
copyobj(allchild(hAx1),hAx_Comb2);
copyobj(allchild(hAx2),hAx_Comb2);

close(hF1)
close(hF2)
hold on 

% plot the XYZ axis of the second structure
% Scale = 2;
% Axis2_0 = [1,0,0;0,1,0;0,0,1];
% Axis2_R = RM*Axis2_0 .* Scale;
% 
% Origin   = bsxfun(@times,COM2,[1,1,1]');
% 
% CMatrix = [1,0,0;0,1,0;0,0,1];
% 
% for i = 1:3
% quiver3(Origin(i,1),Origin(i,2),Origin(i,3),...
%         Axis2_R(i,1),Axis2_R(i,2),Axis2_R(i,3),0,...
%         'LineWidth',5,...
%         'Color',CMatrix(i,:));
% end  
    

hold off

%% figure options
axis image;
rotate3d on
grid on

xlabel('X')
ylabel('Y')
view([0,90])