function AssignParent(app,ParentInfo)
switch ParentInfo{1}
    case 'COSMOSS'
        app.Parent = ParentInfo{2};
        disp('Running Model_TCO from COSMOSS...')
    case 'Comb2'
%         hModel_Comb2 = varargin{2};
%         Comb2_Order  = varargin{3};
% 
%         GUI_data.hModel_Comb2 = hModel_Comb2;
%         GUI_data.Comb2_Order  = Comb2_Order;
% 
%         % Add comb2 order # to GUI title, if necessary
%         TitleStr = hModel_TCO.Name;
%         if ~strcmp(TitleStr(1),'#')
%             hModel_TCO.Name = ['#',int2str(Comb2_Order),', ',TitleStr];
%         end
% 
%         disp('Running Model_TCO as a sub GUI of Comb2...')
end