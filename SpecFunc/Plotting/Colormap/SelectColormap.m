function [CMAP,CMap_List] = SelectColormap(map_Index)
    
CMap_List = {'Jet',...
             'JetWhite',...
             'JetBlack',...
             'JetWhiteWide',...
             'Jet_Pale',...
             'CoolBlack'};

%ColormapName = CMap_List{map_Index};
switch map_Index
    case 0
        % doing nothing for list inquary
        CMAP = [];
    case 1
        CMAP = colormap(CMap_List{map_Index});
    otherwise
        CMAP_tmp = load(CMap_List{map_Index});
        CMAP = CMAP_tmp.MAP;
end