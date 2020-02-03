function UpdateGUIs(hGUIs,Para_S,FieldName)
% This function update each GUI input (each elemets in hGUIs) with the
% corresponding parameters (Para_S, parameter structure) in the 
% corresponding filed.

Tag    =  fieldnames(   Para_S);
Type   = struct2cell(FieldName);
Para_S = struct2cell(   Para_S);

for i = 1:length(Tag)
    if strcmp(Type{i},'String')
        hGUIs.(Tag{i}).(Type{i}) = num2str(Para_S{i});
    else
        hGUIs.(Tag{i}).(Type{i}) = Para_S{i};
    end
end