function obj_New = SD_Copy(obj_SD)
% This method copy the input StructureData as a New one
obj_New = StructureData;

FieldName = fieldnames(obj_SD);
N_Fields = length(FieldName);


for i = 1:N_Fields
    mp = findprop(obj_SD,FieldName{i});
    if ~mp.NonCopyable
        % Check if the property is a NonCopyable one
        obj_New.(FieldName{i}) = obj_SD.(FieldName{i});
    else
        %disp(['Property: ', FieldName{i}, ' is not copied'])
    end
end