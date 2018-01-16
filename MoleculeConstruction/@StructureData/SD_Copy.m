function obj_New = SD_Copy(obj_SD)
obj_New = StructureData;

FieldName = fieldnames(obj_SD);
N_Fields = length(FieldName);

for i = 1:N_Fields
    obj_New.(FieldName{i}) = obj_SD.(FieldName{i});
end