function obj_New = SD_Copy(obj)
obj_New = StructureData;

FieldName = fieldnames(obj);
N_Fields = length(FieldName);

for i = 1:N_Fields
    obj_New.(FieldName{i}) = obj.(FieldName{i});
end