function obj_Framed = SD_SetFrame(obj_SD,Center_Ind,Z_Ind,XZ_Ind)
% Given the desired frame for the molecule to align with, this method put
% the StructureData into aligned orentation.
%% Rotated the molecule with [Orientation] info
% Get Ratational Matrix to rotate the molecule to defined frame 
Vec_Z  = obj_SD.XYZ( Z_Ind(2),:) - obj_SD.XYZ( Z_Ind(1),:);
Z      = Vec_Z/norm(Vec_Z);

Vec_XZ = obj_SD.XYZ(XZ_Ind(2),:) - obj_SD.XYZ(XZ_Ind(1),:);
XZ     = Vec_XZ/norm(Vec_XZ);

Y      = cross(Z,XZ);
Y      = Y/norm(Y);

X      = cross(Y,Z);
X      = X/norm(X);

% Generate rotational matrix that transfer molecule from G09 simulation 
% frame to the defined Molecular Frame (MF)
Old_Frame      = zeros(3,3);
Old_Frame(:,1) = X;
Old_Frame(:,2) = Y;
Old_Frame(:,3) = Z;
New_Frame      = eye(3);
Rot2NF         = Euler_Rot(New_Frame,Old_Frame);

obj_RF = SD_Rot(obj_SD,Rot2NF);

%% Set the center Atom to zero 
% Move Center to [0,0,0]
% Note: Center_Ind can be more than 1 atom, for example, the center of
% Bezene will be defined as the avarage position of C1-C6, thus the 
% Center_index would be 1:6
% 
%   C1=C2
%  /     \
% C6  X  C3
%  \\   //
%   C5-C4
% 
Center = sum(obj_RF.XYZ(Center_Ind,:),1)./length(Center_Ind);
obj_Framed = SD_Trans(obj_RF,-Center);