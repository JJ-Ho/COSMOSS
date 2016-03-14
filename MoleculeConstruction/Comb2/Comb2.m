function Structure = Comb2(StrucData1,StrucData2,GUI_Inputs)

Conc_Scaling = GUI_Inputs.Conc_Scaling;
Trans_X      = GUI_Inputs.Trans_X;
Trans_Y      = GUI_Inputs.Trans_Y;
Trans_Z      = GUI_Inputs.Trans_Z;
Rot_Phi      = GUI_Inputs.Rot_Phi/180*pi;
Rot_Psi      = GUI_Inputs.Rot_Psi/180*pi;
Rot_Theta    = GUI_Inputs.Rot_Theta/180*pi;

TransV = [Trans_X,Trans_Y,Trans_Z];
% RM = R1_ZYZ_0(Rot_Phi,Rot_Psi,Rot_Theta);
RM = Rx(Rot_Phi)*Ry(Rot_Psi)*Rz(Rot_Theta);

%% Shift the center of mass of each structure to origin
Center1 = StrucData1.center;
% COM1 = sum(Center1,1)./size(Center1,1);
COM1 = [0,0,0];

Center2 = StrucData2.center;
% COM2 = sum(Center2,1)./size(Center2,1);
COM2 = [0,0,0];

%% Move the second molecule and merge
% Center 
Center1_0 = bsxfun(@minus,Center1,COM1);
Center2_0 = bsxfun(@minus,Center2,COM2);

XYZ2_R = (RM*Center2_0')';
Center2_T = bsxfun(@plus,XYZ2_R,TransV);

M_Center = [Center1_0; Center2_T];

% mu 
TDV1 = StrucData1.mu;
TDV2 = StrucData2.mu;

TDV2_R = Conc_Scaling .* (RM*TDV2')';

M_TDV = [TDV1;TDV2_R];

% alpha matrix
Raman_Matrix1 = StrucData1.alpha_matrix;
Raman_Matrix2 = StrucData2.alpha_matrix;
Num_Modes2    = StrucData2.Num_Modes;

Raman_Matrix2_R = zeros(size(Raman_Matrix2));
for i=1:Num_Modes2
    Raman_Matrix2_R(i,:,:) = Conc_Scaling .* RM*squeeze(Raman_Matrix2(i,:,:))*RM';
end

M_Raman_Matrix = cat(1,Raman_Matrix1,Raman_Matrix2_R);

% alpha vector
Raman1   = reshape(Raman_Matrix1  ,[],9);
Raman2_R = reshape(Raman_Matrix2_R,[],9);
M_Raman  = reshape(M_Raman_Matrix ,[],9);

% XYZ 
XYZ1 = StrucData1.XYZ;
XYZ2 = StrucData2.XYZ;

XYZ1_0 = bsxfun(@minus,XYZ1,COM1);
XYZ2_0 = bsxfun(@minus,XYZ2,COM2);

XYZ2_R = (RM*XYZ2_0')';
XYZ2_T = bsxfun(@plus,XYZ2_R,TransV);

M_XYZ = [XYZ1_0;XYZ2_T];

%% update Structure1 and 2

StrucData1.center       = Center1_0;
StrucData1.mu           = TDV1;
StrucData1.alpha_matrix = Raman_Matrix1;
StrucData1.alpha        = Raman1;
StrucData1.XYZ          = XYZ1_0;

StrucData2.center       = Center2_T;
StrucData2.mu           = TDV2_R;
StrucData2.alpha_matrix = Raman_Matrix2_R;
StrucData2.alpha        = Raman2_R;
StrucData2.XYZ          = XYZ2_T;

%% Merge each catagory that doesn't need to move
% freq
Freq1 = StrucData1.freq;
Freq2 = StrucData2.freq;
M_Freq = [Freq1;Freq2];

% anharmonicity
Anharm1 = StrucData1.anharm;
Anharm2 = StrucData2.anharm;
M_Anharm = [Anharm1;Anharm2];

% % Atom Serial Number
% EAtomSerNo  = StrucData1.AtomSerNo;
% Shift_Index = size(EXYZ,1);
% FAtomSerNo  = StrucData2.AtomSerNo + Shift_Index;
% M_AtomSerNo = [EAtomSerNo;FAtomSerNo];

%% Output
Structure.center       = M_Center;
Structure.freq         = M_Freq;
Structure.anharm       = M_Anharm;
Structure.mu           = M_TDV;
Structure.alpha        = M_Raman;
Structure.alpha_matrix = M_Raman_Matrix;
Structure.Num_Modes    = size(M_TDV,1);
Structure.XYZ          = M_XYZ;
Structure.FilesName    = 'Comb2';

% Export into Structure so it can be passsed around different GUIs
Structure.StrucData1  = StrucData1;
Structure.StrucData2  = StrucData2;
Structure.StructModel = 5;