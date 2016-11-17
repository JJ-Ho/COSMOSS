%% Generation of conversion matrix that convert Lab XYZ frame to PS frame
% Jones reflection Geometry
% 
% For transimissive signial: 
%           x      y    z
%       p  Cos(a)  0   Sin(a)
%       s   0      1    0
% 
% 
% For reflective signial: 
%           x      y    z
%       p -Cos(a)  0   Sin(a)
%       s   0      1    0
% 
% 130909 JJHo

% define incident angle of Pump, Probe, Visible probe, and Signal
% all angle is defined by angle between surface normal and the light beams.
SaveAsFunc = 'n';

syms A_Pu1 A_Pu2 A_Pr A_Vi A_Si

Pu1  = sym('Pu%d%d' ,[2,3]);
Pu2  = sym('Pu%d%d' ,[2,3]);
Pr   = sym('Pr%d%d' ,[2,3]);
Vi   = sym('Vi%d%d' ,[2,3]);
Si_R = sym('SiR%d%d',[2,3]);
Si_T = sym('SiT%d%d',[2,3]);

%% fill in numbers
Pu1(1,1) = cos(A_Pu1);
Pu1(1,2) = 0;
Pu1(1,3) = sin(A_Pu1);
Pu1(2,1) = 0;
Pu1(2,2) = 1;
Pu1(2,3) = 0;

Pu2(1,1) = cos(A_Pu2);
Pu2(1,2) = 0;
Pu2(1,3) = sin(A_Pu2);
Pu2(2,1) = 0;
Pu2(2,2) = 1;
Pu2(2,3) = 0;

Pr(1,1) = cos(A_Pr);
Pr(1,2) = 0;
Pr(1,3) = sin(A_Pr);
Pr(2,1) = 0;
Pr(2,2) = 1;
Pr(2,3) = 0;

Vi(1,1) = cos(A_Vi);
Vi(1,2) = 0;
Vi(1,3) = sin(A_Vi);
Vi(2,1) = 0;
Vi(2,2) = 1;
Vi(2,3) = 0;

Si_R(1,1) = -cos(A_Si); % reflective geometry
Si_R(1,2) = 0;
Si_R(1,3) = sin(A_Si);
Si_R(2,1) = 0;
Si_R(2,2) = 1;
Si_R(2,3) = 0;

Si_T(1,1) = cos(A_Si); % Transmissive geometry
Si_T(1,2) = 0;
Si_T(1,3) = sin(A_Si);
Si_T(2,1) = 0;
Si_T(2,2) = 1;
Si_T(2,3) = 0;


%% Direct product and permutation of index
Rps3 =           kron(kron(Si_R,Vi),Pu1);% SFG
Rps4 =      kron(kron(kron(Si_T,Pr),Pu2),Pu1);% 2DIR
Rps5 = kron(kron(kron(kron(Si_R,Vi),Pr ),Pu2),Pu1);% 2DSFG

%% 
if strcmp(SaveAsFunc,'y')   
    matlabFunction(Rps3,'file','JonesRef3');
    matlabFunction(Rps4,'file','JonesTrans4');
    matlabFunction(Rps5,'file','JonesRef5');
end