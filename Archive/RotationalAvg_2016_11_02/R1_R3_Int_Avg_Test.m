%% Generation of non-averaged R1 to R3
clear all

%% R1 - R3

tic
syms Phi Psi Theta

Rxx = -sin(Psi)*sin(Phi) + cos(Theta)*cos(Psi)*cos(Phi);
Rxy = -cos(Psi)*sin(Phi) - cos(Theta)*sin(Psi)*cos(Phi);
Rxz =  sin(Theta)*cos(Phi);
Ryx =  sin(Psi)*cos(Phi) + cos(Theta)*cos(Psi)*sin(Phi);
Ryy =  cos(Psi)*cos(Phi) - cos(Theta)*sin(Psi)*sin(Phi);
Ryz =  sin(Theta)*sin(Phi);
Rzx = -sin(Theta)*cos(Psi);
Rzy =  sin(Theta)*sin(Psi);
Rzz =  cos(Theta);

% Kron version for permutation 
R1_ZYZ = [Rxx,Rxy,Rxz; Ryx,Ryy,Ryz; Rzx,Rzy,Rzz];
disp('R1_ZYZ matrix generated...')

R2_ZYZ = kron(R1_ZYZ,R1_ZYZ);
disp('R2_ZYZ matrix generated...')

R3_ZYZ = kron(R1_ZYZ,R2_ZYZ);
disp('R3_ZYZ matrix generated...')

disp('All symbolic matrixies generated...')
toc

%% Integration
tic
R1_ZYZ_1   = int(R1_ZYZ               ,Phi  ,0,2*pi)./(2*pi);
disp('integration of R1_ZYZ_1 is done!')
toc

tic
R2_ZYZ_1   = int(R2_ZYZ               ,Phi  ,0,2*pi)./(2*pi);
disp('integration of R2_ZYZ_1 is done!')
toc

tic
R3_ZYZ_1   = int(R3_ZYZ               ,Phi  ,0,2*pi)./(2*pi);
disp('integration of R3_ZYZ_1 is done!')
toc


%% Translate and svae all symbolic matrixes of RRRRR into matlab functions

tic;matlabFunction(R1_ZYZ_1,'file','R1_ZYZ_1');disp('R1_ZYZ_1 is saved as Matlab function!');toc
tic;matlabFunction(R2_ZYZ_1,'file','R2_ZYZ_1');disp('R2_ZYZ_1 is saved as Matlab function!');toc
tic;matlabFunction(R3_ZYZ_1,'file','R3_ZYZ_1');disp('R3_ZYZ_1 is saved as Matlab function!');toc

% tic;matlabFunction(R1_ZYZ,'file','R1_ZYZ_0');disp('R1_ZYZ_0 is saved as Matlab function!');toc
% tic;matlabFunction(R2_ZYZ,'file','R2_ZYZ_0');disp('R2_ZYZ_0 is saved as Matlab function!');toc
% tic;matlabFunction(R3_ZYZ,'file','R3_ZYZ_0');disp('R3_ZYZ_0 is saved as Matlab function!');toc

disp('R1-R3 Function generation finished!')


