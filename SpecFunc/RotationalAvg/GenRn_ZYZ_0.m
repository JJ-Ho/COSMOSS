%% Generation of non-averaged R1 to R5
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

% Original version for permutation 
R1_ZYZ = [Rxx,Rxy,Rxz; Ryx,Ryy,Ryz; Rzx,Rzy,Rzz];
disp('R1_ZYZ matrix generated...')

R2_ZYZ = kron(R1_ZYZ,R1_ZYZ);
disp('R2_ZYZ matrix generated...')

R3_ZYZ = kron(R1_ZYZ,R2_ZYZ);
disp('R3_ZYZ matrix generated...')

disp('All symbolic matrixies generated...')
toc


%% Translate and svae all symbolic matrixes of RRRRR into matlab functions
tic;matlabFunction(R1_ZYZ,'file','R1_ZYZ_0_no_opt','Optimize',false);disp('R1_ZYZ_0_no_opt is saved as Matlab function!');toc
tic;matlabFunction(R2_ZYZ,'file','R2_ZYZ_0_no_opt','Optimize',false);disp('R2_ZYZ_0_no_opt is saved as Matlab function!');toc
tic;matlabFunction(R3_ZYZ,'file','R3_ZYZ_0_no_opt','Optimize',false);disp('R3_ZYZ_0_no_opt is saved as Matlab function!');toc

tic;matlabFunction(R1_ZYZ,'file','R1_ZYZ_0');disp('R1_ZYZ_0 is saved as Matlab function!');toc
tic;matlabFunction(R2_ZYZ,'file','R2_ZYZ_0');disp('R2_ZYZ_0 is saved as Matlab function!');toc
tic;matlabFunction(R3_ZYZ,'file','R3_ZYZ_0');disp('R3_ZYZ_0 is saved as Matlab function!');toc

disp('R1-R3 Function generation finished!')


%% R4
R4_ZYZ = kron(R1_ZYZ,R3_ZYZ);
disp('R4_ZYZ matrix generated...')

tic;matlabFunction(R4_ZYZ,'file','R4_ZYZ_0_no_opt','Optimize',false);disp('R4_ZYZ_0_no_opt is saved as Matlab function!');toc
tic;matlabFunction(R4_ZYZ,'file','R4_ZYZ_0');disp('R4_ZYZ_0 is saved as Matlab function!');toc

%% R5
R5_ZYZ = kron(R1_ZYZ,R4_ZYZ);
disp('R5_ZYZ matrix generated...')

tic;matlabFunction(R5_ZYZ,'file','R5_ZYZ_0_no_opt','Optimize',false);disp('R5_ZYZ_0_no_opt is saved as Matlab function!');toc
tic;matlabFunction(R5_ZYZ,'file','R5_ZYZ_0');disp('R5_ZYZ_0 is saved as Matlab function!');toc
