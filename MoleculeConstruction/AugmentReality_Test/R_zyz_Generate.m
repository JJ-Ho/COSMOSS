%% Generation of non-averaged R1 to R5
clear all
%%
% Elapsed time is 1.371786 seconds.
tic
syms D_X D_Y D_Z

c1= cos(D_X);
s1= sin(D_X);
c2= cos(D_Y);
s2=-sin(D_Y);
c3= cos(D_Z);
s3=-sin(D_Z);

% Original version for permutation 
R1_zyx=[ c1*c2, c1*s2*s3 -    c3*s1,    s1*s3 + c1*c3*s2;
        c2*s1,    c1*c3 + s1*s2*s3, c3*s1*s2 -    c1*s3;
          -s2,               c2*s3,               c2*c3];

disp('R1_ZYZ matrix generated...')
tic;matlabFunction(R1_zyx,'file','R1_zyx_0');disp('R1_zyx_0 is saved as Matlab function!');toc

R2_zyx = outer(R1_zyx,R1_zyx,0);
R2_zyx = permute(R2_zyx,[1,3,2,4]);
R2_zyx = reshape(R2_zyx,[9,9]);
R2_zyx = simplify(R2_zyx,'steps',100);
disp('R2_ZYZ matrix generated...')
tic;matlabFunction(R2_zyx,'file','R2_zyx_0');disp('R2_zyx_0 is saved as Matlab function!');toc


R3_zyx = outer(R2_zyx,R1_zyx,0);
R3_zyx = permute(R3_zyx,[1,2,5,4,3,6]);
R3_zyx = permute(R3_zyx,[1,2,3,5,4,6]);
R3_zyx = reshape(R3_zyx,[27,27]);
R3_zyx = simplify(R3_zyx,'steps',100);
disp('R3_ZYZ matrix generated...')
tic;matlabFunction(R3_zyx,'file','R3_zyx_0');disp('R3_zyx_0 is saved as Matlab function!');toc

R4_zyx = outer(R3_zyx,R1_zyx,0);
R4_zyx = permute(R4_zyx,[1,2,3,7,5,6,4,8]);
R4_zyx = permute(R4_zyx,[1,2,3,4,5,7,6,8]);
R4_zyx = permute(R4_zyx,[1,2,3,4,6,5,7,8]);
R4_zyx = reshape(R4_zyx,[81,81]);
R4_zyx = simplify(R4_zyx,'steps',100);
disp('R4_ZYZ matrix generated...')
tic;matlabFunction(R4_zyx,'file','R4_zyx_0');disp('R4_zyx_0 is saved as Matlab function!');toc

R5_zyx = outer(R4_zyx,R1_zyx,0);
R5_zyx = permute(R5_zyx,[1,2,3,4,9,6,7,8,5,10]);
R5_zyx = permute(R5_zyx,[1,2,3,4,5,6,7,9,8,10]);
R5_zyx = permute(R5_zyx,[1,2,3,4,5,7,6,8,9,10]);
R5_zyx = reshape(R5_zyx,[243,243]);
R5_zyx = simplify(R5_zyx,'steps',100);
disp('R5_ZYZ matrix generated...')
tic;matlabFunction(R5_zyx,'file','R5_zyx_0');disp('R5_zyx_0 is saved as Matlab function!');toc

disp('All symbolic matrixies generated...')
toc


%% Translate and svae all symbolic matrixes of RRRRR into matlab functions






% disp('Function eneration finished!')

%% Avg over Phi, Psi, Phi&Psi, Phi&Psi&Theta 
% 
% tic
% R1_ZYZ_1   = int(R1_zyx               ,Phi  ,0,2*pi)./(2*pi);
% R1_ZXZ_2   = int(R1_ZXZ               ,Psi  ,0,2*pi)./(2*pi);
% R1_ZXZ_12  = int(R1_ZXZ_1             ,Psi  ,0,2*pi)./(2*pi);
% R1_ZXZ_123 = int(R1_ZXZ_12.*sin(Theta),Theta,0,  pi)./2;
% disp('integration of R1_ZYZ_1 is done!')
% toc
% 
% tic
% R2_ZYZ_1   = int(R2_zyx               ,Phi  ,0,2*pi)./(2*pi);
% R2_ZXZ_2   = int(R2_ZXZ               ,Psi  ,0,2*pi)./(2*pi);
% R2_ZXZ_12  = int(R2_ZXZ_1             ,Psi  ,0,2*pi)./(2*pi);
% R2_ZXZ_123 = int(R2_ZXZ_12.*sin(Theta),Theta,0,  pi)./2;
% disp('integration of R2_ZYZ_1 is done!')
% toc
% 
% tic
% R3_ZYZ_1   = int(R3_zyx               ,Phi  ,0,2*pi)./(2*pi);
% R3_ZXZ_2   = int(R3_ZXZ               ,Psi  ,0,2*pi)./(2*pi);
% R3_ZXZ_12  = int(R3_ZXZ_1             ,Psi  ,0,2*pi)./(2*pi);
% R3_ZXZ_123 = int(R3_ZXZ_12.*sin(Theta),Theta,0,  pi)./2;
% disp('integration of R3_ZYZ_1 is done!')
% toc
% 
% tic
% R4_ZYZ_1   = int(R4_zyx               ,Phi  ,0,2*pi)./(2*pi);
% R4_ZXZ_2   = int(R4_ZXZ               ,Psi  ,0,2*pi)./(2*pi);
% R4_ZXZ_12  = int(R4_ZXZ_1             ,Psi  ,0,2*pi)./(2*pi);
% R4_ZXZ_123 = int(R4_ZXZ_12.*sin(Theta),Theta,0,  pi)./2;
% disp('integration of R4_ZYZ_1 is done!')
% toc
% 
% tic
% R5_ZYZ_1   = int(R5_zyx               ,Phi  ,0,2*pi)./(2*pi);
% R5_ZXZ_2   = int(R5_ZXZ               ,Psi  ,0,2*pi)./(2*pi);
% R5_ZXZ_12  = int(R5_ZXZ_1             ,Psi  ,0,2*pi)./(2*pi);
% R5_ZXZ_123 = int(R5_ZXZ_12.*sin(Theta),Theta,0,  pi)./2;
% disp('integration of R4_ZYZ_1 is done!')
% toc

%% Translate and svae all symbolic matrixes of RRRRR into matlab functions

% tic;matlabFunction(R1_ZYZ_1  ,'file','R1_ZYZ_1'  );disp('R1_ZYZ_1   is saved as Matlab function!');toc
% tic;matlabFunction(R1_ZYZ_2  ,'file','R1_ZYZ_2'  );disp('R1_ZYZ_2   is saved as Matlab function!');toc
% tic;matlabFunction(R1_ZYZ_12 ,'file','R1_ZYZ_12' );disp('R1_ZYZ_12  is saved as Matlab function!');toc
% tic;matlabFunction(R1_ZYZ_123,'file','R1_ZYZ_123');disp('R1_ZYZ_123 is saved as Matlab function!');toc
% 
% tic;matlabFunction(R2_ZYZ_1  ,'file','R2_ZYZ_1'  );disp('R2_ZYZ_1   is saved as Matlab function!');toc
% tic;matlabFunction(R2_ZYZ_2  ,'file','R2_ZYZ_2'  );disp('R2_ZYZ_2   is saved as Matlab function!');toc
% tic;matlabFunction(R2_ZYZ_12 ,'file','R2_ZYZ_12' );disp('R2_ZYZ_12  is saved as Matlab function!');toc
% tic;matlabFunction(R2_ZYZ_123,'file','R2_ZYZ_123');disp('R2_ZYZ_123 is saved as Matlab function!');toc
% 
% tic;matlabFunction(R3_ZYZ_1  ,'file','R3_ZYZ_1'  );disp('R3_ZYZ_1   is saved as Matlab function!');toc
% tic;matlabFunction(R3_ZYZ_2  ,'file','R3_ZYZ_2'  );disp('R3_ZYZ_2   is saved as Matlab function!');toc
% tic;matlabFunction(R3_ZYZ_12 ,'file','R3_ZYZ_12' );disp('R3_ZYZ_12  is saved as Matlab function!');toc
% tic;matlabFunction(R3_ZYZ_123,'file','R3_ZYZ_123');disp('R3_ZYZ_123 is saved as Matlab function!');toc
% 
% tic;matlabFunction(R4_ZYZ_1  ,'file','R4_ZYZ_1'  );disp('R4_ZYZ_1   is saved as Matlab function!');toc
% tic;matlabFunction(R4_ZYZ_2  ,'file','R4_ZYZ_2'  );disp('R4_ZYZ_2   is saved as Matlab function!');toc
% tic;matlabFunction(R4_ZYZ_12 ,'file','R4_ZYZ_12' );disp('R4_ZYZ_12  is saved as Matlab function!');toc
% tic;matlabFunction(R4_ZYZ_123,'file','R4_ZYZ_123');disp('R4_ZYZ_123 is saved as Matlab function!');toc
% 
% tic;matlabFunction(R5_ZYZ_1  ,'file','R5_ZYZ_1'  );disp('R5_ZYZ_1   is saved as Matlab function!');toc
% tic;matlabFunction(R5_ZYZ_2  ,'file','R5_ZYZ_2'  );disp('R5_ZYZ_2   is saved as Matlab function!');toc
% tic;matlabFunction(R5_ZYZ_12 ,'file','R5_ZYZ_12' );disp('R5_ZYZ_12  is saved as Matlab function!');toc
% tic;matlabFunction(R5_ZYZ_123,'file','R5_ZYZ_123');disp('R5_ZYZ_123 is saved as Matlab function!');toc


