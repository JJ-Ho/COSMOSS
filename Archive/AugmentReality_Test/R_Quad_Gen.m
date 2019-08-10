%% Generation of non-averaged R1 to R5
clear all
%%
% Elapsed time is 1.371786 seconds.
tic
syms a b c d

% Original version for permutation
R1_Quad=[ a^2+b^2-c^2-d^2,     2*(b*c-a*d),      2*(b*d+a*c);
              2*(b*c+a*d), a^2-b^2+c^2-d^2,      2*(c*d-a*b);
              2*(b*d-a*c),     2*(c*d+a*b), a^2-b^2-c^2+d^2];

disp('R1_Quad matrix generated...')
tic;matlabFunction(R1_Quad,'file','R1_Quad');disp('R1_Quad is saved as Matlab function!');toc

R2_Quad = outer(R1_Quad,R1_Quad,0);
R2_Quad = permute(R2_Quad,[1,3,2,4]);
R2_Quad = reshape(R2_Quad,[9,9]);
R2_Quad = simplify(R2_Quad,'steps',100);
disp('R2_Quad matrix generated...')
tic;matlabFunction(R2_Quad,'file','R2_Quad');disp('R2_Quad is saved as Matlab function!');toc


R3_Quad = outer(R2_Quad,R1_Quad,0);
R3_Quad = permute(R3_Quad,[1,2,5,4,3,6]);
R3_Quad = permute(R3_Quad,[1,2,3,5,4,6]);
R3_Quad = reshape(R3_Quad,[27,27]);
R3_Quad = simplify(R3_Quad,'steps',100);
disp('R3_Quad matrix generated...')
tic;matlabFunction(R3_Quad,'file','R3_Quad');disp('R3_Quad is saved as Matlab function!');toc

R4_Quad = outer(R3_Quad,R1_Quad,0);
R4_Quad = permute(R4_Quad,[1,2,3,7,5,6,4,8]);
R4_Quad = permute(R4_Quad,[1,2,3,4,5,7,6,8]);
R4_Quad = permute(R4_Quad,[1,2,3,4,6,5,7,8]);
R4_Quad = reshape(R4_Quad,[81,81]);
R4_Quad = simplify(R4_Quad,'steps',100);
disp('R4_Quad matrix generated...')
tic;matlabFunction(R4_Quad,'file','R4_Quad');disp('R4_Quad is saved as Matlab function!');toc

R5_Quad = outer(R4_Quad,R1_Quad,0);
R5_Quad = permute(R5_Quad,[1,2,3,4,9,6,7,8,5,10]);
R5_Quad = permute(R5_Quad,[1,2,3,4,5,6,7,9,8,10]);
R5_Quad = permute(R5_Quad,[1,2,3,4,5,7,6,8,9,10]);
R5_Quad = reshape(R5_Quad,[243,243]);
R5_Quad = simplify(R5_Quad,'steps',100);
disp('R5_Quad matrix generated...')
tic;matlabFunction(R5_Quad,'file','R5_Quad');disp('R5_Quad is saved as Matlab function!');toc

disp('All symbolic matrixies generated...')
toc
