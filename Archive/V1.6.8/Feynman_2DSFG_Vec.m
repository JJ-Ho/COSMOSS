function [Freq,Beta,Index] = Feynman_2DSFG_Vec(F1,F2,A_Ex_01,A_Ex_12,M_Ex_01,M_Ex_12)
% 
% This function generate Feynman pathway of 2DSFG with given polarization.
% 
% Since "kron" function use reshape alot, the overall speed of kron product
% version is slower than for loop version
% 
% Note: Alpha_Ex = [size(H) x 9], index: [xx xy xz yx yy yz zx zy zz]
% 
% Todo: code acceleration, delete small signals.
% 
% ------- Version log -----------------------------------------------------
% 
% Ver. 1.0  140126  modified from Feynman_2DIR_kron; vecterized version
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2014-2016

%% debug
% clear all
% F1 = (1:3);
% F2 = (4:7);
% 
% M_Ex_01 = rand(length(F1),3);
% M_Ex_12 = rand(length(F1),length(F2),3);
% A_Ex_01 = rand(length(F1),9);
% A_Ex_12 = rand(length(F1),length(F2),9);

%% Prepare index
% Find number of modes in 1ex and 2ex
N1 = length(F1);
N2 = length(F2);

% R1 R2 NR1 NR2 index expansion,size N1^2 
[Ib,Ia] = ndgrid(1:N1,1:N1);
% XYZ index expansion, size 3^4
[Jd,Jc,Jb,Ja] = ndgrid(1:3,1:3,1:3,1:9);
% R3 NR3 index expansion, size: N2*N1^2
[Kx,Kb,Ka] = ndgrid(1:N2,1:N1,1:N1);

%[length(Ib(:)),length(Kx(:))] %% Debug

%% Reduce mode number base on the transition intensity
Cut_Off_01 = 0;
Cut_Off_12 = 2E-2;

% 01
Norm_M_01  = sqrt(sum(M_Ex_01.^2,2));

Max_Norm_M_01 = max(Norm_M_01);
Ind_Norm_M_01 = Norm_M_01 > Cut_Off_01 * Max_Norm_M_01;

[Mask_Ib,Mask_Ia] = ndgrid(Ind_Norm_M_01,Ind_Norm_M_01);
Ib = Mask_Ib.*Ib;
Ia = Mask_Ia.*Ia;

Ib = Ib(Ib>0);
Ia = Ia(Ia>0);

% 12
Ind_Norm_M_01 = Norm_M_01 > Cut_Off_12 * Max_Norm_M_01;

Norm_M_12  = sqrt(sum(M_Ex_12.^2,3));
Norm_A_12  = sqrt(sum(A_Ex_12.^2,3));

Max_Norm_M_12 = max(Norm_M_12(:));
Max_Norm_A_12 = max(Norm_A_12(:));

Ind_Norm_M_12 = Norm_M_12 > Cut_Off_12 * Max_Norm_M_12;
Ind_Norm_A_12 = Norm_A_12 > Cut_Off_12 * Max_Norm_A_12;
Ind_Norm_AM_12 = and(Ind_Norm_M_12,Ind_Norm_A_12);

Mask_Kx = repmat(Ind_Norm_AM_12', 1,1,N1);
Mask_Kb = repmat(Ind_Norm_M_01',N2,1,N1);
Mask_Ka = permute(Mask_Kb,[1,3,2]);

Kx = Mask_Kx.*Mask_Kb.*Mask_Ka.*Kx;
Kb = Mask_Kx.*Mask_Kb.*Mask_Ka.*Kb;
Ka = Mask_Kx.*Mask_Kb.*Mask_Ka.*Ka;

Kx(Kx==0) = NaN;
Kb(Kb==0) = NaN;
Ka(Ka==0) = NaN;

% Indice of 2D version Mu_Ex, 
% Note [N,N*(N+1)/2] is the first two Dimension of Mu_Ex
ind_ax = sub2ind([N1,N2],Ka(:),Kx(:));
ind_bx = sub2ind([N1,N2],Kb(:),Kx(:));

ind_ax(isnan(ind_ax)) = [];
ind_bx(isnan(ind_bx)) = [];

ind_xa = ind_ax;
ind_xb = ind_bx;

Kx(isnan(Kx)) = [];
Kb(isnan(Kb)) = [];
Ka(isnan(Ka)) = [];

%[length(Ib(:)),length(Kx(:))] %% Debug

%% Get 2DSFG response from index version kronec product

% R1 R2 NR1 NR2
M_0a = M_Ex_01(Ia(:),:);
M_0b = M_Ex_01(Ib(:),:);
M_a0 = M_0a;
M_b0 = M_0b;

A_a0 = A_Ex_01(Ia(:),:);
A_b0 = A_Ex_01(Ib(:),:);

% Response (GB SE)
R1  = A_b0(:,Ja(:)).*M_0b(:,Jb(:)).*M_a0(:,Jc(:)).*M_0a(:,Jd(:));
R2  = A_b0(:,Ja(:)).*M_a0(:,Jb(:)).*M_0b(:,Jc(:)).*M_0a(:,Jd(:));
NR1 = A_b0(:,Ja(:)).*M_0b(:,Jb(:)).*M_a0(:,Jc(:)).*M_0a(:,Jd(:));
NR2 = A_a0(:,Ja(:)).*M_b0(:,Jb(:)).*M_0b(:,Jc(:)).*M_0a(:,Jd(:));

% R3 NR3
N_0a = M_Ex_01(Ka(:),:);
N_0b = M_Ex_01(Kb(:),:);

% Merge the first two index of Mu_Ex into one => 2D verion Mu_Ex
TD_M_Ex = reshape(M_Ex_12,[],3);
TD_A_Ex = reshape(A_Ex_12,[],9);

% get TDV using linear index of the first two indice of Mu_Ex 
% and ":" to extract the whole vector components 
N_ax = TD_M_Ex(ind_ax,:);
N_bx = TD_M_Ex(ind_bx,:);
N_xa = TD_A_Ex(ind_xa,:);
N_xb = TD_A_Ex(ind_xb,:);

% Response (EA)
R3  = N_xa(:,Ja(:)).*N_bx(:,Jb(:)).*N_0b(:,Jc(:)).*N_0a(:,Jd(:));
NR3 = N_xb(:,Ja(:)).*N_ax(:,Jb(:)).*N_0b(:,Jc(:)).*N_0a(:,Jd(:));

%% Prep for export mode index
I_R1  = [Ib(:),Ib(:),Ia(:),Ia(:)]; % GB
I_R2  = [Ib(:),Ia(:),Ib(:),Ia(:)]; % SE
I_R3  = [Ka(:),Kx(:),Kb(:),Ka(:)]; % EA

I_NR1 = [Ib(:),Ib(:),Ia(:),Ia(:)]; % GB
I_NR2 = [Ia(:),Ib(:),Ib(:),Ia(:)]; % SE
I_NR3 = [Kb(:),Kx(:),Kb(:),Ka(:)]; % EA

%% Generate List of interaction Frequencies
Ea_12 = F1(Ia(:));
Eb_12 = F1(Ib(:));

Freq_R1  = [-Ea_12, zeros(size(Ea_12)), Eb_12];
Freq_R2  = [-Ea_12, Eb_12 - Ea_12     , Eb_12];
Freq_NR1 = [ Ea_12, zeros(size(Ea_12)), Eb_12];
Freq_NR2 = [ Ea_12, Ea_12 - Eb_12     , Ea_12]; % This is a interesting term! only contribute to diagonal!

Ea_3 = F1(Ka(:));
Eb_3 = F1(Kb(:));
Ex_3 = F2(Kx(:));

Freq_R3  = [-Ea_3, Eb_3 - Ea_3, Ex_3 - Ea_3];
Freq_NR3 = [ Ea_3, Ea_3 - Eb_3, Ex_3 - Eb_3];

%whos %% debug

%% Output
Beta.R1  = R1';  % [243 x N^2]
Beta.R2  = R2';  % [243 x N^2]
Beta.R3  = R3';  % [243 x N^3*(N+1)/2]
Beta.NR1 = NR1'; % [243 x N^2]
Beta.NR2 = NR2'; % [243 x N^2]
Beta.NR3 = NR3'; % [243 x N^3*(N+1)/2]

Freq.R1  = Freq_R1;
Freq.R2  = Freq_R2;
Freq.R3  = Freq_R3;
Freq.NR1 = Freq_NR1;
Freq.NR2 = Freq_NR2;
Freq.NR3 = Freq_NR3;

Index.R1  = I_R1;
Index.R2  = I_R2;
Index.R3  = I_R3;
Index.NR1 = I_NR1;
Index.NR2 = I_NR2;
Index.NR3 = I_NR3;