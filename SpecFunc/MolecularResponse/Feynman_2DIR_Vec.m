function [Freq,Beta] = Feynman_2DIR_Vec(N,F1,F2,M_Ex_01,M_Ex_12)
% 
% This function generate Feynman pathway of 2DSFG with given polarization.
% 
% Since "kron" function use reshape alot, the overall speed of kron product
% version is slower than for loop version
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
% N = 3;
% 
% Sort_Ex_Freq = rand((N+1)*(N+2)/2,1)*1000;
% 
% Mask = ones((N+1)*(N+2)/2);
% Mask = Mask .* ~blkdiag(1,ones(N),ones(N*(N+1)/2));
% Mask(1,N+2:(N+1)*(N+2)/2) = 0;
% Mask(N+2:(N+1)*(N+2)/2,1) = 0;
% 
% Mu_Ex = rand((N+1)*(N+2)/2,(N+1)*(N+2)/2,3);
% Mu_Ex = Mu_Ex.* repmat(Mask,[1,1,3]);

%% Prepare index

% R1 R2 NR1 NR2 index expansion,size N^2 
[Ib,Ia] = ndgrid(1:N,1:N);
% XYZ index expansion, size 3^4
[Jd,Jc,Jb,Ja] = ndgrid(1:3,1:3,1:3,1:3);
% R3 NR3 index expansion, size: N^3*(N+1)/2
[Kx,Kb,Ka] = ndgrid(1:N*(N+1)/2,1:N,1:N);

%% Get 2DIR response from index version kronec product

% R1 R2 NR1 NR2
M_0a = M_Ex_01(Ia(:),:);
M_0b = M_Ex_01(Ib(:),:);
M_a0 = M_0a;
M_b0 = M_0b;

% Response (GB SE)
R1  = M_b0(:,Ja(:)).*M_0b(:,Jb(:)).*M_a0(:,Jc(:)).*M_0a(:,Jd(:));
R2  = M_b0(:,Ja(:)).*M_a0(:,Jb(:)).*M_0b(:,Jc(:)).*M_0a(:,Jd(:));
NR1 = M_b0(:,Ja(:)).*M_0b(:,Jb(:)).*M_a0(:,Jc(:)).*M_0a(:,Jd(:));
NR2 = M_a0(:,Ja(:)).*M_b0(:,Jb(:)).*M_0b(:,Jc(:)).*M_0a(:,Jd(:));


% R3 NR3
N_0a = M_Ex_01(Ka(:),:);
N_0b = M_Ex_01(Kb(:),:);

% Merge the first two index of Mu_Ex into one => 2D verion Mu_Ex
TD_M_Ex = reshape(M_Ex_12,[],3);

% Indice of 2D version Mu_Ex, 
% Note [N,N*(N+1)/2] is the first two Dimension of Mu_Ex
ind_ax = sub2ind([N,N*(N+1)/2],Ka(:),Kx(:));
ind_bx = sub2ind([N,N*(N+1)/2],Kb(:),Kx(:));
ind_xa = ind_ax;
ind_xb = ind_bx;

% get TDV using linear index of the first two indice of Mu_Ex 
% and ":" to extract the whole vector components 
N_ax = TD_M_Ex(ind_ax,:);
N_bx = TD_M_Ex(ind_bx,:);
N_xa = TD_M_Ex(ind_xa,:);
N_xb = TD_M_Ex(ind_xb,:);

% Response (EA)
R3  = N_xa(:,Ja(:)).*N_bx(:,Jb(:)).*N_0b(:,Jc(:)).*N_0a(:,Jd(:));
NR3 = N_xb(:,Ja(:)).*N_ax(:,Jb(:)).*N_0b(:,Jc(:)).*N_0a(:,Jd(:));

%% Generate List of interaction Frequencies
Ea_12 = F1(Ia(:))';
Eb_12 = F1(Ib(:))';

Freq_R1  = [-Ea_12 ; zeros(1,N^2)  ; Eb_12]';
Freq_R2  = [-Ea_12 ; Eb_12 - Ea_12 ; Eb_12]';
Freq_NR1 = [ Ea_12 ; zeros(1,N^2)  ; Eb_12]';
Freq_NR2 = [ Ea_12 ; Ea_12 - Eb_12 ; Ea_12]'; % This is a interesting term! only contribute to diagonal!

Ea_3 = F1(Ka(:))';
Eb_3 = F1(Kb(:))';
Ex_3 = F2(Kx(:))';

Freq_R3  = [-Ea_3 ; Eb_3 - Ea_3 ; Ex_3 - Ea_3]';
Freq_NR3 = [ Ea_3 ; Ea_3 - Eb_3 ; Ex_3 - Eb_3]';

%% Output
Beta.R1  = R1';  % [81 x N^2]
Beta.R2  = R2';  % [81 x N^2]
Beta.R3  = R3';  % [81 x N^3*(N+1)/2]
Beta.NR1 = NR1'; % [81 x N^2]
Beta.NR2 = NR2'; % [81 x N^2]
Beta.NR3 = NR3'; % [81 x N^3*(N+1)/2]

% Sparse matrix version
% Response.R1  = sparse(R1)';  % [81 x N^2]
% Response.R2  = sparse(R2)';  % [81 x N^2]
% Response.R3  = sparse(R3)';  % [81 x N^3*(N+1)/2]
% Response.NR1 = sparse(NR1)'; % [81 x N^2]
% Response.NR2 = sparse(NR2)'; % [81 x N^2]
% Response.NR3 = sparse(NR3)'; % [81 x N^3*(N+1)/2]

Freq.R1  = Freq_R1;
Freq.R2  = Freq_R2;
Freq.R3  = Freq_R3;
Freq.NR1 = Freq_NR1;
Freq.NR2 = Freq_NR2;
Freq.NR3 = Freq_NR3;

