function [Beta,dF] = Coupling_Jansen(Structure)
%% debug
% S = Data_COSMOSS.Structure;

%% Prep Variables
load('Jansen_Map.mat')

%% Calculating Dihedral angles
% Calculate diherdrals
%        O(i)                     O(i+1)
%        ||                       ||
% (N)... C(i) -> N(i) -> Ca(i) -> C(i+1) -> N(i+1) - ... (C)
%                   Phi(i)    Psi(i)

[Phi,Psi] = PeptideDihedral(Structure);

% Convert phi,psi to map index, [-180,180] => [0,360]
Phi_Ind = round(Phi + 180);
Psi_Ind = round(Psi + 180);

% if index = 0, move it to 1, [0,360] => [1,360]
Phi_Ind(eq(Phi_Ind,0)) = 1;
Psi_Ind(eq(Psi_Ind,0)) = 1;

% temporary assign nan inde to 1, will delete them later
NaN_Ind = isnan(Phi_Ind);
Phi_Ind(NaN_Ind) = 1;
Psi_Ind(NaN_Ind) = 1;

%% Nearest Neighbor Couping
% Index map
% Couping(i,i+1) <= (psi(i),phi(i))

NNC_Ind = sub2ind(size(Jansen_Map),Psi_Ind,Phi_Ind);
Beta_NN = Jansen_Map(NNC_Ind);
Beta_NN(NaN_Ind,:) = 0; % zeros the non-NN

% create beta matrix
Beta = diag(Beta_NN,1) + diag(Beta_NN,-1);

%% Local Mode Frequency Shift
% Index mapping
% # Amide         1,2,3,4,5,6,...7,8,9,10,...
% # N(Phi,Psi)    x,1,2,3,4,5,...x,7,8, 9,...
% # C(Phi,Psi)    1,2,3,4,5,x,...7,8,9, x,...

dF_N_Ind = [1;NNC_Ind];
dF_N = Jansen_wn(dF_N_Ind);
dF_N([true;NaN_Ind]) = 0;

dF_C_Ind = [NNC_Ind;1];
dF_C = Jansen_wc(dF_C_Ind);
dF_C([NaN_Ind;true]) = 0;

dF = dF_N + dF_C;

