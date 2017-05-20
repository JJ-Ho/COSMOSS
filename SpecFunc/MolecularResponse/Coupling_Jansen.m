function Beta = Coupling_Jansen(Structure)
%% debug
% S = Data_COSMOSS.Structure;

%% Prep Variables
load('Jansen_Map.mat')

%% map phi,psi to NN beta
% Calculate diherdrals
[Phi,Psi] = PeptideDihedral(Structure);

% Convert phi,psi to map index
Phi_Ind = round(Phi + 180);
Psi_Ind = round(Psi + 180);

% if index = 0, move it to 1
Phi_Ind(eq(Phi_Ind,0)) = 1;
Psi_Ind(eq(Psi_Ind,0)) = 1;

% temporary assign nan inde to 1, will delete them later
NaN_Ind = isnan(Phi_Ind);
Phi_Ind(NaN_Ind) = 1;
Psi_Ind(NaN_Ind) = 1;

Index = sub2ind(size(Jansen_Map),Phi_Ind,Psi_Ind);
Beta_NN = Jansen_Map(Index);
Beta_NN(NaN_Ind,:) = 0; % zeros the non-NN

% create beta matrix
Beta = diag(Beta_NN,1) + diag(Beta_NN,-1);
