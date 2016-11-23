function Beta = Coupling_TDC_PBC(S)
% Apply periodic boundary condition on the edge with only considering
% nearest neighbors
%% Reassign variable names from StrucInfo
Num_Modes = S.Num_Modes;
Center    = S.center;
Mu        = S.mu;
Freq      = S.freq;
% for PBC
V1     = S.Vec_1;
V2     = S.Vec_2;
NV1    = S.N_Vec1;
NV2    = S.N_Vec2;

%% Prep for PBC index pairs
%       ^
%   V2 / E3
%     ------
%  E2/     / E4
%   /     /
%   ------> V1
%     E1

E1_Ind = sub2ind([NV1,NV2],           1:NV1, ones(1,NV1)     );
E2_Ind = sub2ind([NV1,NV2],ones(1,NV2)     ,            1:NV2);
E3_Ind = sub2ind([NV1,NV2],           1:NV1, ones(1,NV1).*NV2);
E4_Ind = sub2ind([NV1,NV2],ones(1,NV2).*NV1,            1:NV2);

% Transform the mode index to row index of Center
E24 = sub2ind([Num_Modes,Num_Modes],E2_Ind,E4_Ind); % V1
E42 = sub2ind([Num_Modes,Num_Modes],E4_Ind,E2_Ind); % V1
E13 = sub2ind([Num_Modes,Num_Modes],E1_Ind,E3_Ind); % V2
E31 = sub2ind([Num_Modes,Num_Modes],E3_Ind,E1_Ind); % V2

%% Coupling terms in one exciton Hamiltonian
[StateIndex1,StateIndex2] = ndgrid(1:Num_Modes); % SI1 => I; SI2 => J

ModeDistVec = Center(StateIndex1(:),:) - Center(StateIndex2(:),:); % aka R

% substitute PBC
ModeDistVec(E24,:) = repmat(V1,NV2,1);
ModeDistVec(E42,:) = repmat(V1,NV2,1);
ModeDistVec(E13,:) = repmat(V2,NV1,1);
ModeDistVec(E31,:) = repmat(V2,NV1,1);

ModeDist    = sqrt(sum(abs(ModeDistVec).^2,2)); % Got to have abs, or otherwise will have imaginary part

muIdotJ = dot(Mu(StateIndex1(:),:),Mu(StateIndex2(:),:),2);
RdotmuI = dot(ModeDistVec,Mu(StateIndex1(:),:),2);
RdotmuJ = dot(ModeDistVec,Mu(StateIndex2(:),:),2);

%% Unit conversion

%-BetaPreFactor = 5034*((4.1058/sqrt(1600))*3.144)^2; %From Jenny's
%-mathematica code "TransitionDipoleCoupling.m" The result is in cm-1 unit
%-The input TDV in this case is an unit vector 
% BetaPreFactor = 5034*((4.1058/sqrt(1600))*3.144)^2;

%-From Chem. Phys. Lett. v106, #6, p613. The input TDV is in unit of  
%- sqrt(KM/mole) from G09
BetaPreFactor = 84862/Freq(1)*(0.23344)^2; 

%% Construct Beta
Beta = BetaPreFactor*(muIdotJ./(ModeDist.^3) - 3*RdotmuJ.*RdotmuI./(ModeDist).^5); % Beta in NumLocMode^2 x 3 size
Beta = reshape(Beta,Num_Modes,Num_Modes);

Beta(isnan(Beta)) = 0; % Get ride of diagnal