function Beta = Coupling_Cavity(S)
%% Reassign variable names from StrucInfo
Nmodes    = S.Nmodes;
LocMu     = S.LocMu;
CavityInd = S.Extra.CavityInd;

% substitue cavity's LocMu with unit vactor along X
CavityAxes = zeros(size(CavityInd,1),3);
CavityAxes(:,1) = 1;
LocMu(CavityInd,:) = CavityAxes;

% Coupling terms in one exciton Hamiltonian
[StateIndex1,StateIndex2] = ndgrid(1:Nmodes); % SI1 => I; SI2 => J
muIdotJ = dot(LocMu(StateIndex1(:),:),LocMu(StateIndex2(:),:),2);

%% Unit conversion

%-BetaPreFactor = 5034*((4.1058/sqrt(1600))*3.144)^2; %From Jenny's
%-mathematica code "TransitionDipoleCoupling.m" The result is in cm-1 unit
%-The input TDV in this case is an unit vector 
% BetaPreFactor = 5034*((4.1058/sqrt(1600))*3.144)^2;

%-From Chem. Phys. Lett. v106, #6, p613. The input TDV is in unit of  
%- sqrt(KM/mole) from G09
% BetaPreFactor = 84862/LocFreq(1)*(0.23344)^2; % can be change to avg of LocFreq
preFactors = 1;

 
%% Construct Beta
Beta = preFactors .* muIdotJ; % Need to figure out the preFactor in cm-1 unit
Beta = reshape(Beta,Nmodes,Nmodes);

Beta(isnan(Beta)) = 0; % Get ride of diagnal
Beta = (Beta+Beta')./2; % remove numerical error, make the matrix exactly symmetric

%% Remove non-cavity to molecule couping
MoleculeInd = 1:Nmodes;
MoleculeInd(CavityInd) = [];
Beta(MoleculeInd,MoleculeInd) = 0;