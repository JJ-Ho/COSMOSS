function [Beta,DistM] = Coupling_TDC(S)
%% Reassign variable names from StrucInfo
Nmodes    = S.Nmodes;
LocCenter = S.LocCenter;
LocMu     = S.LocMu;
LocFreq   = S.LocFreq;

% Coupling terms in one exciton Hamiltonian
[StateIndex1,StateIndex2] = ndgrid(1:Nmodes); % SI1 => I; SI2 => J

ModeDistVec = LocCenter(StateIndex1(:),:) - LocCenter(StateIndex2(:),:); % aka R
ModeDist    = sqrt(sum(abs(ModeDistVec).^2,2)); % Got to have abs, or otherwise will have imaginary part

muIdotJ = dot(LocMu(StateIndex1(:),:),LocMu(StateIndex2(:),:),2);
RdotmuI = dot(ModeDistVec,LocMu(StateIndex1(:),:),2);
RdotmuJ = dot(ModeDistVec,LocMu(StateIndex2(:),:),2);

%% Unit conversion

%-BetaPreFactor = 5034*((4.1058/sqrt(1600))*3.144)^2; %From Jenny's
%-mathematica code "TransitionDipoleCoupling.m" The result is in cm-1 unit
%-The input TDV in this case is an unit vector 
% BetaPreFactor = 5034*((4.1058/sqrt(1600))*3.144)^2;

%-From Chem. Phys. Lett. v106, #6, p613. The input TDV is in unit of  
%- sqrt(KM/mole) from G09
BetaPreFactor = 84862/LocFreq(1)*(0.23344)^2; % can be change to avg of LocFreq

%% Construct Beta
Beta = BetaPreFactor*(muIdotJ./(ModeDist.^3) - 3*RdotmuJ.*RdotmuI./(ModeDist).^5); % Beta in NumLocMode^2 x 3 size
Beta = reshape(Beta,Nmodes,Nmodes);

Beta(isnan(Beta)) = 0; % Get ride of diagnal
Beta = (Beta+Beta')./2; % remove numerical error, make the matrix exactly symmetric

%% Construct distance matrix
DistM = reshape(sqrt(sum(ModeDistVec.^2,2)),Nmodes,Nmodes);