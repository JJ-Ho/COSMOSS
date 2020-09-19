function Output= H_handler(SData,Main_GUI_Inputs,ExMode)
% this function handle the Hamiltonian generation from the StrutureData with
% the following tasks:
%   1. Diagonal disorder
%   2. Off-Diagonal disorder
%   3. Generate Two-Exciton part if needed
%   4. Diagnalization

%% Inputs parser
GUI_Inputs_C      = fieldnames(Main_GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(Main_GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = true;

defaultSampling     = 0;
defaultDD_FWHM      = 0;
defaultODD_FWHM     = 0;
defaultP_FlucCorr   = 0;

addOptional(INPUT,'Sampling'    ,defaultSampling);
addOptional(INPUT,'DD_FWHM'     ,defaultDD_FWHM);
addOptional(INPUT,'ODD_FWHM'    ,defaultODD_FWHM);
addOptional(INPUT,'P_FlucCorr'  ,defaultP_FlucCorr);

parse(INPUT,GUI_Inputs_C{:});

Sampling     = INPUT.Results.Sampling;
DD_FWHM      = INPUT.Results.DD_FWHM;
ODD_FWHM     = INPUT.Results.ODD_FWHM;
P_FlucCorr   = INPUT.Results.P_FlucCorr;

%% reassign variable names
N = SData.Nmodes;
SData_tobedrop = SD_Copy(SData); % to avoide later DD/ODD dBeta modifications being accumulated

%% check if apply random sampling
if Sampling
    DD_std  =  DD_FWHM/(2*sqrt(2*log(2)));
    ODD_std = ODD_FWHM/(2*sqrt(2*log(2)));
else
    DD_std  = 0;
    ODD_std = 0;
end

%% Diagonal/Off-Diaginal disorder if any
P_FlucCorr = P_FlucCorr/100; % turn percentage to number within 0~1

% diagonal disorder
Correlation_Dice = rand;
if Correlation_Dice < P_FlucCorr
    dF_DD = DD_std.*(randn.*ones(N,1));
else 
    dF_DD = DD_std.*randn(N,1); 
end

% Off diagonal disorder
dBeta   = ODD_std*randn(N);
dBeta   = (dBeta + dBeta')./2; % symetrize
dBeta(logical(eye(N))) = zeros(N,1); % Get ride of diagnal

% Update the local mode frequencies and the couplings
SData_tobedrop.LocFreq = SData_tobedrop.LocFreq + dF_DD;
SData_tobedrop.Beta    = SData_tobedrop.Beta + dBeta;
OneExH                 = SData_tobedrop.OneExH;

%% Generate Two Exciton block diag part of full Hamiltonianif needed (TwoExOvertoneH & TwoExCombinationH)
TEDIndexBegin = [];
TEDIndexEnd   = [];
TwoExPart     = [];

if strcmp(ExMode,'TwoEx')   
    TwoExH = SData_tobedrop.TwoExH;
    TEDIndexBegin = TwoExH.TEDIndexBegin;
    TEDIndexEnd   = TwoExH.TEDIndexEnd;
    TwoExPart     = TwoExH.TwoExPart;    
end



%% Diagonalize the full hamiltonian
Sort_Ex_F1 = [];
Sort_Ex_V1 = [];
Sort_Ex_F2 = [];
Sort_Ex_V2 = [];

if strcmp(ExMode,'TwoEx')
    FullH = blkdiag(OneExH,TwoExPart);
else
    FullH = OneExH;
end 


% Diagonalization
if ~isobject(SData_tobedrop.LocFreq)
    % note: the eiganvector V_Full(:,i) has been already normalized.
    [V_Full,D_Full] = eig(FullH);
    Ex_Freq = diag(D_Full);

    % sort eiganvalue form small to big and reorder the eiganvectors
    [Sort_Ex_Freq,Indx] = sort(Ex_Freq);
     Sort_Ex_V          = V_Full(:,Indx);

    Sort_Ex_F1 = Sort_Ex_Freq(2:N+1);
    Sort_Ex_V1 = Sort_Ex_V(2:N+1,2:N+1);

    if strcmp(ExMode,'TwoEx')
        
        Sort_Ex_F2 = Sort_Ex_Freq(N+2:end);
        Sort_Ex_V2 = Sort_Ex_V(N+2:end,N+2:end);
    end 
end 
%% export
Output.Sort_Ex_F1    = Sort_Ex_F1;
Output.Sort_Ex_V1    = Sort_Ex_V1;
Output.Sort_Ex_F2    = Sort_Ex_F2;
Output.Sort_Ex_V2    = Sort_Ex_V2;
Output.TEDIndexBegin = TEDIndexBegin;
Output.TEDIndexEnd   = TEDIndexEnd;
Output.TwoExPart     = TwoExPart;
Output.OneExH        = OneExH;
Output.dLocFreq      = SData_tobedrop.LocFreq;
Output.dBeta         = SData_tobedrop.Beta;
