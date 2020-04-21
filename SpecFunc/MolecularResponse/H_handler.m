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
defaultP_FlucCorr   = 100;

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

%% check if apply random sampling
if Sampling
    DD_std  = DD_FWHM./(2*sqrt(2*log(2)));
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
    dF_DD = DD_std.*(randn(1,1).*ones(N,1));
else 
    dF_DD = DD_std.*randn(N,1); 
end

% Off diagonal disorder
dBeta   = ODD_std*randn(N);
dBeta   = (dBeta + dBeta')./2; % symetrize
dBeta(logical(eye(N))) = zeros(N,1); % Get ride of diagnal

% Update the local mode frequencies and the couplings
SData.LocFreq = SData.LocFreq + dF_DD;
SData.Beta    = SData.Beta + dBeta;
OneExH        = SData.OneExH;

%% Polariton 1Ex
%define cavity energy here for now
w_c = 1700;  %cavity energy in units of wavenumber
g = 10;      %cavity molecule coupling term

%Extend First Excited State Hamiltonian..
OneExH = blkdiag(OneExH,w_c);
OneExH(2:N+1,N+2)=g;
OneExH(N+2,2:N+1)=g;


%% Generate Two Exciton block diag part of full Hamiltonianif needed (TwoExOvertoneH & TwoExCombinationH)
TEDIndexBegin = [];
TEDIndexEnd   = [];
TwoExPart     = [];

if strcmp(ExMode,'TwoEx')   
    TwoExH = SData.TwoExH;
    TEDIndexBegin = TwoExH.TEDIndexBegin;
    TEDIndexEnd   = TwoExH.TEDIndexEnd;
    TwoExPart     = TwoExH.TwoExPart;    
    
    %Extend Second Excited State Hamiltonian...
    %energies of one cavity excitation + (one molecular excitation or cavity excitation)
    TwoExCavH= diag(diag(OneExH(2:N+2,2:N+2))+w_c);
    TwoExCavH(1:N,N+1)=g;
    TwoExCavH(N+1,1:N)=g;
    TwoExPart = blkdiag(TwoExPart,TwoExCavH);

    %Fill in cross terms of between cavity and 2nd excited state Hamiltonian
    for i=1:N
       TwoExPart(((N+1)*N/2)-((N-i+2)*(N-i+1)/2)+1:((N+1)*N/2)-((N-i+2)*(N-i+1)/2)+1+(N-i),...
                 ((N+1)*N/2+i):((N+1)*N/2+N))=diag(ones(N-i+1,1))*g;  

       TwoExPart(((N+1)*N/2+i):((N+1)*N/2+N),...
                 ((N+1)*N/2)-((N-i+2)*(N-i+1)/2)+1:((N+1)*N/2)-((N-i+2)*(N-i+1)/2)+1+(N-i))=diag(ones(N-i+1,1))*g;        

    end

    for i=1:N-1
        TwoExPart(((N+1)*N/2)-((N-i+1)*(N-i)/2)+1:((N+1)*N/2)-((N-i+1)*(N-i+0)/2)+N-i,(N+1)*N/2+i)=g*ones(N-i,1);
        TwoExPart((N+1)*N/2+i,((N+1)*N/2)-((N-i+1)*(N-i)/2)+1:((N+1)*N/2)-((N-i+1)*(N-i+0)/2)+N-i)=g*ones(N-i,1)';
    end
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
if ~isobject(SData.LocFreq)
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
