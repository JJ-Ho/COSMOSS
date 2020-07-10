% Run COSMOSS once and leave it open before using this script. In this way
% the necessary paths are registered before executing this code.
%% Construct NCOs
% parameter setups 
vibFreq     = 1650; % vibrational local mode frequencies
couplingVC  = 20; % coupling constant between the vibrational mode and the cavity mode
cavScanFreq = -100:1:100; % scanning range (cm-1)
minCavTrans = 1450; % minimum cavity transmission (cm-1)
FSR         = 200;  % Freq spectral range (cm-1)
N_CavTrans  = 3;    % number of cavity transmissions

% NCO model structural parameters (aka the Table input of the NCO)
N_V = length(vibFreq);    % number of vibrational modes
StructureInputs = zeros(N_V,8);
StructureInputs(:,1)  =  ones(N_V,1); % N_Modes_toggle
StructureInputs(:,2)  = zeros(N_V,1); % Phi_monomer
StructureInputs(:,3)  = zeros(N_V,1); % Psi_monomer
StructureInputs(:,4)  = zeros(N_V,1); % Theta_monomer
StructureInputs(:,5)  = zeros(N_V,1); % X_monomer
StructureInputs(:,6)  = zeros(N_V,1); % Y_monomer
StructureInputs(:,7)  = zeros(N_V,1); % Z_monomer
StructureInputs(:,8)  =  ones(N_V,1); % ScalingFactor
StructureInputs(:,9)  =      vibFreq; % vibrational local mode frequencies
StructureInputs(:,10) = zeros(N_V,1); % vibrational mode anharmonicity

% non-default GUI Inputs for NCO
G1.StructureInputs = StructureInputs; 
G1.CouplingType    = 'Zero';

% non-default GUI Inputs for Cavity
G2.NCModes = N_CavTrans;
G2.FSR     = FSR;
G2.NLFreq  = minCavTrans;

% non-default GUI Inputs for Comb2
GC.CouplingType = 'Constant';
GC.Beta_NN      = couplingVC;

% non-default GUI Inputs for FTIR
FreqRangeMin       = (min(cavScanFreq) + minCavTrans -100);
FreqRangeMax       = (max(cavScanFreq) + minCavTrans + (N_CavTrans-1)*FSR +100);
GFTIR.FreqRange    = FreqRangeMin:FreqRangeMax;
GFTIR.LineWidth_1D = 10;

%% Generate the dispersion plot
Z = zeros(length(cavScanFreq),length(GFTIR.FreqRange));

for i = 1:length(cavScanFreq)

    G2.NLFreq  = minCavTrans + cavScanFreq(i); % scan the Cavity transmissions
    
    SD_NCO        = ConstructNCO(G1);
    SD_Cavity     = ConstructCavity(G2);
    SD_Cavity_NCO = Comb2(SD_NCO,SD_Cavity,GC);
    
    OneD = FTIR_Main(SD_Cavity_NCO,GFTIR);
    [~,Output] = Plot1D('none',OneD,GFTIR);
    Z(i,:) = (Output.Y)';    
end

%% make figure
hF = figure;
hAx = axes('Parent',hF);

X = GFTIR.FreqRange;
Y = cavScanFreq;

imagesc(hAx,X,Y,Z)
hAx.XLabel.String = 'Transmission Spectra (cm^{-1})';
hAx.YLabel.String = 'Cavity detuning (cm^{-1})';
