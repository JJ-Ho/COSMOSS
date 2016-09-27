% Holstein-like Hamiltonian 
% This hamiltonian consist of a dimer system compose of identical monomer. 
% This system only invole 1 electronic state 1 vibrational state on the 1st
% electronic exicted state for each monomer. 
% Ref: Kistler, K. A.; Pochas, C. M.; Yamagata, H. JPCB, 2011, 116, 77?86.

%% options
% molecule parameters
NV = 3; % NUmber of maximun vibrational quata
w0 =1400;
Lambda = sqrt(0.57);
J12 = 100;
D = -100;
w0_0 = 19000;

% Transition dipole in local mode
% mu1 = [1,1,0]./sqrt(2);
% mu2 = [1,-1,0]./sqrt(2);
mu1 = [1,0,0];
mu2 = [1,0,0];

% figure
PlotStick = 1;
LS = 'L';
LineWidth = 0.4*w0;
F_Min = 16000;
F_Max = 30000;

%% Construct dimer
% beta version: run the Model_TCO in stand alone mode then extract
% structural info from exported handles

%% Construct local mode basis

NS = (NV+2)*(NV+1)/2; % number of state

[W1,W2] = ndgrid(NV+1:-1:1,1:NV+1);

WV1 = tril(W1);
WV1 = WV1(WV1>0)-1;

WV2 = tril(W2);
WV2 = WV2(WV2>0)-1;

WV = [WV1;WV2];

%% Construct vibrational ladder operators 
% prep vib mode index in matrix form
[V1_Ind1,V1_Ind2] = ndgrid(WV1);
[V2_Ind1,V2_Ind2] = ndgrid(WV2);

% prep vib quantum # diff matrix
V1_qn_diff = V1_Ind1 - V1_Ind2;
V2_qn_diff = V2_Ind1 - V2_Ind2;

% prep nuclear wavefunction overlape matrix
V1_overlap = ~(V1_qn_diff);
V2_overlap = ~(V2_qn_diff);

% % FC using shifted parabolic approxi.
% V1_overlap = 1./sqrt(factorial(abs(V1_qn_diff))).*exp(-Lambda.^2./2).*Lambda.^(abs(V1_qn_diff));
% V2_overlap = 1./sqrt(factorial(abs(V2_qn_diff))).*exp(-Lambda.^2./2).*Lambda.^(abs(V2_qn_diff)); 

% Creation/Anihilation operator for |1> state 
B1_p_Ind = ~(V1_qn_diff-1);
B1_n_Ind = ~(V1_qn_diff+1);
B1_p = (      sqrt(WV1)*ones(size(WV1))') .* B1_p_Ind .* V2_overlap;
B1_n = (ones(size(WV1))*sqrt(WV1)'      ) .* B1_n_Ind .* V2_overlap;

% Creation/Anihilation operator for |2> state 
B2_p_Ind = ~(V2_qn_diff-1);
B2_n_Ind = ~(V2_qn_diff+1);
B2_p = (      sqrt(WV2)*ones(size(WV2))') .* B2_p_Ind .* V1_overlap;
B2_n = (ones(size(WV2))*sqrt(WV2)'      ) .* B2_n_Ind .* V1_overlap;

%% Construct Hmiltonian
% H = (H1) + (H2) + (H3) + (H4)
% H1: w0(Bn_p.*Bn_n)
% H2: lambda*w0*(Bn_p + Bn_n)|n><n|
% H3: J12*(|1><2| + |2><1|)
% H4: D + w0-0 + Lambda^2*w0

% H1
Bpn = B1_p*B1_n+B2_p*B2_n;
H1 = w0 * blkdiag(Bpn,Bpn);

% H2
H2 = Lambda * w0 * blkdiag(B1_p+B1_n,B2_p+B2_n);

% H3
H3 = reshape(eye(size(H1)),size(H1,1),NS,2); % need replace eye with FC factor matrix
H3 = reshape(H3(:,:,[2,1]),size(H1));
H3 = J12.*H3;

% H4
H4 = (D+w0_0+Lambda^2*w0).*eye(size(H1));

H = blkdiag(0,H1+H2+H3+H4);

%% Diagonalize H
% note: the eiganvector V_Full(:,i) has been already normalized.
[V_Full,D_Full] = eig(H);
Ex_Freq = diag(D_Full);

% sort eiganvalue form small to big and reorder the eiganvectors
[Sort_Ex_Freq,Indx] = sort(Ex_Freq);
 Sort_Ex_V          = V_Full(:,Indx);
 
%% Construct Transiton matrix

% considerign FC bt <g|v1,v2>
V1G_overlap = 1./sqrt(factorial(abs(WV1))).*(-Lambda).^(WV1).*exp(-Lambda.^2./2);
V2G_overlap = 1./sqrt(factorial(abs(WV2))).*(-Lambda).^(WV2).*exp(-Lambda.^2./2); 

Vib_overlap = V1G_overlap.*V2G_overlap;
% Vib_overlap = ones(NS,1);

Trans_Monent1 = bsxfun(@times,Vib_overlap,mu1);
Trans_Monent2 = bsxfun(@times,Vib_overlap,mu2);


Trans_Loc = zeros([size(H),3]);
Trans_Loc(1,   2:  NS+1,:) = Trans_Monent1;
Trans_Loc(1,NS+2:2*NS+1,:) = Trans_Monent2;
Trans_Loc(   2:  NS+1,1,:) = Trans_Monent1;
Trans_Loc(NS+2:2*NS+1,1,:) = Trans_Monent2;


Trans_Ex = zeros([size(H),3]);
% [Improve], maybe able to do more sophisticated
for C_mu=1:3
    Trans_Ex(:,:,C_mu) = Sort_Ex_V'*Trans_Loc(:,:,C_mu)*Sort_Ex_V;
end

%% Make figure
IntM = sum(Trans_Ex(2:end,1,:).^2,3);
IntMx = Trans_Ex(2:end,1,1).^2;
IntMy = Trans_Ex(2:end,1,2).^2;
IntMz = Trans_Ex(2:end,1,3).^2;

freq_OneD = Sort_Ex_Freq(2:end);

Num_Modes = length(freq_OneD);

% Get Frequency axis range
spec_range = F_Min:F_Max;

spec_array1 = bsxfun(@times,ones(Num_Modes,length(spec_range)),spec_range);
spec_array2 = bsxfun(@minus,spec_array1,freq_OneD);

switch LS
    case 'G' % Gaussian
        LineShape = exp(-(spec_array2.^2)./(LineWidth^2));
    case 'L' % Lorentzain 
        LineWidth = LineWidth/2;
        LineShape = LineWidth./((spec_array2.^2)+(LineWidth^2));
end


hF = figure; hold on

%% X 
if eq(PlotStick,1)
    line([freq_OneD';freq_OneD'],[zeros(1,Num_Modes);IntMx'],'Color','b')
end
CVLx = bsxfun(@times,LineShape,IntMx); 
CVLx_Total = sum(CVLx,1);
CVLx_Total = CVLx_Total.*max(abs(IntMx(:)))./max(abs(CVLx_Total));
plot(spec_range,CVLx_Total,'b-')

%% Y
if eq(PlotStick,1)
    line([freq_OneD';freq_OneD'],[zeros(1,Num_Modes);IntMy'],'Color','r')
end

CVLy = bsxfun(@times,LineShape,IntMy); 
CVLy_Total = sum(CVLy,1);
CVLy_Total = CVLy_Total.*max(abs(IntMy(:)))./max(abs(CVLy_Total)); 
plot(spec_range,CVLy_Total,'r-')

%% total
CVL = bsxfun(@times,LineShape,IntM); 
CVL_Total = sum(CVL,1);
CVL_Total = CVL_Total.*max(abs(IntM(:)))./max(abs(CVL_Total)); 
plot(spec_range,CVL_Total,'k:','linewidth',2)


hold off