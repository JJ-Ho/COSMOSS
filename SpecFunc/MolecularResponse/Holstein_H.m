% Holstein-like Hamiltonian 
% This hamiltonian consist of a dimer system compose of identical monomer. 
% This system only invole 1 electronic state 1 vibrational state on the 1st
% electronic exicted state for each monomer. 
% Ref: Kistler, K. A.; Pochas, C. M.; Yamagata, H. JPCB, 2011, 116, 77?86.

%% Construct local mode basis
NV = 1; % NUmber of maximun vibrational quata
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

w0 =1600;
Lambda = 0.57;
J12 = 300;
D = 0;
w0_0 = 10000;

% H1
H1 = w0 * blkdiag(B1_p*B1_n,B2_p*B2_n);

% H2
H2 = Lambda * w0 * blkdiag(B1_p+B1_n,B2_p+B2_n);

% H3
H3 = reshape(eye(size(H1)),size(H1,1),NS,2); % need replace eye with FC factor matrix
H3 = reshape(H3(:,:,[2,1]),size(H1));
H3 = J12.*H3;

% H4
H4 = (D+w0_0+Lambda^2*w0).*eye(size(H1));

H = blkdiag(0,H1+H2+H3+H4);




