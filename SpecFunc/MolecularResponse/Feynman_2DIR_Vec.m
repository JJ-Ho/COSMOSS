function Response = Feynman_2DIR_Vec(N,Sort_Ex_Freq,Mu_Ex)
% 
% This function generate Feynman pathway of 2DSFG with given polarization.
% 
% Since "kron" function use reshape alot, the overall speed of kron product
% version is slower than for loop version
% 
% Todo: code acceleration, delete small signals.
% 
% ------- Version log -----------------------------------------------------
% 
% Ver. 1.0  140126  modified from Feynman_2DIR_kron; vecterized version
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2014

%% debug
% N = NumLocMode;

%% Prepare index

% R1 R2 NR1 NR2 index expansion,size N^2 
[Ib,Ia] = ndgrid(2:N+1,2:N+1);
% XYZ index expansion, size 3^4
[Jd,Jc,Jb,Ja] = ndgrid(1:3,1:3,1:3,1:3);
% R3 NR3 index expansion, size: N^3*(N+1)/2
[Kx,Kb,Ka] = ndgrid(2+N:(N+1)*(N+2)/2,2:N+1,2:N+1);

%% Get 2DIR response from index version kronec product

% R1 R2 NR1 NR2
M_0a = squeeze(Mu_Ex(1,Ia(:),:));
M_0b = squeeze(Mu_Ex(1,Ib(:),:));
M_a0 = squeeze(Mu_Ex(Ia(:),1,:));
M_b0 = squeeze(Mu_Ex(Ib(:),1,:));

% Response (GB SE)
R1  = M_b0(:,Ja(:)).*M_0b(:,Jb(:)).*M_a0(:,Jc(:)).*M_0a(:,Jd(:));
R2  = M_b0(:,Ja(:)).*M_a0(:,Jb(:)).*M_0b(:,Jc(:)).*M_0a(:,Jd(:));
NR1 = M_b0(:,Ja(:)).*M_0b(:,Jb(:)).*M_a0(:,Jc(:)).*M_0a(:,Jd(:));
NR2 = M_a0(:,Ja(:)).*M_b0(:,Jb(:)).*M_0b(:,Jc(:)).*M_0a(:,Jd(:));


% R3 NR3
N_0a = squeeze(Mu_Ex(1,Ka,:));
N_0b = squeeze(Mu_Ex(1,Kb,:));

% Merge the first two index of Mu_Ex into one => 2D verion Mu_Ex
TD_Mu_Ex = reshape(Mu_Ex,[],3);

% Indice of 2D version Mu_Ex, 
% Note [(N+1)*(N+2)/2,(N+1)*(N+2)/2] is the first two Dimension of Mu_Ex

ind_ax = sub2ind([(N+1)*(N+2)/2,(N+1)*(N+2)/2],Ka(:),Kx(:));
ind_bx = sub2ind([(N+1)*(N+2)/2,(N+1)*(N+2)/2],Kb(:),Kx(:));
ind_xa = sub2ind([(N+1)*(N+2)/2,(N+1)*(N+2)/2],Kx(:),Ka(:));
ind_xb = sub2ind([(N+1)*(N+2)/2,(N+1)*(N+2)/2],Kx(:),Kb(:));

% get TDV using linear index of the first two indice of Mu_Ex 
% and ":" to extract the whole vector components 
N_ax = TD_Mu_Ex(ind_ax,:);
N_bx = TD_Mu_Ex(ind_bx,:);
N_xa = TD_Mu_Ex(ind_xa,:);
N_xb = TD_Mu_Ex(ind_xb,:);

% Response (EA)
R3  = N_xa(:,Ja(:)).*N_bx(:,Jb(:)).*N_0b(:,Jc(:)).*N_0a(:,Jd(:));
NR3 = N_xb(:,Ja(:)).*N_ax(:,Jb(:)).*N_0b(:,Jc(:)).*N_0a(:,Jd(:));

%% Generate List of interaction Frequencies

Ea_12 = Sort_Ex_Freq(Ia(:))';
Eb_12 = Sort_Ex_Freq(Ib(:))';

Freq_R1  = [-Ea_12 ; zeros(1,N^2)  ; Eb_12]';
Freq_R2  = [-Ea_12 ; Eb_12 - Ea_12 ; Eb_12]';
Freq_NR1 = [ Ea_12 ; zeros(1,N^2)  ; Eb_12]';
Freq_NR2 = [ Ea_12 ; Ea_12 - Eb_12 ; Ea_12]'; % This is a interesting term! only contribute to diagonal!

Ea_3 = Sort_Ex_Freq(Ka(:))';
Eb_3 = Sort_Ex_Freq(Kb(:))';
Ex_3 = Sort_Ex_Freq(Kx(:))';

Freq_R3  = [-Ea_3 ; Eb_3 - Ea_3 ; Ex_3 - Ea_3]';
Freq_NR3 = [ Ea_3 ; Ea_3 - Eb_3 ; Ex_3 - Eb_3]';

%% Output
% Response.R1  = R1;
% Response.R2  = R2;
% Response.R3  = R3;
% Response.NR1 = NR1;
% Response.NR2 = NR2;
% Response.NR3 = NR3;

% Frequency.Freq_R1  = Freq_R1;
% Frequency.Freq_R2  = Freq_R2;
% Frequency.Freq_R3  = Freq_R3;
% Frequency.Freq_NR1 = Freq_NR1;
% Frequency.Freq_NR2 = Freq_NR2;
% Frequency.Freq_NR3 = Freq_NR3;

Response.R1  = [Freq_R1,R1];
Response.R2  = [Freq_R2,R2];
Response.R3  = [Freq_R3,R3];
Response.NR1 = [Freq_NR1,NR1];
Response.NR2 = [Freq_NR2,NR2];
Response.NR3 = [Freq_NR3,NR3];

%% Old for loop version
% 
% NumP = 3^4;
% N12Path = N^2;
% N3Path = N^3*(N+1)/2;
% 
% 
% RT1  = zeros(N12Path,NumP);
% RT2  = zeros(N12Path,NumP);
% NRT1 = zeros(N12Path,NumP);
% NRT2 = zeros(N12Path,NumP);
% 
% RT3  = zeros(N3Path,NumP);
% NRT3 = zeros(N3Path,NumP);
% 
% tic
% for Na=1:N
%     for Nb=1:N
%         
%             Ia = Na+1;
%             Ib = Nb+1;
%             
% %             Ea = Sort_Ex_Freq(Ia);
% %             Eb = Sort_Ex_Freq(Ib);
%             
%             ISave = Nb+(Na-1)*N;
%             
% %             R1(ISave,1:3) =[-Ea,0    ,Eb];
% %             R2(ISave,1:3) =[-Ea,Eb-Ea,Eb];
% %             NR1(ISave,1:3)=[ Ea,0    ,Eb];
% %             NR2(ISave,1:3)=[ Ea,Ea-Eb,Ea];
% %             
% 
% %             A_a0 = AAAA(Ia,1,:);
% %             A_b0 = AAAA(Ib,1,:);
%             
%             M_a0 = Mu_Ex(Ia,1,:);
%             M_b0 = Mu_Ex(Ib,1,:);
%             
%             M_0a = Mu_Ex(1,Ia,:);
%             M_0b = Mu_Ex(1,Ib,:);
% 
%             
%             RT1(ISave,:)  = kron(kron(kron(M_b0(:),M_0b(:)),M_a0(:)),M_0a(:));
%             RT2(ISave,:)  = kron(kron(kron(M_b0(:),M_0a(:)),M_b0(:)),M_0a(:));
%             
%             NRT1(ISave,:) = kron(kron(kron(M_b0(:),M_0b(:)),M_a0(:)),M_0a(:));
%             NRT2(ISave,:) = kron(kron(kron(M_a0(:),M_0b(:)),M_b0(:)),M_0a(:));
%             
%         
%         for Nx=1:N*(N+1)/2
%             
%             Ix = Nx+1+N;
%             Ex = Sort_Ex_Freq(Ix);
%             ISave = Nx+(Nb-1)*N*(N+1)/2+(Na-1)*N^2*(N+1)/2;
%             
% %             R3(ISave,1:3) =[-Ea,Eb-Ea,Ex-Ea];
% %             NR3(ISave,1:3)=[ Ea,Ea-Eb,Ex-Eb];
%             
%             
% %             A_ax = AAAA(Ia,Ix,:);
% %             A_bx = AAAA(Ib,Ix,:);
%             
%             M_ax = Mu_Ex(Ia,Ix,:);
%             M_bx = Mu_Ex(Ib,Ix,:);
% 
%             M_xa = Mu_Ex(Ix,Ia,:);
%             M_xb = Mu_Ex(Ix,Ib,:);
% 
%             RT3(ISave,:)  = kron(kron(kron(M_xa(:),M_bx(:)),M_0b(:)),M_0a(:));
%             NRT3(ISave,:) = kron(kron(kron(M_xb(:),M_ax(:)),M_0b(:)),M_0a(:));
%             
%         end
%     end
% end
% toc


% % Define optimized kron product function for old for-loop version
% function Product = kron(T1,T2)
%     Product_matrix = T1(:)*T2(:)';
%     Product = Product_matrix(:)';
% 
