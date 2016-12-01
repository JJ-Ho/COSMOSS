function Output=Feynman_2DSFG_kron(N,Sort_Ex_Freq,alpha_Ex,mu_Ex)
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
% Ver. 1.0  131201  Use build-in function 'kron' to do tensor product
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2013

% %% debug
% N = LocModeNum;

%%
NumP = 3^5;
N12Path = N^2;
N3Path = N^3*(N+1)/2;

R1  = zeros(N12Path,NumP);
R2  = zeros(N12Path,NumP);
NR1 = zeros(N12Path,NumP);
NR2 = zeros(N12Path,NumP);

R3  = zeros(N3Path,NumP);
NR3 = zeros(N3Path,NumP);

% R1_freq = zeros(N12Path,3);
% R2_freq = zeros(N12Path,3);
% NR1_freq = zeros(N12Path,3);
% NR2_freq = zeros(N12Path,3);
% 
% R3_freq = zeros(N3Path,3);
% NR3_freq = zeros(N3Path,3);
% 

for Na=1:N
    for Nb=1:N
        
            Ia = Na+1;
            Ib = Nb+1;
            
            %Ea = Sort_Ex_Freq(Ia);
            %Eb = Sort_Ex_Freq(Ib);
            
            ISave = Nb+(Na-1)*N;
            
            %R1(ISave,1:3) =[-Ea,0    ,Eb];
            %R2(ISave,1:3) =[-Ea,Eb-Ea,Eb];
            %NR1(ISave,1:3)=[ Ea,0    ,Eb];
            %NR2(ISave,1:3)=[ Ea,Ea-Eb,Ea];
            

            A_a0 = alpha_Ex(Ia,1,:);
            A_b0 = alpha_Ex(Ib,1,:);
            
            M_a0 = mu_Ex(Ia,1,:);
            M_b0 = mu_Ex(Ib,1,:);
            
            M_0a = mu_Ex(1,Ia,:);
            M_0b = mu_Ex(1,Ib,:);

            
            R1(ISave,:)  = kron(kron(kron(A_b0(:),M_0b(:)),M_a0(:)),M_0a(:));
            R2(ISave,:)  = kron(kron(kron(A_b0(:),M_0a(:)),M_b0(:)),M_0a(:));
            
            NR1(ISave,:) = kron(kron(kron(A_b0(:),M_0b(:)),M_a0(:)),M_0a(:));
            NR2(ISave,:) = kron(kron(kron(A_a0(:),M_0b(:)),M_b0(:)),M_0a(:));
            
        
        for Nx=1:N*(N+1)/2
            
            Ix = Nx+1+N;
            %Ex = Sort_Ex_Freq(Ix);
            ISave = Nx+(Nb-1)*N*(N+1)/2+(Na-1)*N^2*(N+1)/2;
            
            %R3(ISave,1:3) =[-Ea,Eb-Ea,Ex-Ea];
            %NR3(ISave,1:3)=[ Ea,Ea-Eb,Ex-Eb];
            
            
            A_ax = alpha_Ex(Ia,Ix,:);
            A_bx = alpha_Ex(Ib,Ix,:);
            
            M_xa = mu_Ex(Ix,Ia,:);
            M_xb = mu_Ex(Ix,Ib,:);

            R3(ISave,:)  = kron(kron(kron(A_ax(:),M_xb(:)),M_b0(:)),M_0a(:));
            NR3(ISave,:) = kron(kron(kron(A_bx(:),M_xa(:)),M_b0(:)),M_0a(:));
            
        end
    end
end

%% Generate List of interaction Frequencies

% R1 R2 NR1 NR2 index expansion,size N^2 
[Ib,Ia] = ndgrid(2:N+1,2:N+1);
% R3 NR3 index expansion, size: N^3*(N+1)/2
[Kx,Kb,Ka] = ndgrid(2+N:(N+1)*(N+2)/2,2:N+1,2:N+1);

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
Output.R1=R1';
Output.R2=R2';
Output.R3=R3';
Output.NR1=NR1';
Output.NR2=NR2';
Output.NR3=NR3';

Output.Freq_R1  = Freq_R1;
Output.Freq_R2  = Freq_R2;
Output.Freq_R3  = Freq_R3;
Output.Freq_NR1 = Freq_NR1;
Output.Freq_NR2 = Freq_NR2;
Output.Freq_NR3 = Freq_NR3;



% Define optimized kron product function
function Product = kron(T1,T2)
    Product_matrix = T1(:)*T2(:)';
    Product = Product_matrix(:)';

