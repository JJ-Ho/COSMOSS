function Output=Feynman_2DIR_kron(N,Sort_Ex_Freq,mu_Ex)
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
% Ver. 1.0  131209  modified from Feynman_2DSFG_kron
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2013

% %% debug
% N = LocModeNum;

%%
NumP = 3^4;
N12Path = N^2;
N3Path = N^3*(N+1)/2;

R1  = zeros(N12Path,NumP+3);
R2  = zeros(N12Path,NumP+3);
NR1 = zeros(N12Path,NumP+3);
NR2 = zeros(N12Path,NumP+3);

R3  = zeros(N3Path,NumP+3);
NR3 = zeros(N3Path,NumP+3);


for Na=1:N
    for Nb=1:N
        
            Ia = Na+1;
            Ib = Nb+1;
            
            Ea = Sort_Ex_Freq(Ia);
            Eb = Sort_Ex_Freq(Ib);
            
            ISave = Nb+(Na-1)*N;
            
            R1(ISave,1:3) =[-Ea,0    ,Eb];
            R2(ISave,1:3) =[-Ea,Eb-Ea,Eb];
            NR1(ISave,1:3)=[ Ea,0    ,Eb];
            NR2(ISave,1:3)=[ Ea,Ea-Eb,Ea];
            
            M_a0 = mu_Ex(Ia,1,:);
            M_b0 = mu_Ex(Ib,1,:);
            
            M_0a = mu_Ex(1,Ia,:);
            M_0b = mu_Ex(1,Ib,:);

            
            R1(ISave,4:end)  = kron(kron(kron(M_b0(:),M_0b(:)),M_a0(:)),M_0a(:));
            R2(ISave,4:end)  = kron(kron(kron(M_b0(:),M_0a(:)),M_b0(:)),M_0a(:));
            
            NR1(ISave,4:end) = kron(kron(kron(M_b0(:),M_0b(:)),M_a0(:)),M_0a(:));
            NR2(ISave,4:end) = kron(kron(kron(M_a0(:),M_0b(:)),M_b0(:)),M_0a(:));
            
        
        for Nx=1:N*(N+1)/2
            
            Ix = Nx+1+N;
            Ex = Sort_Ex_Freq(Ix);
            ISave = Nx+(Nb-1)*N*(N+1)/2+(Na-1)*N^2*(N+1)/2;
            
            R3(ISave,1:3) =[-Ea,Eb-Ea,Ex-Ea];
            NR3(ISave,1:3)=[ Ea,Ea-Eb,Ex-Eb];
            
            
            M_ax = mu_Ex(Ia,Ix,:);
            M_bx = mu_Ex(Ib,Ix,:);

            M_xa = mu_Ex(Ix,Ia,:);
            M_xb = mu_Ex(Ix,Ib,:);

            R3(ISave,4:end)  = kron(kron(kron(M_ax(:),M_xb(:)),M_b0(:)),M_0a(:));
            NR3(ISave,4:end) = kron(kron(kron(M_bx(:),M_xa(:)),M_b0(:)),M_0a(:));
            
        end
    end
end


%% Output
Output.R1=R1;
Output.R2=R2;
Output.R3=R3;
Output.NR1=NR1;
Output.NR2=NR2;
Output.NR3=NR3;


% Define optimized kron product function
function Product = kron(T1,T2)
    Product_matrix = T1(:)*T2(:)';
    Product = Product_matrix(:)';

