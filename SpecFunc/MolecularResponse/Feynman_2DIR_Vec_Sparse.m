function [Grid,Freq,Int,Index,CutOff] = Feynman_2DIR_Vec_Sparse(PCutOff,FreqRange,EJR,F1,F2,M01,M12)
% 
% This function generate Feynman pathway of 2DIR and bin them into 1cm-1
% bucket
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
% Copyright Jia-Jung Ho, 2014-2016

%% debug
% GI = ParseGUI_Main(Data_COSMOSS.hGUIs);
% PCutOff = GI.PCutOff;
% FreqRange = GI.FreqRange;
% 
% R_Avg = LabFrameAvg('Isotropic',4); 
% J = JonesTrans4(pi/2,pi/2,pi/2,pi/2);
% E = EPolar4(0,0,0,0);
% EJR = E*J*R_Avg;
% 
% TwoDIR = Data_COSMOSS.TwoDIR;
% F1 = TwoDIR.H.Sort_Ex_F1;
% F2 = TwoDIR.H.Sort_Ex_F2;
% M01 = TwoDIR.Mu.M_Ex_01;
% M12 = TwoDIR.Mu.M_Ex_12;

%% Rund up frequency for sparse accumulation
% Round off so the sparse matrix accumulation will bin the response into
% 1cm-1 bucket. For finer bin modify this part
F1 = round(F1);
F2 = round(F2);

% Merge the first two index of Mu_Ex into one => 2D verion Mu_Ex
M12 = reshape(M12,[],3);

%% Prepare index
% Find number of modes in 1ex and 2ex
N1 = length(F1);
N2 = length(F2);

% R1 R2 NR1 NR2 index expansion,size N^2 
[Ib,Ia] = ndgrid(1:N1,1:N1);
Ia = Ia(:);
Ib = Ib(:);

% XYZ index expansion, size 3^4
[Jd,Jc,Jb,Ja] = ndgrid(1:3,1:3,1:3,1:3);
Ja = Ja(:);
Jb = Jb(:);
Jc = Jc(:);
Jd = Jd(:);

% R3 NR3 index expansion, size: N^3*(N+1)/2
[Kx,Kb,Ka] = ndgrid(1:N2,1:N1,1:N1);
Ka = Ka(:);
Kb = Kb(:);
Kx = Kx(:);

% Indexes of 2D version Mu_Ex, 
% Note [N,N*(N+1)/2] is the first two Dimension of Mu_Ex
Iax = sub2ind([N1,N2],Ka,Kx);
Ibx = sub2ind([N1,N2],Kb,Kx);

%% Reduce mode number base on the transition intensity
Norm_M_01 = sqrt(sum(M01.^2,2));
Norm_M_12 = sqrt(sum(M12.^2,2));

% Response intensity for (GB SE)
MI_0a = Norm_M_01(Ia);
MI_0b = Norm_M_01(Ib);
MI_a0 = MI_0a;
MI_b0 = MI_0b;

I_R1  = MI_b0.*MI_0b.*MI_a0.*MI_0a;
I_R2  = MI_b0.*MI_a0.*MI_0b.*MI_0a;
I_NR1 = MI_b0.*MI_0b.*MI_a0.*MI_0a;
I_NR2 = MI_a0.*MI_b0.*MI_0b.*MI_0a;

% Response intensity for EA
NI_0a = Norm_M_01(Ka);
NI_0b = Norm_M_01(Kb);

NI_ax = Norm_M_12(Iax);
NI_bx = Norm_M_12(Ibx);
NI_xa = NI_ax; % ind_xa = ind_ax
NI_xb = NI_bx; % ind_xb = ind_bx

I_R3  = NI_xa.*NI_bx.*NI_0b.*NI_0a;
I_NR3 = NI_xb.*NI_ax.*NI_0b.*NI_0a;

% Apply CutOff
R1_Max  = max(I_R1);
R2_Max  = max(I_R2);
R3_Max  = max(I_R3);
NR1_Max = max(I_NR1);
NR2_Max = max(I_NR2);
NR3_Max = max(I_NR3);

PI_R1  = I_R1  >= (PCutOff *  R1_Max);
PI_R2  = I_R2  >= (PCutOff *  R2_Max);
PI_R3  = I_R3  >= (PCutOff *  R3_Max);
PI_NR1 = I_NR1 >= (PCutOff * NR1_Max);
PI_NR2 = I_NR2 >= (PCutOff * NR2_Max);
PI_NR3 = I_NR3 >= (PCutOff * NR3_Max);

%% Get 2DIR response from index version kronecker product(GB SE)
% R1
Ib_R1   = Ib(PI_R1);
Ia_R1   = Ia(PI_R1);
S_R1    = ( M01(Ib_R1 ,Ja).*M01(Ib_R1 ,Jb).*M01(Ia_R1 ,Jc).*M01(Ia_R1 ,Jd) )*EJR';

Ea_1    = F1(Ia_R1);
Eb_1    = F1(Ib_R1);
Freq_R1 = [ Ea_1, zeros(size(Ea_1)), Eb_1]; % [pump,diff,probe]

% R2
Ib_R2   = Ib(PI_R2);
Ia_R2   = Ia(PI_R2);
S_R2    = ( M01(Ib_R2 ,Ja).*M01(Ia_R2 ,Jb).*M01(Ib_R2 ,Jc).*M01(Ia_R2 ,Jd) )*EJR';

Ea_2    = F1(Ia_R2);
Eb_2    = F1(Ib_R2);
Freq_R2 = [ Ea_2, Eb_2-Ea_2, Eb_2]; % [pump,diff,probe]

% NR1
Ib_NR1   = Ib(PI_NR1);
Ia_NR1   = Ia(PI_NR1);
S_NR1    = ( M01(Ib_NR1,Ja).*M01(Ib_NR1,Jb).*M01(Ia_NR1,Jc).*M01(Ia_NR1,Jd) )*EJR';

Ea_4     = F1(Ia_NR1);
Eb_4     = F1(Ib_NR1);
Freq_NR1 = [ Ea_4, zeros(size(Ea_4)), Eb_4]; % [pump,diff,probe]

% NR2
Ib_NR2   = Ib(PI_NR2);
Ia_NR2   = Ia(PI_NR2);
S_NR2    = ( M01(Ia_NR2,Ja).*M01(Ib_NR2,Jb).*M01(Ib_NR2,Jc).*M01(Ia_NR2,Jd) )*EJR';

Ea_5     = F1(Ia_NR2);
Eb_5     = F1(Ib_NR2);
Freq_NR2 = [ Ea_5, Ea_5-Eb_5, Ea_5]; % [pump,diff,probe]
% This is a interesting term! only contribute to diagonal!

%% Prep for R3/NR3
% estimate of largest array size and break it down to several for loop
MEM_CutOff = 3e-1; %[GB]
Ele_Max = round(MEM_CutOff/(81 * 8 / 1e9)) + 1; % number of elements to reach MEM_CufOff

%% R3
Iax_R3 = Iax(PI_R3);
Ibx_R3 = Ibx(PI_R3);
Ib_R3  =  Kb(PI_R3);
Ia_R3  =  Ka(PI_R3);

N_R3 = numel(Ia_R3);

if N_R3 > Ele_Max
    % Add NaNs so I can reshape indexes
    Padding_L = Ele_Max - mod(N_R3,Ele_Max);
    Padding_NaN = nan(Padding_L,1);
    
    Ia_R3_L  = reshape([ Ia_R3; Padding_NaN],Ele_Max,[]);
    Ib_R3_L  = reshape([ Ib_R3; Padding_NaN],Ele_Max,[]);
    Iax_R3_L = reshape([Iax_R3; Padding_NaN],Ele_Max,[]);
    Ibx_R3_L = reshape([Ibx_R3; Padding_NaN],Ele_Max,[]);

    Loop_N = size(Ia_R3_L,2);
    % display how many loop is going to be done
    MEM = numel(Ia_R3_L) * 81 * 8 / 1e9; % roughly unit in GB
    disp('--------------------------------------')
    disp(['Memory cut-off set at ', num2str(MEM_CutOff), 'GB...'])
    disp(['Memory size of R3 is about ', sprintf('%4.1f',MEM), 'GB...'])
    disp(['Will run R3 calculation into ', num2str(Loop_N), ' Loops to reduce memory load...'])
    disp('--------------------------------------')

    % pre-allocate matrix and run the loop
    S_R3_L = zeros(Ele_Max,Loop_N-1);
    for LN = 1:Loop_N -1 
        % note pre-allocate matrxi expansion to save calculation time
        P4 = M12(Iax_R3_L(:,LN),:);
        P3 = M12(Ibx_R3_L(:,LN),:);
        P2 = M01( Ib_R3_L(:,LN),:);
        P1 = M01( Ia_R3_L(:,LN),:);
        
        S_R3_L(:,LN) = ( P4(:,Ja).*P3(:,Jb).*P2(:,Jc).*P1(:,Jd) )*EJR';
    end    
    S_R3_L = S_R3_L(:);

    % deall with the last loop
    Last_Ind   = (1:mod(N_R3,Ele_Max))';
    Ia_R3_End  =  Ia_R3_L(Last_Ind,Loop_N);
    Ib_R3_End  =  Ib_R3_L(Last_Ind,Loop_N);
    Iax_R3_End = Iax_R3_L(Last_Ind,Loop_N);
    Ibx_R3_End = Ibx_R3_L(Last_Ind,Loop_N);

    P4 = M12(Iax_R3_End,:);
    P3 = M12(Ibx_R3_End,:);
    P2 = M01( Ib_R3_End,:);
    P1 = M01( Ia_R3_End,:);
    
    S_R3_End = ( P4(:,Ja).*P3(:,Jb).*P2(:,Jc).*P1(:,Jd) )*EJR';
        
    S_R3 = [ S_R3_L; S_R3_End];    
    
else
    P4 = M12(Iax_R3,:);
    P3 = M12(Ibx_R3,:);
    P2 = M01( Ib_R3,:);
    P1 = M01( Ia_R3,:);
    S_R3 = ( P4(:,Ja).*P3(:,Jb).*P2(:,Jc).*P1(:,Jd) )*EJR';
end

Ea_3 = F1(Ia_R3);
Eb_3 = F1(Ib_R3);
Ex_3 = F2(Kx(PI_R3));
Freq_R3 = [ Ea_3 , Eb_3 - Ea_3 , Ex_3 - Ea_3 ];

%% NR3
Iax_NR3 = Iax(PI_NR3);
Ibx_NR3 = Ibx(PI_NR3);
Ib_NR3  =  Kb(PI_NR3);
Ia_NR3  =  Ka(PI_NR3);

NR3 = numel(Ia_NR3);

if NR3 > Ele_Max
    % Add NaNs so I can reshape indexes
    Padding_L = Ele_Max - mod(NR3,Ele_Max);
    Padding_NaN = nan(Padding_L,1);
    
    Ia_NR3_L  = reshape([ Ia_NR3; Padding_NaN],Ele_Max,[]);
    Ib_NR3_L  = reshape([ Ib_NR3; Padding_NaN],Ele_Max,[]);
    Iax_NR3_L = reshape([Iax_NR3; Padding_NaN],Ele_Max,[]);
    Ibx_NR3_L = reshape([Ibx_NR3; Padding_NaN],Ele_Max,[]);

    Loop_N = size(Ia_NR3_L,2);
    % display how many loop is going to be done
    MEM = numel(Ia_NR3_L) * 81 * 8 / 1e9; % roughly unit in GB
    disp('--------------------------------------')
    disp(['Memory cut-off set at ', num2str(MEM_CutOff), 'GB...'])
    disp(['Memory size of NR3 is about ', sprintf('%4.1f',MEM), 'GB...'])
    disp(['Will run NR3 calculation into ', num2str(Loop_N), ' Loops to reduce memory load...'])
    disp('--------------------------------------')

    % pre-allocate matrix and run the loop
    S_NR3_L = zeros(Ele_Max,Loop_N-1);
    for LN = 1:Loop_N -1 
        % note pre-allocate matrxi expansion to save calculation time
        P4 = M12(Ibx_NR3_L(:,LN),:);
        P3 = M12(Iax_NR3_L(:,LN),:);
        P2 = M01( Ib_NR3_L(:,LN),:);
        P1 = M01( Ia_NR3_L(:,LN),:);
        
        S_NR3_L(:,LN) = ( P4(:,Ja).*P3(:,Jb).*P2(:,Jc).*P1(:,Jd) )*EJR';
    end    
    S_NR3_L = S_NR3_L(:);

    % deall with the last loop
    Last_Ind    = (1:mod(NR3,Ele_Max))';
    Ia_NR3_End  =  Ia_NR3_L(Last_Ind,Loop_N);
    Ib_NR3_End  =  Ib_NR3_L(Last_Ind,Loop_N);
    Iax_NR3_End = Iax_NR3_L(Last_Ind,Loop_N);
    Ibx_NR3_End = Ibx_NR3_L(Last_Ind,Loop_N);

    P4 = M12(Ibx_NR3_End,:);
    P3 = M12(Iax_NR3_End,:);
    P2 = M01( Ib_NR3_End,:);
    P1 = M01( Ia_NR3_End,:);
    
    S_NR3_End = ( P4(:,Ja).*P3(:,Jb).*P2(:,Jc).*P1(:,Jd) )*EJR';
        
    S_NR3 = [ S_NR3_L; S_NR3_End];    
    
else
    P4 = M12(Ibx_NR3,:);
    P3 = M12(Iax_NR3,:);
    P2 = M01( Ib_NR3,:);
    P1 = M01( Ia_NR3,:);
    
    S_NR3 = ( P4(:,Ja).*P3(:,Jb).*P2(:,Jc).*P1(:,Jd) )*EJR';
end

Ea_6 = F1(Ia_NR3);
Eb_6 = F1(Ib_NR3);
Ex_6 = F2(Kx(PI_NR3));
Freq_NR3 = [ Ea_6 , Ea_6 - Eb_6 , Ex_6 - Eb_6 ];

%% Prep for export mode index
Ind_R1  = [Ib_R1,Ib_R1,Ia_R1,Ia_R1]; % GB
Ind_R2  = [Ib_R2,Ia_R2,Ib_R2,Ia_R2]; % SE
Ind_R3  = [Iax_R3,Ibx_R3,Ib_R3,Ia_R3]; % EA

Ind_NR1 = [Ib_NR1,Ib_NR1,Ia_NR1,Ia_NR1]; % GB
Ind_NR2 = [Ia_NR2,Ib_NR2,Ib_NR2,Ia_NR2]; % SE
Ind_NR3 = [Iax_NR3,Ibx_NR3,Ib_NR3,Ia_NR3]; % EA

%% Construct sparse matrix version of Responses
% This part also deal with binning
SparseMax = max([max(F1), max(F2(Kx(:)) - F1(Ka(:))), max(FreqRange)]); 

R1  = sparse( Freq_R1(:,1), Freq_R1(:,3),S_R1 ,SparseMax,SparseMax);
R2  = sparse( Freq_R2(:,1), Freq_R2(:,3),S_R2 ,SparseMax,SparseMax);
R3  = sparse( Freq_R3(:,1), Freq_R3(:,3),S_R3 ,SparseMax,SparseMax);
NR1 = sparse(Freq_NR1(:,1),Freq_NR1(:,3),S_NR1,SparseMax,SparseMax);
NR2 = sparse(Freq_NR2(:,1),Freq_NR2(:,3),S_NR2,SparseMax,SparseMax);
NR3 = sparse(Freq_NR3(:,1),Freq_NR3(:,3),S_NR3,SparseMax,SparseMax);

%% Output
% MinF = min(FreqRange);
% MaxF = max(FreqRange);
% 
% Grid.R1  = full( R1(MinF:MaxF,MinF:MaxF)); % [Pumpx Probe]
% Grid.R2  = full( R2(MinF:MaxF,MinF:MaxF)); % [Pumpx Probe]
% Grid.R3  = full( R3(MinF:MaxF,MinF:MaxF)); % [Pumpx Probe]
% Grid.NR1 = full(NR1(MinF:MaxF,MinF:MaxF)); % [Pumpx Probe]
% Grid.NR2 = full(NR2(MinF:MaxF,MinF:MaxF)); % [Pumpx Probe]
% Grid.NR3 = full(NR3(MinF:MaxF,MinF:MaxF)); % [Pumpx Probe]

Grid.R1  = R1;  % [Pumpx Probe]
Grid.R2  = R2;  % [Pumpx Probe]
Grid.R3  = R3;  % [Pumpx Probe]
Grid.NR1 = NR1; % [Pumpx Probe]
Grid.NR2 = NR2; % [Pumpx Probe]
Grid.NR3 = NR3; % [Pumpx Probe]

Int.R1  = S_R1;
Int.R2  = S_R2;
Int.R3  = S_R3;
Int.NR1 = S_NR1;
Int.NR2 = S_NR2;
Int.NR3 = S_NR3;

Freq.R1  = Freq_R1;
Freq.R2  = Freq_R2;
Freq.R3  = Freq_R3;
Freq.NR1 = Freq_NR1;
Freq.NR2 = Freq_NR2;
Freq.NR3 = Freq_NR3;

Index.R1  = Ind_R1;
Index.R2  = Ind_R2;
Index.R3  = Ind_R3;
Index.NR1 = Ind_NR1;
Index.NR2 = Ind_NR2;
Index.NR3 = Ind_NR3;

CutOff.PI_R1   = PI_R1;
CutOff.PI_R2   = PI_R2;
CutOff.PI_R3   = PI_R3;
CutOff.PI_NR1  = PI_NR1;
CutOff.PI_NR2  = PI_NR2;
CutOff.PI_NR3  = PI_NR3;
CutOff.PCutOff = PCutOff;
