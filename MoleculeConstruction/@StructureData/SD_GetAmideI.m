function Obj_SD_AmideI = SD_GetAmideI(obj_SD)
% This function recongnize the amide I group by searching C-O-N-CA atom
% name sequence. 
% Copyright Jia-Jung Ho, 2013-2020
% obj_SD = Data_PDB_AmideI.app.PDB;
%% Redefine the variable names
XYZ       = obj_SD.XYZ;
AtomName  = obj_SD.AtomName;
ChainID   = obj_SD.Extra.ChainID;

% [C,O,N,CA]
Atom_C  = strcmp(AtomName,'C'); %size(Atom_C)
Atom_O  = strcmp(AtomName,'O'); %size(Atom_O)
Atom_N  = strcmp(AtomName,'N'); %size(Atom_N)
Atom_CA = strcmp(AtomName,'CA'); %size(Atom_CA)
Atom = [Atom_C,Atom_O,Atom_N,Atom_CA];
AtomInd = (1:obj_SD.NAtoms)';

%% Check if there is ChainIDs
if ~isspace([ChainID{:}])
    
    Chain_Ind_TF = ~strcmp(ChainID(2:end),ChainID(1:end-1));
    ChainEnd     = AtomInd(Chain_Ind_TF);
    ChainStart   = ChainEnd+1;
    ChainEnd     = [ChainEnd;obj_SD.NAtoms];
    ChainStart   = [1;ChainStart];
    N_Chain      = length(ChainStart);
else
    N_Chain = 1;
    ChainStart = 1;
    ChainEnd = obj_SD.NAtoms;
end

%% Find AmideI mode in each Chain and determine their atom index and "AmideIAtomSerNo"
AmideIAtomSerNo = [];

for NC = 1:N_Chain
    AtomInd = (ChainStart(NC):ChainEnd(NC))';
    ChainAtom = Atom(ChainStart(NC):ChainEnd(NC),:);
    
    % delete  irrelevent lines with no above atoms 
    AtomInd (~any(ChainAtom,2))   = [];
    ChainAtom(~any(ChainAtom,2),:) = [];
    % Remove first few rows not start with 'C'
    Head   = 'y';
    while strcmp(Head,'y')
        if eq(ChainAtom(1,:), [1,0,0,0])
            Head = 'n';
        else
            ChainAtom(1,:) = [];
             AtomInd(1,:) = [];
        end
    end

    % Remove last few rows not end with 'CA'
    Tail   = 'y';
    while strcmp(Tail,'y')
        if eq(ChainAtom(end,:), [0,0,0,1])
            Tail = 'n';
        else
            ChainAtom(end,:) = [];
             AtomInd(end,:) = [];
        end
    end
    % Remove any lines that does not have all four atoms
    % ex:
    %   i [1,0,0,0]     [1,0,0,0]
    %   j [0,1,0,0]  vs [0,1,0,0]  by clapse the ijk dimemsion
    %   k [0,0,1,0]     [0,1,0,0]
    %   l [0,0,0,1]     [0,0,0,1]
    %        ||            ||
    %        \/            \/
    %     [1,1,1,1]     [1,2,0,1]  

    AtomX = sum(permute(reshape(ChainAtom,4,[],4),[3,2,1]),3);
    Mask = ones(size(AtomX));
    Array_Mask = logical(bsxfun(@times,all(~(AtomX - Mask)),[1,1,1,1]'));
    IndX = AtomInd(Array_Mask(:));  
    Tmp_AmideIAtomSerNo = reshape(IndX,4,[])';
    AmideIAtomSerNo = [AmideIAtomSerNo;Tmp_AmideIAtomSerNo];
end

%% Define Aminde I modes coordinate system
% Find the XYZ coordinate of Amide I related Atoms
Nmodes = size(AmideIAtomSerNo,1);
XYZ_AmideI_Atoms = XYZ(AmideIAtomSerNo(:),:);
XYZ_AmideI_Atoms = reshape(XYZ_AmideI_Atoms,Nmodes,[],3); % [C_xyz O_xyz N_xyz CA_xyz; C_xyz, O_xyz, N_xyz, CA_xyz ; ...]

% Define local coordinate system
Vec_CO = XYZ_AmideI_Atoms(:,2,:)-XYZ_AmideI_Atoms(:,1,:);
Vec_CO = reshape(Vec_CO,Nmodes,3);
Vec_CO = bsxfun(@rdivide,Vec_CO,sqrt(sum(abs(Vec_CO).^2,2))); % normaliz CO vectors

Vec_CN = XYZ_AmideI_Atoms(:,3,:)-XYZ_AmideI_Atoms(:,1,:);
Vec_CN=  reshape(Vec_CN,Nmodes,3); 
Vec_CN = bsxfun(@rdivide,Vec_CN,sqrt(sum(abs(Vec_CN).^2,2))); % normaliz CN vectors

AmideICenter = reshape(XYZ_AmideI_Atoms(:,1,:),Nmodes,3) + Vec_CO.*0.665 + Vec_CN.*0.256; % center of amide I mode, ref from Jenny's mathematica code

% Define Lab frame coordinate of each mode
Z_Sim = Vec_CO;
X_Sim = cross(Vec_CN,Z_Sim,2);
X_Sim = bsxfun(@rdivide,X_Sim,sqrt(sum(abs(X_Sim).^2,2))); % normalize
Y_Sim = cross(Z_Sim,X_Sim,2);
Y_Sim = bsxfun(@rdivide,Y_Sim,sqrt(sum(abs(Y_Sim).^2,2))); % normalize

XYZ_Sim = [X_Sim(:); Y_Sim(:); Z_Sim(:)]';
XYZ_Sim = reshape(XYZ_Sim,Nmodes,3,[]); 

%% Calculate the transition dipoles (mu) and Raman tensors (alpha) in the Lab frame
% Transition dipole Angle
mu_angle = -27.5*pi/180; % from J. Laaser's paper (10.1021/jp408064b)
% mu_angle = -20*pi/180; % from D. Strasfeld's paper (10.1021/jp9072203)
% mu_angle = -25*pi/180; % from M. H. Cho's paper
% mu_angle = -10*pi/180; % follow J. Skinner's paper (doi:10.1063/1.4882059)

mu_vec0 = 16.1016*[0,0,-1]';       % this scale factor was calculated from DFT simulation, check the 'NMA' note 
mu_Mol  = (Rx(mu_angle)*mu_vec0)'; % 1 by 3

% Raman tensor
alpha_angle = -34*pi/180;    % counter clockwise roation about x axis from Jenny's paper (10.1021/jp408064b)
alpha0 = 5*[ 0.05, 0  , 0;   
             0   , 0.2, 0; 
             0   , 0  , 1 ];

alpha_Mol = Rx(alpha_angle) * alpha0 * Rx(alpha_angle)';

mu_Sim    = zeros(Nmodes,3);
alpha_Sim = zeros(Nmodes,3,3);

for ii=1:Nmodes
    mu_Sim(ii,:)      = (squeeze(XYZ_Sim(ii,:,:))*mu_Mol')';
    alpha_Sim(ii,:,:) = squeeze(XYZ_Sim(ii,:,:))*alpha_Mol'*squeeze(XYZ_Sim(ii,:,:))';
end    

% take vectorize Alpha
% alpha_Sim = [N x 3 x3 ]
% for a signle mode, the alpha_Sim: 
% [ XX, XY, XZ ]
% [ YX, YY, YZ ]
% [ ZX, ZY, ZZ ]
% Following the kron convention that the following E J L R beta used, the
% alpha vector need to be in this index order
% [XX,XY, XZ, YX, YY, YZ, ZX, ZY, ZZ]'
% eventhough the Raman tensor we encounter in IR resonance SFG is always 
% symmetric, I am being exta caucious here to make the indexing right.
alpha = reshape(permute(alpha_Sim,[1,3,2]),[Nmodes,9]);

%% Define Mode frequency and anharmonicity
% Define mode frequeny with default values
NLFreq = 1645;
Anharm = 20;

AmideIFreq   = ones(Nmodes,1).* NLFreq;
AmideIAnharm = ones(Nmodes,1).* Anharm;

%% Output Structure
Obj_SD_AmideI = SD_Copy(obj_SD);

Obj_SD_AmideI.LocCenter = AmideICenter;
Obj_SD_AmideI.LocMu     = mu_Sim;
Obj_SD_AmideI.LocAlpha  = alpha; % raman tensor vector form [N x 9]

Obj_SD_AmideI.LocFreq   = AmideIFreq;
Obj_SD_AmideI.LocAnharm = AmideIAnharm;

Obj_SD_AmideI.Extra.AmideIAtomSerNo = AmideIAtomSerNo;

