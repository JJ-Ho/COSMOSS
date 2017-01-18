function Conn = Connectivity_New(Atom,XYZ)
%% Connectivity
%  
% This script determines atomic connectivities from gamess xyz output
% 
% ------- Version log -----------------------------------------------------
%
% Ver. 2.0  170118  Add covalent radii data so I can determine bonding more 
%                   accurately  
% 
% Ver. 1.0  140605  Modified from Molecule_Plot.m. Start version log
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2014

%% debug
% XYZ = Data_TwoDGrid.Structure.XYZ;
% Atom = Data_TwoDGrid.Structure.G09_Output.Atom_Name;

%% Load covalent radii data
% data from Mathematica, cmd:
% A = Table[{ElementData[z, "Abbreviation"], ElementData[z, "CovalentRadius"]}, {z, 118}]
% Only 96 elements have data
load('AtomProperties.mat')

N_Atom = length(Atom);
RC = zeros(N_Atom,1); % [Nx1]

for i = 1:N_Atom
    RC(i) = Radii(strcmp(Atom(i),AtomName));
end

%% generate connectivity matrix ndgrid version

A_Num=size(XYZ,1);

% generate index for substracting
[n m] = ndgrid(1:A_Num);
% Generate bond length matrix, /100 frm pm to Ang, *1.1 scaling facotr
B = reshape((RC(m(:)) + RC(n(:)) )./100.*1.1,A_Num,A_Num);
% substract atom's position according to index generate above
T1 = XYZ(m(:),:)-XYZ(n(:),:);
% Get atom distance
T2 = sqrt(sum(T1.^2,2));
% reshape a (anum^2)x1 matrix to anum x anum matrix
Distance_matrix = reshape(T2,A_Num,A_Num);
% Check Bondlength
Bonded = Distance_matrix < B;
% grep lower triangular part in order to avoid over counting
lower = tril(Bonded,-1);
% Transform into bonding index metrix
[a b] = find(lower>0);
C_index = [a,b];
Conn = false(A_Num);
Conn(C_index(:,1)+(C_index(:,2)-1)*A_Num) = true;
Conn = Conn|Conn';

