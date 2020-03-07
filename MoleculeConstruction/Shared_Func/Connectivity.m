function Conn = Connectivity(Atom,XYZ,varargin)
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

%% Input parser
N = size(XYZ,1);

if nargin > 2
    INPUT = inputParser;
    INPUT.KeepUnmatched = 1;

    % Default values
    defaultBondLength = 1.6;

    % Add optional inputs to inputparser object
    addOptional(INPUT,'BondLength',defaultBondLength);

    parse(INPUT,varargin{:});

    BondLength = INPUT.Results.BondLength;
    RC = ones(N,1).* (BondLength/2*100); % turn unit into pm
    BL_Scale = 1.1; 
else
    %% Load covalent radii data
    % data from Mathematica, cmd:
    % A = Table[{ElementData[z, "Abbreviation"], ElementData[z, "CovalentRadius"]}, {z, 118}]
    % Only 96 elements have data
    load('AtomProperties.mat')

    RC = zeros(N,1); % [Nx1]

    for i = 1:N
        STR = Atom{i}; % this extra step if for PDB atom names which have multiple charaters, ex: CA => C atom
        RC(i) = Radii(strcmp(STR(1),AtomName));
    end
    BL_Scale = 1.5;
end

%% generate connectivity matrix ndgrid version
% Define bonded if dist. < (radii A + radii B)* BL_Scale

% generate index for substracting
[n m] = ndgrid(1:N);
% Generate bond length matrix, /100 frm pm to Ang, *1.1 scaling facotr
B = reshape((RC(m(:)) + RC(n(:)) )./100.*BL_Scale,N,N);
% substract atom's position according to index generate above
T1 = XYZ(m(:),:)-XYZ(n(:),:);
% Get atom distance
T2 = sqrt(sum(T1.^2,2));
% reshape a (anum^2)x1 matrix to anum x anum matrix
Distance_matrix = reshape(T2,N,N);
% Check Bondlength
Bonded = Distance_matrix < B;
% grep lower triangular part in order to avoid over counting
lower = tril(Bonded,-1);
% Transform into bonding index metrix
[a b] = find(lower>0);
C_index = [a,b];
Conn = false(N);
Conn(C_index(:,1)+(C_index(:,2)-1)*N) = true;
Conn = Conn|Conn';

