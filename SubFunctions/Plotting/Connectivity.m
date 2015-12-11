function Conn = Connectivity(XYZ,varargin)
%% Connectivity
%  
% This script determines atomic connectivities from gamess xyz output
% 
% ------- Version log -----------------------------------------------------
%
% Ver. 1.0  140605  Modified from Molecule_Plot.m. Start version log
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2014

%% debug
% clear all
% Read input.geo

%% Input parser

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultBondLength = 1.6;

% Add optional inputs to inputparser object
addOptional(INPUT,'BondLength',defaultBondLength);

parse(INPUT,varargin{:});

BondLength = INPUT.Results.BondLength  ;

%% generate connectivity matrix ndgrid version

A_Num=size(XYZ,1);

% generate index for substracting
[n m] = ndgrid(1:A_Num);
% substract atom's position according to index generate above
T1 = XYZ(m(:),:)-XYZ(n(:),:);
% Get atom distance
T2 = sqrt(sum(T1.^2,2));
% reshape a (anum^2)x1 matrix to anum x anum matrix
Distance_matrix = reshape(T2,A_Num,A_Num);
% grep lower triangular part in order to avoid over counting
lower = tril(Distance_matrix,-1);
% Transform into bonding index metrix
[a b] = find(lower< BondLength & lower>0);
C_index = [a,b];
Conn = false(A_Num);
Conn(C_index(:,1)+(C_index(:,2)-1)*A_Num) = true;
Conn = Conn|Conn';

