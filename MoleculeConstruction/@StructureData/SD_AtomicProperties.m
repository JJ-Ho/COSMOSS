function AP = SD_AtomicProperties(obj_SD)
%% Load Atom Properties
% The 'AtomicProperties.mat' was extracted from Mathematica with cmd:
% ´MatrixForm[Table[{ElementData[z, "AtomicNumber"], ElementData[z, "Abbreviation"], ElementData[z, "AtomicRadius"], ElementData[z, "AtomicWeight"]}, {z, 118}]]
% and with some string clean-up using vim.
% load : AtomName
%        AtomMass,[amu]
%        AtomNumber 
%        AtomRadii,[pm]
load('AtomicProperties.mat')

%% Match Atom Name with its properties
NAtoms = obj_SD.NAtoms;
Name   =  cell(NAtoms,1);
Mass   = zeros(NAtoms,1);
Number = zeros(NAtoms,1);
Radii  = zeros(NAtoms,1);

for i = 1:NAtoms
    % this extra step taking care of PDB atom names, which have 
    % multiple charaters, eg: CA => C atom
    AtomName_Str = obj_SD.AtomName{i};
    if ~any(strcmp(AtomName_Str,AtomName))
        AtomName_Str = AtomName_Str(1);
    end
    AtomIndex = strcmp(AtomName_Str,AtomName);

    Name(i)   = AtomName(AtomIndex);
    Mass(i)   = AtomMass(AtomIndex);
    Number(i) = AtomNumber(AtomIndex);
    Radii(i)  = AtomRadii(AtomIndex);
end

AP.Name   = Name;
AP.Mass   = Mass;
AP.Number = Number;
AP.Radii  = Radii;