function [Phi,Psi] = PeptideDihedral(Structure)
%% Debug
% Structure = Data_COSMOSS.Structure;

%% Prep parameters
Ind_CONCa = Structure.Extra.AmideIAtomSerNo; %[C,O,N,CA]
XYZ       = Structure.XYZ;

%% prep bond vectors
%     O(i)                     O(i+1)
%     ||                       ||
% ... C(i) -> N(i) -> Ca(i) -> C(i+1) -> N(i+1) - ...
%          V1      V2       V3        V4

V1 = XYZ(Ind_CONCa(    :,3),:) - XYZ(Ind_CONCa(      :,1),:); % Vec:   N(i) <-   C(i)
V2 = XYZ(Ind_CONCa(    :,4),:) - XYZ(Ind_CONCa(      :,3),:); % Vec:  Ca(i) <-   N(i)
V3 = XYZ(Ind_CONCa(2:end,1),:) - XYZ(Ind_CONCa(1:end-1,4),:); % Vec: C(i+1) <-  Ca(i)
V4 = XYZ(Ind_CONCa(    :,3),:) - XYZ(Ind_CONCa(      :,1),:); % Vec: N(i+1) <- C(i+1)

% N mode -> N-1 dihedrals 
V1(end,:) = [];
V2(end,:) = [];
V4(  1,:) = []; % shift index to the next amide group

%% normalize all vectors
V1 = bsxfun(@rdivide,V1,sqrt(sum(V1.^2,2)));
V2 = bsxfun(@rdivide,V2,sqrt(sum(V2.^2,2)));
V3 = bsxfun(@rdivide,V3,sqrt(sum(V3.^2,2)));
V4 = bsxfun(@rdivide,V4,sqrt(sum(V4.^2,2)));

%% Phi angle
N_12 = cross(V1,V2,2);
N_23 = cross(V2,V3,2);

X1 = dot(cross(N_12,N_23,2),V2,2);
Y1 = dot(N_12,N_23,2);
Phi = atan2d(X1,Y1);

%% Psi angle
N_34 = cross(V3,V4,2);

X2 = dot(cross(N_23,N_34,2),V3,2);
Y2 = dot(N_23,N_34,2);
Psi = atan2d(X2,Y2);

%% Debug plot
% ramachandran('2lmq.pdb')
% hF = gcf;
% hAx = hF.Children;
% 
% hold(hAx,'on')
%     scatter(hAx,Phi,Psi,12,'r','filled')
% hold(hAx,'off')



