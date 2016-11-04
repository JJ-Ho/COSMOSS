function [Xout,Yout,Zout]=gplot3(A,xyz,varargin)

%   John Gilbert, 1991.
%   Modified 1-21-91, LS; 2-28-92, 6-16-92 CBM.
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.12 $  $Date: 2002/04/15 04:13:43 $

%   modified by Johnson Chuang 2006

%% Debug
% A =Conn1;
% xyz = Carbon_Pos;

%% Main
xyz = double(xyz);

[i,j] = find(A);
[ignore, p] = sort(max(i,j));
i = i(p);
j = j(p);

% Create a long, NaN-separated list of line segments,
% rather than individual segments.

X = [ xyz(i,1) xyz(j,1) repmat(NaN,size(i))]';
Y = [ xyz(i,2) xyz(j,2) repmat(NaN,size(i))]';
Z = [ xyz(i,3) xyz(j,3) repmat(NaN,size(i))]';
X = X(:);
Y = Y(:);
Z = Z(:);

if nargout==0,
    if nargin<3,
        plot3(X, Y, Z)
    else
        plot3(X, Y, Z, varargin{:});
    end
else
    Xout = X;
    Yout = Y;
    Zout = Z;
end
