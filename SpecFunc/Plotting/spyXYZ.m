function spyXYZ(hAx,X,Y,Z)
% spyXY Visualize sparsity pattern with X Y inputs.
% This code is modified from the original spy by JJH
% Copyright 1984-2013 The MathWorks, Inc. 
[m,n] = size(Z);
markersize = 0; 
if markersize == 0
   units = get(gca,'units');
   set(gca,'units','points');
   pos = get(gca,'position');
   markersize = max(4,min(14,round(6*min(pos(3:4))/max(m+1,n+1))));
   set(gca,'units',units);
end

[I,J,V] = find(Z);
MS = abs(V)./max(abs(V)).*markersize*5;

X_j = X(J);
Y_i = Y(I);

scatter(hAx,X_j,Y_i,MS,V,'filled')

