function out=rotate_nma(fpar,nma_s,nma_m)

[a,b] = size(nma_s);

T = [fpar(1) fpar(2) fpar(3)];

Rx = [1 0 0; 0 cos(fpar(4)) -sin(fpar(4)); 0 sin(fpar(4)) cos(fpar(4))];
Ry = [cos(fpar(5)) 0 sin(fpar(5)); 0 1 0; -sin(fpar(5)) 0 cos(fpar(5))];
Rz = [cos(fpar(6)) -sin(fpar(6)) 0; sin(fpar(6)) cos(fpar(6)) 0; 0 0 1];

for n = 1:a
    rot(n,:) = (Rx*Ry*Rz*(nma_s(n,:)' + T'))';
end

out = sum(sum((nma_m - rot).^2));