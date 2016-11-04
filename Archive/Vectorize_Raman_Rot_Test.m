%% Tried to test if vectorize way for Raman tensor ratation would be faster
% The for loop is way faster!

Mode_Num = 100;
Raman = rand(3,3,Mode_Num);
Rot = rand(3);

%% For loop
tic
Raman_Rot = zeros(size(Raman));
for i = 1:Mode_Num
    Raman_Rot(:,:,i) = Rot*squeeze(Raman(:,:,i))*Rot';
end
toc

%% vectorized way to do apply rotational matrix on Raman Matrix [3,3,N]
tic
bsxfun(@times,Rot,ones(3,3,Mode_Num));
C_Rot = num2cell(bsxfun(@times,Rot,ones(3,3,Mode_Num)),[1,2]);
Rot_blkdiag = blkdiag(C_Rot{:});

bsxfun(@times,Rot',ones(3,3,Mode_Num));
C_Rot_inv = num2cell(bsxfun(@times,Rot',ones(3,3,Mode_Num)),[1,2]);
Rot_inv_blkdiag=blkdiag(C_Rot_inv{:});

C_Raman = num2cell(Raman,[1,2]);
Raman_blkdiag = blkdiag(C_Raman{:});

Raman_Rot_blkdiag = Rot_blkdiag * Raman_blkdiag * Rot_inv_blkdiag;

% inverse of blkdiag
[I,~,J] = ndgrid(1:3,1,1:3);
K = double(3*(0:Mode_Num-1)); % if use 3, K will be int32 and the following bsxfun will be confused...
I = bsxfun(@plus,I,K);
J = bsxfun(@plus,J,K);

Raman_Rot = Raman_Rot_blkdiag(sub2ind(size(Raman_Rot_blkdiag),I,J));
Raman_Rot = reshape(Raman_Rot,3,3,[]);
toc

