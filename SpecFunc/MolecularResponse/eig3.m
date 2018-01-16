function [VX,DX] = eig3(X)
% the eig3 can do eig for multiple square matrix at once
% given input X, size M == [M,N,N]. There are M number of N-by-N spare
% matrix, eig3 put these M matrix on the block diagonal and to eig at once.
% The output VX is the eigenvecotr of each square matrix, with size = [M,N,N]
% The output DX is the eignevalues of each square matrix, with size = [M,N]
% The code utilize sparse matrix to save memory

%% Debug
% X = rand(3,3);
% X = X + X';
% X = reshape(X,1,3,3);
% X = repmat(X,7,1,1);

%% check input size 
if ~eq(size(X,2),size(X,3))
    disp('The input has to be a NxNxM, 3D array!')
    return
end

%% Main
M = size(X,1);
N = size(X,2);

VX = zeros(M,N,N);
DX = zeros(M,N);


for i = 1:M
    [V,D] = eig(squeeze(X(i,:,:)));
    VX(i,:,:) = V;
    DX(i,:) = diag(D);
end

%% Vectorized version
%  this doesn't work since the eigen value order was messed up
% %% Main
% N = size(X,1);
% M = size(X,3);
% 
% [I,J]   = ndgrid(1:N);
% Ind     = repmat([I(:),J(:)],1,1,M);
% K       = reshape(0:N:N*(M-1),1,1,M);
% Ind_All = bsxfun(@plus,Ind,K);
% Ind1    = reshape(Ind_All(:,1,:),[],1);
% Ind2    = reshape(Ind_All(:,2,:),[],1);
% 
% %%
% X_S_blkdiag = sparse(Ind1,Ind2,X(:));
% X_S_blkdiag = full(X_S_blkdiag);
% [V,D] = eig(X_S_blkdiag);
% 
% % reorgainze the eigenvector and eigenvalues into designed shape 
% DX = reshape(diag(D),N,M);
% % VX = reshape(V(logical(T)),N,N,M);