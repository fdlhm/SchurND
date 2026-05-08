% compute $\mathbf A\square_j\mathbf X$ in older MATLAB versions
function X=multNDalt(A,X,j)
if isscalar(A) % test whether A is a scalar
    X=A*X; % trivial case
elseif j==1 % multiplication along the 1st dimension is simpler
    sizeX=size(X); % size of X
    X=reshape(X,sizeX(j),[]); % reshape X
    X=A*X; % compute the product
    X=reshape(X,sizeX); % reshape the product
else
    sizeX=size(X); % size of X
    N=length(sizeX); % number of dimensions
    permutation=[j 2:j-1 1 j+1:N];
    X=permute(X,permutation); % swap the first and jth dimensions
    sizeXpermuted=size(X); % size of X, after permutation
    X=reshape(X,sizeX(j),[]); % reshape X
    X=A*X; % compute the product
    X=reshape(X,sizeXpermuted); % reshape the product
    X=permute(X,permutation); % swap the first and jth dimensions
end