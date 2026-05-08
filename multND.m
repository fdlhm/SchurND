% multND.m
% If you use it, please cite the corresponding paper:
% Carlota M. Cuesta, Francisco de la Hoz,
% A non-recursive Schur-Decomposition Algorithm for $N$-Dimensional Matrix Equations,
% arXiv:2412.15840, (2026).
%
% compute $\mathbf A\square_j\mathbf X$ in MATLAB R2020b or newer
function X=multND(A,X,j)
if isscalar(A) % test whether A is a scalar
    X=A*X; % trivial case
else
    sizeX=size(X);
    X=reshape(X,prod(sizeX(1:j-1)),sizeX(j),prod(sizeX(j+1:end)));
    X=pagemtimes(X,A.');
    X=reshape(X,sizeX);
end
