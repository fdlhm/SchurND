% sylvesterNDtest.m
% If you use it, please cite the corresponding paper:
% Carlota M. Cuesta, Francisco de la Hoz,
% A non-recursive Schur-Decomposition Algorithm for $N$-Dimensional Matrix Equations,
% arXiv:2412.15840, (2024).
%
clear
rng(1)
NN=[2 9 33 74 231]; % size of X
N=length(NN); % Number of dimensions
NN=[NN,ones(1,2-N)]; % The solver also works in 1D
AA=cell(1,N); % matrices $\mathbf A_j$
for m=1:N
    AA{m}=randn(NN(m))+1i*randn(NN(m));
end
X=randn(NN)+1i*randn(NN); % solution $\mathbf X$
B=zeros(NN); % right-hand side $\mathbf B$
for m=1:N
    B=B+multND(AA{m},X,m);
end
% numerical approximation using sylvesterND, and elapsed time
tic, Xnum1=sylvesterND(AA,B); toc
% numerical approximation using laplace_merge, and elapsed time
tic, Xnum2=laplace_merge(AA,B); toc
% numerical approximation using laplace_recursive, and elapsed time
tic, Xnum3=laplace_recursive(AA,B); toc
% maximum discrepancies between the exact solution
% and the numerical approximations
disp(max(abs(Xnum1(:)-X(:))))
disp(max(abs(Xnum2(:)-X(:))))
disp(max(abs(Xnum3(:)-X(:))))
% maximum discrepancies between the the numerical approximations
disp(max(abs(Xnum1(:)-Xnum2(:))))
disp(max(abs(Xnum1(:)-Xnum3(:))))
disp(max(abs(Xnum2(:)-Xnum3(:))))
