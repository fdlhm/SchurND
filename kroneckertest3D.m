% kroneckertest3D.m
% If you use it, please cite the corresponding paper:
% Carlota M. Cuesta, Francisco de la Hoz,
% A non-recursive Schur-Decomposition Algorithm for $N$-Dimensional Matrix Equations,
% arXiv:2412.15840, (2026).
clear
rng(1)
NN=[36,36,36]; % size of X
A1=rand(NN(1))+1i*rand(NN(1)); % $\mathbf A_1$
A2=rand(NN(2))+1i*rand(NN(2)); % $\mathbf A_2$
A3=rand(NN(3))+1i*rand(NN(3)); % $\mathbf A_3$
X=rand(NN)+1i*rand(NN); % solution $\mathbf X$
B=multND(A1,X,1)+multND(A2,X,2)+multND(A3,X,3); % right-hand side $\mathbf B$
A123=kron(kron(A3,eye(NN(2))),eye(NN(1)))+kron(kron(eye(NN(3)),A2),eye(NN(1)))+kron(kron(eye(NN(3)),eye(NN(2))),A1); % $\mathbf A_3\oplus \mathbf A_2\oplus\mathbf A_1$
% numerical approximation using Kronecker sums, and elapsed time
tic, Xnum1=reshape(A123\B(:),NN); toc
% numerical approximation using sylvesterND, and elapsed time
tic, Xnum2=sylvesterND({A1,A2,A3},B); toc
norm(X(:)-Xnum1(:),inf)
norm(X(:) - Xnum2(:),inf)