% sylvesterNDtest.m
% If you use it, please cite the corresponding paper:
% Carlota M. Cuesta, Francisco de la Hoz,
% A non-recursive Schur-Decomposition Algorithm for $N$-Dimensional Matrix Equations,
% arXiv:2412.15840, (2024).
%
clear
NN=[2 9 33 74 231]; % size of X
AA=cell(1,length(NN)); % matrices $\mathbf A_j$
for m=1:length(NN)
    AA{m}=rand(NN(m))+1i*rand(NN(m));
end
X=rand(NN)+1i*rand(NN); % solution $\mathbf X$
B=zeros(NN); % right-hand side $\mathbf B$
for m=1:length(NN)
    B=B+multND(AA{m},X,m);
end
tic
Xnum=sylvesterND(AA,B); % computes $\mathbf X_{num}$
toc % elapsed time
disp(max(abs(Xnum(:)-X(:)))) % maximum discrepancy