% evolNDtest.m
% If you use it, please cite the corresponding paper:
% Carlota M. Cuesta, Francisco de la Hoz,
% A non-recursive Schur-Decomposition Algorithm for $N$-Dimensional Matrix Equations,
% arXiv:2412.15840, (2024).
%
clear
t=0.1; % final time
NN=[2 3 4 5 6 7 8]; % size of $\mathbf X$
N=length(NN); % $N$
AA=cell(1,length(NN)); % matrices $\mathbf A_j$
for m=1:length(NN)
    AA{m}=rand(NN(m))+1i*rand(NN(m));
end
B=rand(NN)+1i*rand(NN); % $\mathbf B$
X0=rand(NN)+1i*rand(NN); % initial data $\mathbf X_0$
tic
X_Sylv=evolND(AA,B,X0,t); % compute the solution using Theorem 1
toc,tic % elapsed time
mmax=4000; % number of time steps
dt=t/mmax; % $\Delta t$
X_RK=X0; % compute the solution using a fourth-order Runge-Kutta
for m=1:mmax
    k1X=B;
    for j=1:N
        k1X=k1X+multND(AA{j},X_RK,j);
    end
    Xaux=X_RK+.5*dt*k1X;
    k2X=B;
    for j=1:N
        k2X=k2X+multND(AA{j},Xaux,j);
    end
    Xaux=X_RK+.5*dt*k2X;
    k3X=B;
    for j=1:N
        k3X=k3X+multND(AA{j},Xaux,j);
    end
    Xaux=X_RK+dt*k3X;
    k4X=B;
    for j=1:N
        k4X=k4X+multND(AA{j},Xaux,j);
    end
    X_RK=X_RK+dt*(k1X+2*k2X+2*k3X+k4X)/6;
end
toc % elapsed time
disp(norm(X_RK(:)-X_Sylv(:),inf)) % maximum discrepancy