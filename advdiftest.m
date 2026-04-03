% advdiftest.m
% If you use it, please cite the corresponding paper:
% Carlota M. Cuesta, Francisco de la Hoz,
% A non-recursive Schur-Decomposition Algorithm for $N$-Dimensional Matrix Equations,
% arXiv:2412.15840, (2024).
%
clear
tic
t=1; % time at which the numerical solution will be computed
N=6; % number of dimensions $N$
M=16; % number of nodes at each dimension
[x,DD]=herdif(M,2,1.4); % create $\mathbf x$ and the differentiation matrices
D1=DD(:,:,1); % first-order differentiation matrix
D2=DD(:,:,2); % second-order differentiation matrix
AA=cell(1,N); % stores the matrices $\mathbf A_j$
XX=cell(1,N); % stores the spatial mesh
[XX{:}]=ndgrid(x); % creates the spatial mesh
B=-exp(-XX{1}.^2); % creates $\mathbf B$
for m=2:N
    B=B.*exp(-XX{m}.^2);
end
clear XX % we do not need XX any more
U0=-2*B; % initial data
Uexact=-(1+exp(t))*B; % exact solution at time $t$
A=D2+2*diag(x)*D1+((2*N+1)/N)*eye(M); % creates $\mathbf A_j$
for m=1:N
    AA{m}=A;
end
U=real(evolND(AA,B,U0,t)); % computes the solution at time $t$
norm(U(:)-Uexact(:),inf) % error
toc % elapsed time