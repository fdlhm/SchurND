% sylvesterND.m
% If you use it, please cite the corresponding paper:
% Carlota M. Cuesta, Francisco de la Hoz,
% A non-recursive Schur-Decomposition Algorithm for $N$-Dimensional Matrix Equations,
% arXiv:2412.15840, (2024).
%
function X=sylvesterND(AA,X)
N=length(AA);
sizeX=size(X);
M=min(N,length(sizeX));
UU=cell(1,M);
TT=cell(1,N);
NN=zeros(1,N);
for j=1:N
    NN(j)=size(AA{j},1);
end
cumprodN=cumprod([1 NN(1:M)]);
cumprodNinv=cumprod([1 NN(M:-1:2)]);
for j=1:M
    if NN(j)==1
        TT{j}=AA{j};
    else
        [UU{j},TT{j}]=schur(AA{j},'complex');
        X=reshape(X,cumprodN(j),NN(j),cumprodNinv(M+1-j));
        X=pagemtimes(X,conj(UU{j}));
    end
    AA{j}=[];
end
for j=M+1:N
    TT{j}=AA{j};
    AA{j}=[];    
end
clear AA
ii=NN;
for index=cumprodN(M+1):-1:1
    num=X(index);
    den=0;
    for j=1:N
        Tj=TT{j};
        den=den+Tj(ii(j),ii(j));
        for k=1:NN(j)-ii(j)
            num=num-X(index+k*cumprodN(j))*Tj(ii(j),k+ii(j));
        end
    end
    X(index)=num/den;
    k=1;
    while k<=M
        if ii(k)>1
            ii(k)=ii(k)-1;
            break;
        end
        ii(k)=NN(k);
        k=k+1;
    end
end
clear TT Tj
for j=1:M
    if NN(j)>1
        X=reshape(X,cumprodN(j),NN(j),cumprodNinv(M+1-j));
        X=pagemtimes(X,'none',UU{j},'transpose');
        UU{j}=[];
    end
end
X=reshape(X,sizeX);
