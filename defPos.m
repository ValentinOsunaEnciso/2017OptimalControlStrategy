% Test: Matrix Definite Positive, Valentín Osuna-Enciso, CUTONALA, 2014.
% A simetric matrix of size (n x n) is definite positive if:
% I)    All the upper left determinants of size (k x k) are positive, or
%       det(A_k)>0, for all 1<=k<=n.
function DP=defPos(A)
    n=size(A,1);
    signos=zeros(1,n);
    signos(1,1)=A(1,1);
    for ind=2:n
        signos(1,ind)=det(A(1:ind,1:ind));
    end
    if(signos>0)
        DP=true;
    else
        DP=false;
    end
end