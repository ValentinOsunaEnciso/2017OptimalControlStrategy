function f = FUNCION_CSA(X,Np,d)
    f=zeros(Np,1);
    for ind=1:Np
        f(ind,1)=1*X(ind,1).*sin(4*pi.*X(ind,1))-1*X(ind,2).*sin(4*pi.*X(ind,2)+pi)+1;
    end
end