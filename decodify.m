% Decodify bitstrings, CUTONALA, José Valentín Osuna Enciso, 2014.
function X = decodify(B,bit,d,dec,r)
    % X		-> real value (precision: 6)
    % B		-> binary string (length: bit)
    bits=bit-1; s=size(B,1); aux=0:1:bits; X=[];
    for ind=0:d-1
        temp1 = fliplr(B(:,(bit*ind)+1:bit*(ind+1)));                  
        temp2 = ones(s(1),1)*aux;
        x1 = sum((temp1.*2.^temp2)');
        X = [X; r(2,ind+1)-x1.*((r(2,ind+1)-r(1,ind+1))/(dec-1));];
    end
end