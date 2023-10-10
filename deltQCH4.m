%Función para calcular el incremento de metano (%)

function deltQCH4=incremento1(times1,QCH4n1,S2eq,S2pert1,R1,R2,R3,t_deb)

%persistent QCH4eq
%persistent deltaS2in
%persistent deltaQCH4total

if times1<=t_deb           % Avant la perturbation
   %QCH4eq=QCH4n;
   deltaS2in=S2pert1-S2eq; %constante
   deltQCH4=R1*R2*deltaS2in/R3; %constante (100%)
   %deltaQCH4=QCH4n-QCH4eq; %variable
else
   deltaS2in=S2pert1-S2eq; %constante
   deltQCH4=R1*R2*deltaS2in/R3;
   %incQCH4=deltaQCH4actuel*100/deltaQCH4total;
   %deltQCH4=deltaQCH4;
   %incQCH4=deltaQCH4actuel*deltaQCH4total/100;
   %incQCH4=deltaQCH4total*100/deltaQCH4actuel;
end;
   