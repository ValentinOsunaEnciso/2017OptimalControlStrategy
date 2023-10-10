% Función para calcular la COJ/X2

function COJ_X2=cojX2(D,S2ini,A2,X2,t_sim,T_bgn)


% t_sim=15000;
% T_bgn=200;
% D=.1;
% S2ini=.1;
% A2=.1;
% X2=.1;

if t_sim<=T_bgn             %Before perturbation
   COJ_X2=0;
elseif t_sim>T_bgn
   COJ_X2=(D*S2ini*A2)/(X2);  %After perturbation
else
   COJ_X2=0 ;
end;