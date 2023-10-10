
%global C1 C2 C3 DELTA1 DELTA2%k
%k=0;
%C´s values  (Limits of COJ_X2 sets)
C1=1.6;%1.65;%1.5;%1.7;%1.4;%%1.4;
C2=1.7;%1.75;%1.6;%2.3;%1.7;%%1.98;%1.7;
C3=1.8;%1.85;%1.8;%1.7;%1.8;%1.99;%5.5;
%C4=3.1;%3.5;%20;

%DELTAS values (Limits of deltaQCH4 sets)
DELTA1=.004;%0.004;%.0005;%%3E-4;%40;
DELTA2=.009;%0.009;%.001;%%3.6E-4;%60;
%DELTA3=7E-4;%200;
%DELTA4=200;


%Values  for binc_k action
Bset=Bic_0;            %%%  Valor deseado B*
binc_min=-0.9*Bset;%-0.45*Bset;
binc_max=12*Bset;
%binc_sat=-0.35*Bset;

alpha=20;%10;
delta=0.05;%0.017;
K2c=(1/(1+alpha))*(log((binc_max-binc_min)/-binc_min)/log(Bset/(Bset-delta)));
K1c=alpha*K2c;

%Values  for D_k action
Deq=D_0;
Dmax=1.5*Deq;
Dsat=1000;
D_anti=10E-4;

alphad=100;%0.01;
deltad=0.1;%0.08;

K2d=(1/(1+alphad))*(log(Dmax/Deq)/log(Bset/(Bset+deltad)));
K1d=alphad*K2d;

