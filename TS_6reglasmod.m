%TAKAGI-SUGENO RULES 
%Fuzzy sets form: Triangle
%Fuzzified Variables: COJ_X2 and incQCH4 
%Rules number: 6 (3x2)
%Logic operation: and
%Date: 29august2002

function [action]=Rules(binc,D,COJ_X2,incQCH4,pH,De)

%global C1 C2 C3 DELTA1 DELTA2 %DELTA3

C1=1.6;%1.5;%=1.65;%1.7;%1.4;%%1.4;
C2=1.7;%=1.75;%2.3;%1.7;%%1.98;%1.7;
C3=1.8;%=1.85;%1.7;%1.8;%1.99;%5.5;
%C4=3.1;%3.5;%20;

%DELTAS values (Limits of deltaQCH4 sets)
DELTA1=0.004;%.015;%.0005;%%3E-4;%40;
DELTA2=0.009;%.1;%.001;%%3.6E-4;%60;

md=[]; mb=[];
m=zeros(1,4);

%C´s values  (Limits of COJ_X2 sets)
%C1=2;
%C2=6;
%C3=10;
%C4=20;

%DELTAS values (Limits of deltaQCH4 sets)
%DELTA1=70;
%DELTA2=130;
%DELTA3=200;

%RULE 1  (BOUCLE OUVERTE)
%Si COJ/X2 est basse et DQCH4 est bas alors u=0; (BO)
if (0<=COJ_X2 & COJ_X2<C1) & (0<=incQCH4 & incQCH4<=DELTA1)
   mC=memCOJ_X2_3E(C1,C2,C3,COJ_X2);
   mQ=memQCH4_2E(DELTA1,DELTA2,incQCH4);   
   m(1)=min(mC,mQ);       %And rule
else
   m(1)=0;
end;
% 
% %RULE 2  (BOUCLE OUVERTE)
% %Si COJ/X2 est basse et DQCH4 est haut alors u=0; (BO)
if (0<=COJ_X2 & COJ_X2<=C2) & (incQCH4>DELTA1)
   mC=memCOJ_X2_3E(C1,C2,C3,COJ_X2);
   mQ=memQCH4_2E(DELTA1,DELTA2,incQCH4);   
   m(2)=min(mC,mQ);       %And rule
else
   m(2)=0;
end;

%RULE 3 (BOUCLE FERMEE)
%Si COJ/X2 est moyenne et DQCH4 est bas alors u=binc; (BF)
if (C1<=COJ_X2 & COJ_X2<=C3) & (0<=incQCH4 & incQCH4<=DELTA1)
   mC=memCOJ_X2_3E(C1,C2,C3,COJ_X2);
   mQ=memQCH4_2E(DELTA1,DELTA2,incQCH4);   
   m(1)=min(mC,mQ);       %And rule
else
   m(1)=0;
end;

%RULE 4 (BOUCLE FERMEE)
%Si COJ/X2 est moyenne et DQCH4 est haut alors u=binc; (BF)
if (C1<=COJ_X2 & COJ_X2<=C3) & (incQCH4>DELTA1)
   mC=memCOJ_X2_3E(C1,C2,C3,COJ_X2);
   mQ=memQCH4_2E(DELTA1,DELTA2,incQCH4);   
   m(2)=min(mC,mQ);       %And rule
else
   m(2)=0;
end;

%RULE 5 (BOUCLE FERMEE)
%Si COJ/X2 est haute et DQCH4 est bas alors u=D; (BF)
if (C3<COJ_X2) & (0<=incQCH4 & incQCH4<=DELTA1)
   mC=memCOJ_X2_3E(C1,C2,C3,COJ_X2);
   mQ=memQCH4_2E(DELTA1,DELTA2,incQCH4);   
   m(3)=min(mC,mQ);       %And rule
else
   m(3)=0;
end;

%RULE 6 (BOULCE OUVERTE)
%Si COJ/X2 est haute et DQCH4 est haut alors u=0; (BO)
if (C3<=COJ_X2) & (incQCH4>DELTA1)
   mC=memCOJ_X2_3E(C1,C2,C3,COJ_X2);
   mQ=memQCH4_2E(DELTA1,DELTA2,incQCH4);   
   m(4)=min(mC,mQ);       %And rule
else
   m(4)=0;
end;

m_binc=m(1)+m(2);
m_D=m(3)+m(4);

%Defuzzyfication
if m_D ~= 0
   D_action=m_D*D/sum(m);
else
   D_action=De;
end;

if m_binc ~= 0
   binc_action=m_binc*binc/sum(m);
  else
   binc_action=0;
end;

md=[md;D_action];
mb=[mb;binc_action];
%ol_action=0;

action=[1*binc_action .48*D_action];
%action=[.8*binc_action 1*D_action]; %.48
