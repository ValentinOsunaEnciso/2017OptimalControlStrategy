%In this program the member level is calculated
%For the COJ_X2 variable.

%Fuzzy sets form: triangle.
%Last modification: 29august2002

function mC=memCOJ_X2(C1,C2,C3,COJ_X2)

v_mC=zeros(1,3);

% FIRST SET   'BASSE'
  if 0<=COJ_X2 & COJ_X2<C1        
     v_mC(1)=1;%COJ_X2/C1;
  elseif C1<=COJ_X2 & COJ_X2<C2
     v_mC(1)=(C2-COJ_X2)/(C2-C1);
  else
     v_mC(1)=0;
  end;

% SECOND SET   'MOYENNE'
  if C1<=COJ_X2 & COJ_X2<C2
     v_mC(2)=(COJ_X2-C1)/(C2-C1);
  elseif C2<=COJ_X2 & COJ_X2<C3
     v_mC(2)=(C3-COJ_X2)/(C3-C2);
  else
     v_mC(2)=0;         
  end;

% THIRD SET   'HAUTE'
  if C2<=COJ_X2 & COJ_X2<C3
     v_mC(3)=(COJ_X2-C2)/(C3-C2);
  elseif C3<=COJ_X2
     v_mC(3)=1;%(C4-COJ_X2)/(C4-C3);
  else
     v_mC(3)=0;
  end;
  
mC=max(v_mC);      