%In this program the member level is calculated
%For the incQCH4 variable.

%Fuzzy sets form: triangle.
%Last modification: 29august2002

function mQ=memQCH4_2E(DELTA1,DELTA2,incQCH4)

v_mQ=zeros(1,2);

% FIRST SET   'BAS'
  if 0<=incQCH4 & incQCH4<DELTA1        
     v_mQ(1)=1;%incQ/DELTA1;
  elseif DELTA1<=incQCH4 & incQCH4<DELTA2
     v_mQ(1)=(DELTA2-incQCH4)/(DELTA2-DELTA1);
  else
     v_mQ(1)=0;
  end;

% SECOND SET   'HAUT'
  if DELTA1<=incQCH4 & incQCH4<DELTA2
     v_mQ(2)=(incQCH4-DELTA1)/(DELTA2-DELTA1);
  elseif incQCH4>=DELTA2;%DELTA2<incQCH4 & incQCH4<DELTA3   
     v_mQ(2)=1;%(DELTA3-incQ)/(DELTA3-DELTA2);
  else
     v_mQ(2)=0;
  end;
  
mQ=max(v_mQ);      