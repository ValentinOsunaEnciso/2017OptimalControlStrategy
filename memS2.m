%In this program the member level is calculated
%For the COJ_X2 variable.

%Fuzzy sets form: triangle.
%Last modification: 29august2002

function ms=memS2(x4in,X2pt5,S2pt5,ICpt5,QCH4pt5,X2_1,S2_1,IC_1,QCH4_1,X2,S2,IC,QCH41,X2k2,x4r,x5r,QCH4k2, X2k,S2r,ICr,QCH427)


LS1=.07;
LS2=.11;
LS3=.14;
LS4=.175;
LS5=.21;
LS6=.26;



v_ms=zeros(1,5);

% FIRST SET   'BASSE'
  if LS1<=x4in & x4in<LS2        
     v_ms(1)=1;%COJ_X2/C1;
  elseif LS2<=x4in & x4in<LS3
     v_ms(1)=(LS3-x4in)/(LS3-LS2);
  else
     v_ms(1)=0;
  end;

% SECOND SET   'MOYENNE'
  if LS2<=x4in & x4in<LS3
     v_ms(2)=(x4in-LS2)/(LS3-LS2);
  elseif LS3<=x4in & x4in<LS4
     v_ms(2)=(LS4-x4in)/(LS4-LS3);
  else
     v_ms(2)=0;         
  end;

% THIRD SET   'HAUTE'
if LS3<=x4in & x4in<LS4
     v_ms(3)=(x4in-LS3)/(LS4-LS3);
  elseif LS4<=x4in & x4in<LS5
     v_ms(3)=(LS5-x4in)/(LS5-LS4);
  else
     v_ms(3)=0;         
  end;


% FOURTH SET
if LS4<=x4in & x4in<LS5
     v_ms(4)=(x4in-LS4)/(LS5-LS4);
  elseif LS5<=x4in & x4in<LS6
     v_ms(4)=(LS6-x4in)/(LS6-LS5);
  else
     v_ms(4)=0;         
  end;

% FIFTH SET
  if LS5<=x4in & x4in<LS6
     v_ms(5)=(x4in-LS5)/(LS6-LS5);
  elseif LS6<=x4in
     v_ms(5)=1;%(C4-COJ_X2)/(C4-C3);
  else
     v_ms(5)=0;
  end;
  
  X2ref=[X2pt5 X2_1 X2 X2k2 X2k]';
  S2ref=[S2_1 S2_1 S2 x4r S2r]'; %se repite S2_1 porque la primera parte (t<200 h)de la trayectoria es menor comparada con las demas
  ICref=[IC_1 IC_1 IC x5r ICr]';
QCH4ref=[QCH4pt5 QCH4_1 QCH41 QCH4k2 QCH427]';
       kk=[7 6 3.5 2.6 2.5]';
%      kk=[7 6 2.58 2.62]';
   
%jj=v_ms
  
ms=v_ms*[X2ref S2ref ICref QCH4ref kk];
