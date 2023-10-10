%% Differential Evolution, best/1/bin, Valentin Osuna-Enciso, CUCEI-UDG,
% CUTONALA, Junio, 2014
function Fitness_best=DE_DA  
  N=15000;
  [r,d]=INICIA_DA; 
  corridas=[];
  for ind5=1:1 %corrida
  %INICIALIZACION DE: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Np=5;          %Tamano de la poblacion.90
  F=0.2;         %Factor de escalamiento .25
  Cr=0.8;         %Probabilidad de cruza .8
  maxIter=5;       %Numero maximo de iteraciones
  k=0;            %Contador de iteraciones
  for ind1=1:Np
      defi1=0; defi2=0;
      while ~(defi1 && defi2)
        x(ind1,:)=r(1,:)+(r(2,:)-r(1,:)).*rand(1,d); 
        Pu1=[x(ind1,1),x(ind1,2),x(ind1,3);x(ind1,2),x(ind1,4),x(ind1,5);x(ind1,3),x(ind1,5),x(ind1,6)];
        Pd2=[x(ind1,7),x(ind1,8),x(ind1,9);x(ind1,8),x(ind1,10),x(ind1,11);x(ind1,9),x(ind1,11),x(ind1,12)];
        defi1=defPos(Pu1); defi2=defPos(Pd2);
      end
      Fitness(ind1)=FUNCION_DA3(Pu1,Pd2,N,1); 
  end
  EvaluacioneFO=Np;
  Fitness_best= Fitness(1); x_best=x(1,:); EvalFO=Np;
  disp(sprintf('Iterac=%d, Fitnes=%e, Eval:%d\n',k,Fitness_best,EvalFO));  
  while k<maxIter 
   for ind1=1:Np
       r1=randi(Np);r2=randi(Np);
       while(r1==r2==ind1)
          r1=randi(Np); r2=randi(Np);       %Generados sean diferentes.
       end
       v(1,:)=x_best+F*(x(r1,:)-x(r2,:));   %Mutacion
       u(1,:)=x(ind1,:);
       j_rand=randi(d);
       for ind2=1:d                         %Genero vector de prueba
          if(rand()<Cr || ind2==j_rand) && v(1,j_rand)>0     %Cruza
             u(1,j_rand)=v(1,j_rand);
          end
       end
       Pu1=[u(1,1),u(1,2),u(1,3);u(1,2),u(1,4),u(1,5);u(1,3),u(1,5),u(1,6)];
       Pd2=[u(1,7),u(1,8),u(1,9);u(1,8),u(1,10),u(1,11);u(1,9),u(1,11),u(1,12)];
       defi1=defPos(Pu1); defi2=defPos(Pd2);
       if (defi1 && defi2) % Si son definidas positivas, evalúa      
        temp=FUNCION_DA3(Pu1,Pd2,N,1);
       else
        temp=Fitness(ind1);   
       end
       if(temp<Fitness(ind1))
          x(ind1,:)=u(1,:);
          Fitness(ind1)=temp;
          if (temp<Fitness_best)
             x_best=u(1,:);
             Fitness_best=temp;
          end
       end
   end
   EvaluacioneFO=EvaluacioneFO+Np;
   k=k+1; EvalFO=EvalFO+Np;
   CONVERGENCIA(k)=Fitness_best;   
   disp(sprintf('Iterac=%d, Fitnes=%e, Eval:%d\n',k,Fitness_best,EvaluacioneFO));   
  end
  %% Muestra resultados de la evolucion de DE:
  semilogy(CONVERGENCIA,'r','LineWidth',3),figure,  
  t=[1:N]*.1;
  tt=[0:N]*.1;
  Pu1=[x_best(1,1),x_best(1,2),x_best(1,3);x_best(1,2),x_best(1,4),x_best(1,5);x_best(1,3),x_best(1,5),x_best(1,6)];
  Pd2=[x_best(1,7),x_best(1,8),x_best(1,9);x_best(1,8),x_best(1,10),x_best(1,11);x_best(1,9),x_best(1,11),x_best(1,12)];
  [states]=FUNCION_DA3(Pu1, Pd2, N,0);
  x5in=states(1,:);
  D=states(2,:);
  x11=states(3,:);
  x21=states(4,:);
  x31=states(5,:);
  Y1N=states(6,:);
  X2con=states(7,:);
  S2con=states(8,:);
  ICcon=states(9,:);
  QCH4con=states(10,:);
  figure(6);
  subplot(2,1,1),
  plot(t,x5in,'g','LineWidth',2),xlabel('Time (h)','fontsize',20),ylabel('b_i_n_c (mol/L)','fontsize',20),grid on;
  set(gca,'fontsize',18)
  v=[0 1500 0.04 .1 ]; axis(v);
  subplot(2,1,2),
  plot(t,D,'g','LineWidth',2),xlabel('Time (h)','fontsize',20),ylabel('D_i_n (1/h)','fontsize',20),grid on;
  set(gca,'fontsize',18)
  v=[0 1500 0 .2 ]; axis(v);
  print -dtiff -r300 binc_Dinc
  rf1=X2con;
  rf2=S2con;
  rf3=ICcon;
  rf4=QCH4con;   
  figure(1)
  plot(t,x11,'--g',t,rf1,'g','LineWidth',2),xlabel('Time (h)','fontsize',20),ylabel('X_2 (mol/L)','fontsize',20),grid on;
  legend('Neural estimated','Reference trajectory', 1);
  set(gca,'fontsize',18)
  v=[0 1500 0 .02 ]; axis(v);
  print -dtiff -r300 X2c
  figure(2)
  plot(t,x21,'--g',t,rf2,'g','LineWidth',2),xlabel('Time (h)','fontsize',20),ylabel('S_2 (mol/L)','fontsize',20),grid on;
  legend('Neural estimated','Reference trajectory', 1);
  set(gca,'fontsize',18)
  v=[0 1500 0 .02 ]; axis(v);
  print -dtiff -r300 S2c
  figure(3)
  plot(t,x31,'--g',t,rf3,'g','LineWidth',2),xlabel('Time (h)','fontsize',20),ylabel('IC (mol/L)','fontsize',20),grid on;
  legend('Neural estimated','Reference trajectory', 1);
  set(gca,'fontsize',18)
  v=[0 1500 .06 .11 ]; axis(v);
  print -dtiff -r300 ICc
  figure(4)
  plot(t,Y1N,'--g',t,rf4,'g','LineWidth',2),xlabel('Time (h)','fontsize',20),ylabel('YCH_4 (mol/h)','fontsize',20),grid on;
  legend('Neural estimated','Reference trajectory', 1);
  set(gca,'fontsize',18)
  v=[0 1500 0 .05 ]; axis(v);
  print -dtiff -r300 metanoc

  %disp(sprintf('Corrida=%d, Fitnes=%e, Evaluaciones:%d\n',ind5,Fitness_best,EvaluacioneFO));
  %semilogy(CONVERGENCIA,':k','LineWidth',2)%,title('DE'),hold on,figure,bar(x_best,'k','LineWidth',2)
  corridas(ind5)=Fitness_best; 
  end
  %plot(CONVERGENCIA(1:k-1),'k','LineWidth',2),title('DE'),hold on,figure,bar(x_best,'k','LineWidth',2)
  %disp(sprintf('DE:Func:%d,EvalFO:%d,Best:%.3e\n',funcion,EvalFO,Fitness_best));
  %end
end
