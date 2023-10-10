% Clonal Selection Algorithm, CUTONALA, José Valentín Osuna-Enciso, Julio,
% 2014.
function Fitness_best = CSA_DA 
    Np=20;                     % Population size   
    [r,d]=INICIA_DA;          % Limits of search space,Problem's dimension 
    bits=22;                    % Bits per individual per dimension
    L=bits*d;                   % Bits per individual
    dec=2^(bits);               % Decimal per individual         
    %B=cadeia(Np,bits*d,0,0,0);  % Initial population, binary format       
    maxIter = 5;               % Maximum iteration number (Generations)
    pm=0.1; per=0.0; fat=0.4;  % Parameters: MutationProbability,??,Clones
    pma=pm; itpm=maxIter; pmr=0.8;% Hypermutation parameters                
    k=0;                        % Iterations
    N=15000;                    % Number of steps in the reactor
    for ind1=1:Np
            defi1=0; defi2=0;
            while ~(defi1 && defi2)
                B(ind1,:)=cadeia(1,bits*d,0,0,0);  % Binary
                X(ind1,:) = decodify(B(ind1,:),bits,d,dec,r)';% Decimal
                Pu1=[X(ind1,1),X(ind1,2),X(ind1,3);X(ind1,2),X(ind1,4),...
                    X(ind1,5);X(ind1,3),X(ind1,5),X(ind1,6)];
                Pd2=[X(ind1,7),X(ind1,8),X(ind1,9);X(ind1,8),X(ind1,10),...
                    X(ind1,11);X(ind1,9),X(ind1,11),X(ind1,12)];
                defi1=defPos(Pu1); defi2=defPos(Pd2);
            end
    end
    % Generations
    while k <=  maxIter
        for ind1=1:Np            
            X(ind1,:) = decodify(B(ind1,:),bits,d,dec,r)';% Decimal
            Pu1=[X(ind1,1),X(ind1,2),X(ind1,3);X(ind1,2),X(ind1,4),...
                    X(ind1,5);X(ind1,3),X(ind1,5),X(ind1,6)];
            Pd2=[X(ind1,7),X(ind1,8),X(ind1,9);X(ind1,8),X(ind1,10),...
                    X(ind1,11);X(ind1,9),X(ind1,11),X(ind1,12)];
            fit(ind1) = FUNCION_DA3(Pu1,Pd2,N,1);% Function being optimized 
        end
        T=[];                         % Immune response, maturation
        [fx,ind] = sort(fit);   % Np best individuals (minimization)
        valx = X(ind(1:end),:); 
        % Reproduction
        [T,pcs] = reprod(Np,fat,Np,ind,B,T);
        % Hypermutation
        M = rand(size(T,1),L) <= pm;
        T = T - 2 .* (T.*M) + M;
        T(pcs,:) = B(fliplr(ind(end-Np+1:end)),:);
        % New Re-Selection (Multi-peak solution)
        Tm = decodify(T,bits,d,dec,r)';         % Mutated population
        for in2=1:size(Tm,1)
            Pu1=[Tm(in2,1),Tm(in2,2),Tm(in2,3);Tm(in2,2),Tm(in2,4),Tm(in2,5);Tm(in2,3),Tm(in2,5),Tm(in2,6)];
            Pd2=[Tm(in2,7),Tm(in2,8),Tm(in2,9);Tm(in2,8),Tm(in2,10),Tm(in2,11);Tm(in2,9),Tm(in2,11),Tm(in2,12)];
            defi1=defPos(Pu1); defi2=defPos(Pd2);
            if (defi1 && defi2) % Si son definidas positivas, evalúa      
                fit2(in2)=FUNCION_DA3(Pu1,Pd2,N,1);
            else
                fit2(in2)=fx(end);  
            end
        end
        pcs = [0 pcs];
        for i=1:Np,
            [out(i),bcs(i)] = min(fit2(pcs(i)+1:pcs(i+1)));  % Minimization
            bcs(i) = bcs(i) + pcs(i);
        end;
        B(fliplr(ind(end-Np+1:end)),:) = T(bcs,:);
        % Editing (Repertoire shift)
        nedit = round(per*Np); k = k + 1;
        B(ind(1:nedit),:) = cadeia(nedit,L,0,0,0);
        pm = pmcont(pm,pma,pmr,k,itpm); valfx = min(fx);
        %vpm = [vpm pm]; %vfx = [vfx valfx]; %vmfit = [vmfit mean(fit)];
        disp(sprintf('It.: %d, pm: %.4f, f(x,y): %2.6e',k,pm,valfx));
    end; % end while
    x_best = valx(1,:); Fitness_best = min(fx);
    t=[1:N]*.1;    
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
    plot(t,x5in,'b','LineWidth',2),xlabel('Time (h)','fontsize',20),ylabel('b_i_n_c (mol/L)','fontsize',20),grid on;
    set(gca,'fontsize',18)
    v=[0 1500 0.04 .1 ]; axis(v);
    subplot(2,1,2),
    plot(t,D,'b','LineWidth',2),xlabel('Time (h)','fontsize',20),ylabel('D_i_n (1/h)','fontsize',20),grid on;
    set(gca,'fontsize',18)
    v=[0 1500 0 .2 ]; axis(v);
    print -dtiff -r300 binc_Dinc
    rf1=X2con;
    rf2=S2con;
    rf3=ICcon;
    rf4=QCH4con;   
    figure(1)
    plot(t,x11,'--b',t,rf1,'b','LineWidth',2),xlabel('Time (h)','fontsize',20),ylabel('X_2 (mol/L)','fontsize',20),grid on;
    legend('Neural estimated','Reference trajectory');
    set(gca,'fontsize',18)
    v=[0 1500 0 .02 ]; axis(v);
    print -dtiff -r300 X2c
    figure(2)
    plot(t,x21,'--b',t,rf2,'b','LineWidth',2),xlabel('Time (h)','fontsize',20),ylabel('S_2 (mol/L)','fontsize',20),grid on;
    legend('Neural estimated','Reference trajectory');
    set(gca,'fontsize',18)
    v=[0 1500 0 .02 ]; axis(v);
    print -dtiff -r300 S2c
    figure(3)
    plot(t,x31,'--b',t,rf3,'b','LineWidth',2),xlabel('Time (h)','fontsize',20),ylabel('IC (mol/L)','fontsize',20),grid on;
    legend('Neural estimated','Reference trajectory');
    set(gca,'fontsize',18)
    v=[0 1500 .06 .11 ]; axis(v);
    print -dtiff -r300 ICc
    figure(4)
    plot(t,Y1N,'--b',t,rf4,'b','LineWidth',2),xlabel('Time (h)','fontsize',20),ylabel('YCH_4 (mol/h)','fontsize',20),grid on;
    legend('Neural estimated','Reference trajectory');
    set(gca,'fontsize',18)
    v=[0 1500 0 .05 ]; axis(v);
    print -dtiff -r300 metanoc
end
