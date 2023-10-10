%**************************************************************************
%         MODELO REACTOR ANAEROBICO DISCRETO CON CONTROL OPTIMO INVERSO
%**************************************************************************
function [fitness]=FUNCION_DA3(Pu1,Pd2, N,fitne)
    %load referencias_RHONO;
    load X2pt5
    load S2pt5
    load ICpt5
    load QCH4pt5
    load x2_1
    load s2_1
    load ic_1
    load qch4_1
    load X21P5
    load S21P5
    load IC1P5
    load metano
    load x3r
    load x4r_modificado
    load x5r_modificado
    load y1r
    load X2k27
    load S2r
    load ICr
    load QCH4_27
    T=0.092;
    k=1;                %Iteration, k=1 : 15000
    %%%%%%%%%% PARAMETER FOR AD MODEL (Kinetic parameters):%%%%%%%%%%%%%%%%
    mu1max = 0.1845;    %Afinidad de la biomasa por el sustrato
    ks1 = 0.26;         %Constante de saturación de crecimiento
    kin1 = 16.3333e-4;  %Constante de inhibición por exceso de sustrato
    mu2max = .0188;     %Afinidad de la biomasa por el sustrato
    ks2 = 2.18e-5;      %Constante de saturación de crecimiento
    kin2 = 8.22e-3;     %Constante de inhibición por exceso de sustrato    
    %%%%%%%%%% PARAMETER FOR AD MODEL (Equation parameters):%%%%%%%%%%%%%%%
    kd1 = 0.035;        %Tasa de mortalidad de la biomasa X1 
    kd2 = 0.0085;       %Tasa de mortalidad de la biomasa X2
    R1=.99;             %.99  1.05  .38          
    R2=.99;             %.99  .96  .36
    R3=345;             %350 SENSIBLE PARA OBSERVADOR 345
    R4=0.0666;          %.0663 SENSIBLE PARA OBSERVADOR .0666
    R5=0.0005;          %.0005
    R6=5;               %5       
    ka=1.7e-5; %Constante de equilibrio ácido-base del medio entre S- y HS
    kb=1.7e-7; %Constante de equilibrio ácido-base del medio entre B y CO2d
    kh=0.065;           %Constante de Henry
    Pt=1;               %Presión atmosférica    
    %%%%%%%%%%% CONDICIONES INICIALES DE LOS ESTADOS BIOLÓGICOS %%%%%%%%%%%
    pH_0 = 7;     %% pH inicial
    D_0 = 0.1;    %% Valor inicial Tasa de dilución
    x2in_0 = 10;  %% Valor inicial de entrada de sustrato Sc
    x4in_0 = 0.07;%% Valor inicial de entrada de sustrato S2
    Hplus_0 = 10^(-pH_0);    %% Ión Hidrógeno inicial
    HS_0=min(roots([kd2/kin2 (kd2-mu2max) kd2*ks2]));%FormaAcida(noIonizada) del sustrato S2
    x40 = HS_0*(ka+Hplus_0)/Hplus_0;    %S2 Sustrato rápidamente degradable
    x20=kd1*ks1/(mu1max-kd1-((kd1/kin1)*HS_0));%Sc SustratoLentamenteDegrad
    x10 = D_0*(x2in_0-x20)/(R6*kd1);    %Xc Microorganismos acidogénicos
    x30=(R4*kd1*x10+D_0*(x4in_0-x40))/(R3*kd2);%X2 MIcroorganismoMetanogéni
    %De las ecuaciones (1) y (2)de equilibrio químico entre ácidos y bases:
    Smin_0 = x40*ka/(ka+Hplus_0);%S- Forma base (ionizada) del sustrato S2
    % A partir de la ecuación de segundo grado para CO2D:      
    CO2D_0 = min(roots([1 -(kh*Pt*D_0+R1*R3*kd2*x30+R2*R3*kd2*x30+...
        R5*kd1*x10+Smin_0*D_0)/D_0 (R2*R3*kd2*x30+R5*kd1*x10+...
        Smin_0*D_0)*kh*Pt/D_0]));% CO2D  Dióxido de carbono disuelto
    % De la ecuación (4) de equilibrio entre bicarbonato y CO2:
    Bic_0 = CO2D_0*kb/Hplus_0;      % B  Bicarbonato  
    % De la ecuación (5) de cationes presentes en el reactor:
    x60 = Bic_0+Smin_0;
    % De la ecuación (3) de equilibrio de generación de carbono:
    x50 = Bic_0+CO2D_0;  % IC Carbono Inorgánico
    x1(k)=x10; x2(k)=x20; x3(k)=x30; x4(k)=x40; x5(k)=x50; x6(k)=x60;
    X0=[x1(k) x2(k) x3(k) x4(k) x5(k) x6(k)];    
    %%%%%%%% PARAMETER FOR NEURAL OBSERVER:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w1(:,k)=[x30;x30;x30;x30;x30];      %[w11, w12, w13, w14, w15]
    w2(:,k)=[x40;x40;x40;x40;x40;x40];  %[w21, w22, w23, w24, w25, w26]
    w3(:,k)=[x50;x50;x50;x50;x50];      %[w31, w32, w33, w34, w35]
    p1=1500;p2=1000;p3=1500;q1=1.5;q2=1.5;q3=0.2;%Matrices d Covarianza PyQ
    P1(:,:,k)=p1*eye(5);     Q1=q1*eye(5);   R11=1.5*[10  0;0 10];
    P2(:,:,k)=p2*eye(6);     Q2=q2*eye(6);   R22=1.5*[10 0;0 10];
    P3(:,:,k)=p3*eye(5);     Q3=q3*eye(5);   R33=1.5*[1  0;0  1];
    gm1=[0.12 0.12]; gm2=[0.09 0.09]; gm3=[0.09 0.09];
    eta1=2; eta2=1; eta3=10;
    X2_(k)=x30*1.5;S2_(k)=x40*1.5;IC_(k)=x50*1.5;%CondicionInicialObservado      
    %%%%%%%%%%%% PARAMETER FOR CONTROLLERS:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gs1=1600; gs2=500;
    Rc1=0.6; Rc2=0.4;
    pu(1:2000)=3; pd(1:2000)=3;        
    eval Regulacion %Donde esta? OSUNA   
    B(1)=Bic_0;
for k=1:N
        a(k)=x6(k);
        b(k)=(ka+kb)*x6(k)-ka*x4(k)-kb*x5(k);
        c(k)=-ka*kb*(x4(k)+x5(k)-x6(k));
        AB=[a(k) b(k) c(k)];
        Hplus(k)=max(roots(AB));
        pH1=-log10(Hplus);
        B(k+1)=kb*x5(k)/(kb+Hplus(k));
        HS(k)=(Hplus(k)*x4(k))/(Hplus(k)+ka);
        mu1(k)=(mu1max*x2(k))/(ks1+x2(k)+(x2(k)*HS(k)/kin1));
        mu2(k)=(mu2max*HS(k))/(ks2+HS(k)+(HS(k)^2/kin2));
        CO2d(k)=Hplus(k)*x5(k)/(kb+Hplus(k));
        PCO2(k)=CO2d(k)/kh;
        lbd(k)=PCO2(k)/(Pt-PCO2(k));
        HS2(k)=(Hplus(k)*S2_(k))/(Hplus(k)+ka);
        mu22(k)=(mu2max*HS2(k))/(ks2+HS2(k)+(HS2(k)^2/kin2));
        lbd2(k)=lbd(k);
        HS3(k)=HS(k);
        mu23(k)=mu2(k);
        lbd3(k)=Hplus(k)*IC_(k)/(Pt*kh*kb+Pt*kh*Hplus(k)-Hplus(k)*IC_(k));
        x2in(k)=x2in_0;
        x6in(k)=x60;
        Y1(k)=R1*R3*mu2(k)*x3(k);
        Y2(k)=lbd(k)*R2*R3*mu2(k)*x3(k);
        YCH4_(k)=1*R1*R3*mu2(k)*X2_(k);
        YCO2_(k)=lbd(k)*R2*R3*mu2(k)*X2_(k);
        e1=[Y1(k)-YCH4_(k) Y2(k)-YCO2_(k)]';
        Y1N2(k)=R1*R3*mu22(k)*X2_(k);
        Y2N2(k)=lbd2(k)*R2*R3*mu22(k)*X2_(k);
        e2=[Y1(k)-Y1N2(k) Y2(k)-Y2N2(k)]';
        Y1N3(k)=R1*R3*mu23(k)*X2_(k);
        Y2N3(k)=lbd3(k)*R2*R3*mu23(k)*X2_(k);
        e3=[Y1(k)-Y1N3(k) Y2(k)-Y2N3(k)]';
        %derivada con respecto a X2_
        dY11N(k)=R1*R3*mu2(k);
        dY12N(k)=lbd(k)*R2*R3*mu2(k);
        %derivada con respecto a S2_                        
        dY21N(k)=(R1*R3*mu2max*kin2)*Hplus(k)*(ka+Hplus(k))*X2_(k)*...
            (2*ka*Hplus(k)*kin2*ks2-S2_(k)^2*Hplus(k)^2+ka^2*kin2*ks2+...
            Hplus(k)^2*kin2*ks2)/(S2_(k)*ka*Hplus(k)*kin2+...
            S2_(k)^2*Hplus(k)^2+S2_(k)*Hplus(k)^2*kin2+...
            2*ka*Hplus(k)*kin2*ks2+ka^2*kin2*ks2+Hplus(k)^2*kin2*ks2)^2;
        dY22N(k)=(R2*R3*mu2max*kin2)*Hplus(k)^2*IC_(k)*(ka+Hplus(k))*...
            S2_(k)*(1/(kb*kh*Pt-Hplus(k)*IC_(k)+kh*Hplus(k)*Pt))*...
            (2*ka*Hplus(k)*kin2*ks2-S2_(k)^2*Hplus(k)^2+ka^2*kin2*ks2+...
            Hplus(k)^2*kin2*ks2)/(S2_(k)*ka*Hplus(k)*kin2+...
            S2_(k)^2*Hplus(k)^2+S2_(k)*Hplus(k)^2*kin2+...
            2*ka*Hplus(k)*kin2*ks2+ka^2*kin2*ks2+Hplus(k)^2*kin2*ks2)^2;
        dY31N(k)=0;
        dY32N(k)=(R2*R3*mu2max*kin2*kh*Pt)*Hplus(k)^2*X2_(k)*...
            (ka+Hplus(k))*(kb+Hplus(k))/((Hplus(k)*IC_(k)-kb*kh*Pt-kh*...
            Hplus(k)*Pt)^2*(S2_(k)*ka*Hplus(k)*kin2+S2_(k)^2*...
            Hplus(k)^2*kin2+2*ka*Hplus(k)*kin2*ks2+ka^2*kin2*ks2+...
            Hplus(k)^2*kin2*ks2));
        z11=tanh(X2_(k));z12=tanh(X2_(k))^2;z13=tanh(IC_(k));z14=z12;
        z15=z12; z1=[z11;z12;z13;z14;z15];
        dz11=(1/cosh(X2_(k)))^2;dz12=2*tanh(X2_(k))*dz11;dz13=0;dz14=dz12;
        dz15=dz12; dz1=[dz11;dz12;dz13;dz14;dz15];
        z21=tanh(S2_(k));z22=tanh(S2_(k))^2;z23=tanh(IC_(k));z24=z22;
        z25=z22; z26=z22; z2=[z21;z22;z23;z24;z25;z26];
        dz21=(1/cosh(S2_(k)))^2; dz22=2*tanh(S2_(k))*dz21;dz23=0;dz24=dz22; 
        dz25=dz22; dz26=dz22; dz2=[dz21;dz22;dz23;dz24;dz25;dz26];
        z31=tanh(IC_(k));z32=tanh(IC_(k))^2; z33=tanh(X2_(k));z34=z32;
        z35=z32; z3=[z31;z32;z33;z34;z35];
        dz31=(1/cosh(IC_(k)))^2;dz32=2*tanh(IC_(k))*dz31;dz33=0;dz34=dz32; 
        dz35=dz32; dz3=[dz31;dz32;dz33;dz34;dz35];
        if k<2000 % ESTABILIZACIÓN.
            D(k)=D_0;
            binc(k)=0;
            D0(k)=D_0;   
            A2=0;
            S2pert1=x4in_0;
            x4in(k)=x4in_0;
            ms=memS2(x4in(k),X2pt5(k),S2pt5(k),ICpt5(k),QCH4pt5(k),...
                X2_1(k),S2_1(k),IC_1(k),QCH4_1(k),X2(k),S2(k),IC(k),...
                QCH41(k),X2k2(k),x4r(k),x5r(k),QCH4k2(k),X2k(k),S2r(k),...
                ICr(k),QCH427(k)); 
            x2ref2(k)=ms(:,1);
            s2ref2(k)=ms(:,2);
            ICref2(k)=ms(:,3);
            QCH4ref2(k)=ms(:,4);
        else     % DESPUÉS DE LA ESTABILIZACIÓN.      
            pH=pH1; %pH_0;
            De=D_0;
            S2ini=x4in_0;
            S2eq=x4in_0;
            X22(k)=x3(k);
            QCH4n1=Y1(k);
            t_sim=k;
            times1=k;
            t_deb=2000;
            T_bgn=2000;
            A2=1.9;
            S2pert1=x4in_0+x4in_0*A2;
            x4in(k)=S2pert1;
            ms=memS2(x4in(k),X2pt5(k),S2pt5(k),ICpt5(k),QCH4pt5(k),...
                X2_1(k),S2_1(k),IC_1(k),QCH4_1(k),X2(k),S2(k),IC(k),...
                QCH41(k),X2k2(k),x4r(k),x5r(k),QCH4k2(k),X2k(k),S2r(k),...
                ICr(k),QCH427(k)); 
            x2ref(k-1999)=ms(:,1);
            X2con=cat(2,x2ref2, x2ref);
            s2ref(k-1999)=ms(:,2);
            S2con=cat(2,s2ref2, s2ref);
            ICref(k-1999)=ms(:,3);
            ICcon=cat(2,ICref2, ICref);
            QCH4ref(k-1999)=ms(:,4);
            QCH4con=cat(2,QCH4ref2, QCH4ref);
            kk1=ms(5);
            kk=kk1;  % de 1.8 a 1.9 se requiere que kk<=2.6  para s2in=2.7 
          % el valor es kk=2.3, el contrl D entra en 199.9 y B entra en 201
          % la salida de COJX2 y D son proporcionales y muy grandes para
          % S2in=2.7
            ref1=1*[X2con(k);S2con(k);ICcon(k)];
            ref2=ref1;
            b11(:,k)=w1(4,k)'*z14;
            b12(:,k)=w1(5,k)'*z15;
            b21(:,k)=w2(4,k)'*z24;
            b22(:,k)=w2(6,k)'*z26;
            b31(:,k)=w3(4,k)'*z34;
            b32(:,k)=w3(5,k)'*z35;
            f11(:,k)=w1(1:3,k)'*z1(1:3)+gm1*e1;% 
            f12(:,k)=w2(1:3,k)'*z2(1:3)+w2(5,k)'*z22*x4in(k)+gm2*e2;%
            f13(:,k)=w3(1:3,k)'*z3(1:3)+gm3*e3;%
            f1=[f11(:,k);f12(:,k);f13(:,k)];
            f2=[f11(:,k);f12(:,k);f13(:,k)];
            g11=[b11(:,k);b21(:,k);b31(:,k)];
            g21=[b12(:,k);b22(:,k);b32(:,k)];
            fu=f1-ref1;
            fd=f2-ref2;
            R=0.4;
            if A2<=1   
                D(k)=D_0;
                binc(k)=0;
            else
                pu(:,k+1)=pu(:,k)+13*gs1*fu'*Pu1*g11*Rc1^2*g11'*fu/...
                    (2*Rc1+.01*pu(:,k)*g11'*Pu1*g11)^3;                                                       
                pmu=pu(:,k+1)*Pu1;
                u1=-kk*inv(Rc1+2.1*g11'*pmu*g11)*g11'*pmu*fu;
                pd(:,k+1)=pd(:,k)+.1*gs2*fd'*Pd2*g21*Rc2^2*g21'*fd/...
                    (.2*Rc2+.1*pd(:,k)*g21'*Pd2*g21)^3;
                pmd=pd(:,k+1)*Pd2;
                u2=-.2*inv(Rc2+.2*g21'*pmd*g21)*g21'*pmd*fd;
                Dinp(k-1999)=1*u1; 
                Dinp2=cat(2,D0, Dinp);
                D(k)=(1*Dinp2(k));
                Bic(k)=(1*u2);  
                binc(k)=(1*(Bic(k)));
                %FUNCIONES TAKAGI--SUGENO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                deltQ=deltQCH4(times1,QCH4n1,S2eq,S2pert1,R1,R2,R3,t_deb);
                CJX2(k)=COJX2(D(k),S2ini,A2,X22(k),t_sim,T_bgn);
                [vy]=TS_6reglasmod(binc(k),D(k),CJX2(k),deltQ,pH_0,De);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Df(k-1999)=1*vy(:,2);
                Df2=cat(2,D0, Df);
                binc(k)=1*vy(:,1);
                if Df2(k)>.1
                    Df2(k)=.1;
                end
                D(k)=Df2(k);
            end
        end
        Dcl(k)=D(k);
        Bcl(k)=x60+binc(k);
        Din(k)=Dcl(k);
        bin(k)=Bcl(k);
        dx11(:,k+1)=w1(1:3,k)'*dz1(1:3)+w1(4,k)'*dz14*Din(k)+w1(5,k)'*...
            dz15*bin(k)+z1.*[1;1;1;Din(k);bin(k)]; 
        dx21(:,k+1)=w2(1:3,k)'*dz2(1:3)+w2(4,k)'*dz24*Din(k)+w2(5,k)'*...
            dz25*x4in(k)+w2(6,k)'*dz26*bin(k)+z2.*...
            [1;1;1;Din(k);x4in(k);bin(k)];
        dx31(:,k+1)=w3(1:3,k)'*dz3(1:3)+w3(4,k)'*dz34*Din(k)+w3(5,k)'*...
            dz35*bin(k)+z3.*[1;1;1;Din(k);bin(k)];
        H1=[dx11(:,k)*dY11N(k)  dx11(:,k)*dY12N(k)];
        H1=real(H1);
        H2=[dx21(:,k)*dY21N(k)  dx21(:,k)*dY22N(k)];
        H2=real(H2);
        H3=[dx31(:,k)*dY31N(k)  dx31(:,k)*dY32N(k)];
        H3=real(H3);
        %EKF1
        M1=inv(R11+H1'*P1(:,:,k)*H1);
        K1=P1(:,:,k)*H1*M1;
        P1(:,:,k+1)=P1(:,:,k)-K1*H1'*P1(:,:,k)+Q1;
        w1(:,k+1)=w1(:,k)+eta1*K1*e1;
        %EKF2
        M2=inv(R22+H2'*P2(:,:,k)*H2);
        K2=P2(:,:,k)*H2*M2;
        P2(:,:,k+1)=P2(:,:,k)-K2*H2'*P2(:,:,k)+Q2;
        w2(:,k+1)=w2(:,k)+eta2*K2*e2;
        %EKF3
        M3=inv(R33+H3'*P3(:,:,k)*H3);
        K3=P3(:,:,k)*H3*M3;
        P3(:,:,k+1)=P3(:,:,k)-K3*H3'*P3(:,:,k)+Q3;
        w3(:,k+1)=w3(:,k)+eta3*K3*e3;
        x1(k+1)=x1(k)+((mu1(k)-kd1)*x1(k))*T;
        x2(k+1)=x2(k)+(-R6*mu1(k)*x1(k)+Din(k)*(x2in(k)-x2(k)))*T;
        x3(k+1)=x3(k)+((mu2(k)-kd2)*x3(k))*T;
        x4(k+1)=x4(k)+(-R3*mu2(k)*x3(k)+R4*mu1(k)*x1(k)+...
            Din(k)*(x4in(k)-x4(k)))*T;
        x5(k+1)=x5(k)+(R2*R3*mu2(k)*x3(k)+R5*mu1(k)*x1(k)-...
            lbd(k)*R1*R3*mu2(k)*x3(k)+Din(k)*(bin(k)-x5(k)))*T;
        x6(k+1)=x6(k);
% OBSERVADOR: 
X2_(k+1)=w1(1:3,k)'*z1(1:3)+w1(4,k)'*z14*Din(k)+w1(5,k)'*z15*bin(k)+gm1*e1;
S2_(k+1)=w2(1:3,k)'*z2(1:3)+w2(4,k)'*z24*Din(k)+w2(5,k)'*z25*x4in(k)+...
                                                w2(6,k)'*z26*bin(k)+gm2*e2;
IC_(k+1)=w3(1:3,k)'*z3(1:3)+w3(4,k)'*z34*Din(k)+w3(5,k)'*z35*bin(k)+gm3*e3;
end
if fitne
 fitness=sum(mahal(X2_(1:N)',X2con(1:N)')); %Mahalanobis distance
else
fitness=[bin; D; X2_(1:N);S2_(1:N);IC_(1:N); YCH4_;X2con; S2con; ICcon; QCH4con];
end