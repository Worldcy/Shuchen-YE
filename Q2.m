clear;
clc;
close all;

h=1e-5;      %Step length
t=0:h:0.5;   %generate a vector of the independent variable t

%% Create an array of calculation results x, y, z
N=length(t);
E=ones(1,N);
S=10*ones(1,N);
ES=zeros(1,N);
P=zeros(1,N);

%% Fourth-order Runge-Kutta iteration
for i=2:N
    t_n=t(i-1);
    E_n=E(i-1);
    S_n=S(i-1);
    ES_n=ES(i-1);
    P_n=P(i-1);
    
    kE1_1=100*E_n*S_n;
    kE1_2=600*ES_n;
    kE1_3=150*ES_n;
    kS1_1=100*E_n*S_n;
    kS1_2=600*ES_n;
    kES1_1=100*E_n*S_n;
    kES1_2=600*ES_n;
    kES1_3=150*ES_n;
    kP1=150*ES_n;
    
    kE2_1=100*(E_n+kE1_1*h/2)*(S_n+kS1_1*h/2);
    kE2_2=600*(ES_n+kES1_2*h/2);
    kE2_3=150*(ES_n+kES1_3*h/2);
    kS2_1=100*(E_n+kE1_1*h/2)*(S_n+kS1_1*h/2);
    kS2_2=600*(ES_n+kES1_2*h/2);
    kES2_1=100*(E_n+kE1_1*h/2)*(S_n+kS1_1*h/2);
    kES2_2=600*(ES_n+kES1_2*h/2);
    kES2_3=150*(ES_n+kES1_3*h/2);
    kP2=150*(ES_n+kES1_2*h/2);
    
    kE3_1=100*(E_n+kE2_1*h/2)*(S_n+kS2_1*h/2);
    kE3_2=600*(ES_n+kES2_2*h/2);
    kE3_3=150*(ES_n+kES2_3*h/2);
    kS3_1=100*(E_n+kE2_1*h/2)*(S_n+kS2_1*h/2);
    kS3_2=600*(ES_n+kES2_2*h/2);
    kES3_1=100*(E_n+kE2_1*h/2)*(S_n+kS2_1*h/2);
    kES3_2=600*(ES_n+kES2_2*h/2);
    kES3_3=150*(ES_n+kES2_3*h/2);
    kP3=150*(ES_n+kES2_2*h/2);
    
    kE4_1=100*(E_n+kE3_1*h/2)*(S_n+kS3_1*h/2);
    kE4_2=600*(ES_n+kES3_2*h/2);
    kE4_3=150*(ES_n+kES3_3*h/2);
    kS4_1=100*(E_n+kE3_1*h/2)*(S_n+kS3_1*h/2);
    kS4_2=600*(ES_n+kES3_2*h/2);
    kES4_1=100*(E_n+kE3_1*h/2)*(S_n+kS3_1*h/2);
    kES4_2=600*(ES_n+kES3_2*h/2);
    kES4_3=150*(ES_n+kES3_3*h/2);
    kP4=150*(ES_n+kES3_2*h/2);
    
    E(i)=E_n-h/6*(kE1_1+2*kE2_1+2*kE3_1+kE4_1)+h/6*(kE1_2+2*kE2_2+2*kE3_2+kE4_2)+h/6*(kE1_3+2*kE2_3+2*kE3_3+kE4_3);
    S(i)=S_n-h/6*(kS1_1+2*kS2_1+2*kS3_1+kS4_1)+h/6*(kS1_2+2*kS2_2+2*kS3_2+kS4_2);
    ES(i)=ES_n+h/6*(kES1_1+2*kES2_1+2*kES3_1+kES4_1)-h/6*(kES1_2+2*kES2_2+2*kES3_2+kES4_2)-h/6*(kES1_3+2*kES2_3+2*kES3_3+kES4_3);
    P(i)=P_n+h/6*(kP1+2*kP2+2*kP3+kP4);
end
for i=2:N
    T(i-1)=P(i)-P(i-1);
end
a=find(T==max(T))-1;
%% Drawing
figure();
hold on;
plot(t,E,'r');
plot(t,S,'g');
plot(t,ES,'b');
plot(t,P,'m');
text(t(a),P(a),'*','color','r','FontSize',20);
text(t(a),P(a),['(',num2str(t(a)),',',num2str(P(a)),')']);
legend('E','S','ES','P');
xlabel('t');
title('Substances change over time');
hold off;

(P(a+1)-P(a))/(t(a+1)-t(a))
