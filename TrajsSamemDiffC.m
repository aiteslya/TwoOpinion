%This script generates time series for two opinion system coupled with SIR
%infectious disease. These time series are outputed as ease as a phase
%curves
clc;
clear variables;
close all;
format long;
%set up the parameters

thetaA=5;
thetaB=5;
betaA=0.8;
betaB=2;
gammaA=1;
gammaB=1;
k=1.6;
theta=5;
pB=0.4;
omega=0.9;

pA0=0.5;
pA1=1;
Carr=[5,10,20,40,80,100];

m=10;
%set up integration options
Atol=1e-12;
RelTol=1e-10;
opts = odeset('RelTol',RelTol,'AbsTol',Atol);
calc_err=5e-11;
%%

%set up the parameters

figc=1;

ccounter=1;
for C=Carr
       %calculate initial conditions - stable co-existence of
       %opinions
       parsNa=[k,thetaA,thetaB,pA0,pB];
       sa=NA(parsNa);
       sb=1-sa;
       init=[sa,0,0,sb-6e-8,6e-8,0];
       pars=[C,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
       T=500;
       [t,y]=ode45(@(t,y)TwoOpAssSIR(t,y,pars),[0,T], init,opts);
       infect=y(:,2)+y(:,5);
       na=y(:,1)+y(:,2)+y(:,3);
       pA=pA0+(pA0-pA1)*(m+1)*(1./(1+m*infect)-1)/m;
       figure(1);subplot(1,3,1);
       plot(t,na,'LineWidth',4);hold on;
       xlabel('Time','interpreter','latex')
       ylabel('$$n_{a}$$','interpreter','latex')
       set(gca,'FontSize',25);
       figure(1);subplot(1,3,2);
       plot(t,infect,'LineWidth',4);hold on;
       xlabel('Time','interpreter','latex')
       ylabel('$$i$$','interpreter','latex')
       set(gca,'FontSize',25);
       figure(1);subplot(1,3,3);
       plot(t,pA,'LineWidth',4);hold on;
       xlabel('Time','interpreter','latex')
       ylabel('$$p_{a}(t)$$','interpreter','latex')
       set(gca,'FontSize',25);
end
