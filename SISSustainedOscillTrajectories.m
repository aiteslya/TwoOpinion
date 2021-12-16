%this script produces time series trajetcories of existence of periodic solutions,
%and subsequent  for an SIS coupled with opinion
%competition searching for sustained oscillations in betA-betaB-m space
%prepare space
clc;
clear variables;
close all;
format long;

%set parameters
pB=0.4;
pA0=0.283;
pA1=0.6;
C=10;
m=75;
gammaA=1;
gammaB=1;
k=1.6;
thetaA=5;
thetaB=5;
betaAarr=[0.39,0.42,0.46,0.5];
betaB=5.5;
omega=0;
% set up the integration settings
Atol=1e-12;
RelTol=1e-10;
opts = odeset('RelTol',RelTol,'AbsTol',Atol);
%tolerance for calculation of the equilibrium
calc_err=1e-7;
T=2e2;
fig_c=1;
%need to calculate sa and sb to use as init

parsNa=[k,thetaA,thetaB,pA0,pB];
sa=NA(parsNa);
sb=1-sa;
init=[sa,0,sb-6e-8,6e-8];

for betaA=betaAarr
    pars=[C,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
    [t,y]=ode45(@(t,y)TwoOpAssSIS(t,y,pars),[0,T], init,opts);
    infect=y(:,2)+y(:,4);
    na=y(:,1)+y(:,2);
    figure(fig_c);
    yyaxis left
    plot(t,infect,'LineWidth',4);
    ylabel('Prevalence, $$i$$','interpreter','latex','FontWeight','bold');
    ylim([0,5e-1]);
    hold on;
    yyaxis right
    plot(t,na,'LineWidth',4);
    ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','interpreter','latex','FontWeight','bold');
    ylim([0,1]);
    xlabel('Time','interpreter','latex','FontWeight','bold');
    
    title(['Infection rate of $$N_{a}$$ individuals $$\beta_{a}=$$',num2str(betaA,2)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);

    fig_c=fig_c+1;
end
    