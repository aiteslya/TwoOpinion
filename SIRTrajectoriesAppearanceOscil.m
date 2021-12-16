%this code generates time series of infected and na individuals in different regions of
%betaB space
% for assortativity 0, c=80 and m=90 and pA0=0.283 and betaA=0.2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare workspace
close all;
clear variables;
clc;
format long;

%set up integration options
Atol=1e-12;
RelTol=1e-10;
opts = odeset('RelTol',RelTol,'AbsTol',Atol);
calc_err=1e-9;
fig_c=1;
%
%set up the parameters
pB=0.4;
pA0=0.4;
pA1=1;

k=1.6;
thetaA=5;
thetaB=5;
betaA=0.2;
betaBarr=[3.34237,3.6,4,4.15,4.3];
C=40;
m=75;
gammaA=1;
gammaB=1;
omega=0;

%need to calculate sa and sb to set as initial conditions
parsNa=[k,thetaA,thetaB,pA0,pB];
sa=NA(parsNa);
sb=1-sa;
init=[sa,0,0,sb-6e-8,6e-8,0];
T=50;

for counter=1:numel(betaBarr)
    betaB=betaBarr(counter);
    pars=[C,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
    [t,y]=ode45(@(t,y)TwoOpAssSIR(t,y,pars),[0,T], init,opts);
    infect=y(:,2)+y(:,5);
    na=y(:,1)+y(:,2)+y(:,3);
    
    figure(fig_c);
    yyaxis left
    plot(t,infect,'LineWidth',4);
    ylabel('Prevalence, $$i$$','interpreter','latex','FontWeight','bold');
    ylim([0,1e-2]);
    hold on;
    yyaxis right
    plot(t,na,'LineWidth',4);
    ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','interpreter','latex','FontWeight','bold');
    ylim([0,1]);
    xlabel('Time','interpreter','latex','FontWeight','bold');
    
    title(['Infection rate of $$N_{b}$$ individuals $$\beta_{b}=$$',num2str(betaB,2)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);

    fig_c=fig_c+1;
end

%second set of figure in betaA
betaB=4.15;
betaAarr=[0.01,0.1,0.15,0.2,0.25];

for counter=1:numel(betaAarr)
    betaA=betaAarr(counter);
    pars=[C,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
    [t,y]=ode45(@(t,y)TwoOpAssSIR(t,y,pars),[0,T], init,opts);
    infect=y(:,2)+y(:,5);
    na=y(:,1)+y(:,2)+y(:,3);
    
    figure(fig_c);
    yyaxis left
    plot(t,infect,'LineWidth',4);
    ylabel('Prevalence, $$i$$','interpreter','latex','FontWeight','bold');
    ylim([0,1e-2]);
    hold on;
    yyaxis right
    plot(t,na,'LineWidth',4);
    ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','interpreter','latex','FontWeight','bold');
    ylim([0,1]);
    xlabel('Time','interpreter','latex','FontWeight','bold');
    
    title({'Transmission rate to health-positive';['individuals $$\beta_{a}=$$',num2str(betaA,2)]},'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);

    fig_c=fig_c+1;
end
