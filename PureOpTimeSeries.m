%This script generates time series for the competition of two opinions in
%the disease free population. These are the final figures that will go into
%the manuscript
%%
%prepare the state space
clc;
clear variables;
close all;
format long;

%set up integration options
Atol=1e-12;
RelTol=1e-10;
opts = odeset('RelTol',RelTol,'AbsTol',Atol);

%set up the colors
num_traj=12;
orange=[249,87,0]/255;
%spanish blue;
blue=[0,164,204]/255;
colors_p = [linspace(orange(1),blue(1),num_traj)',linspace(orange(2),blue(2),num_traj)', linspace(orange(3),blue(3),num_traj)'];

%% Holling type I
% set parameters
pB=0.1;
pA=1.5*pB;
c=10;
thetaA=0;
thetaB=0;
k=1;
omega=0;

pars=[c,pA,pB,thetaA,thetaB,k,omega];
nainit=linspace(0,1,num_traj);
nbinit=1-nainit;
T=20;
for counter=1:num_traj
   init=[nainit(counter), nbinit(counter)];
   [t,y]=ode45(@(t,y)TwoOpAss(t,y,pars),[0,T], init,opts);
   figure(1);plot(t,y(:,1),'LineWidth',4,'color',colors_p(counter,:));hold on;
end
xlabel('Time','interpreter','latex');
ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','interpreter','latex');
title(['$$k=',num2str(k),'$$, $$\theta_{a}=\theta_{b}=',num2str(thetaA),'$$, $$p_{a}/p_{b}=',num2str(pA/pB),'$$, H1'],'interpreter','latex');
set(gca,'FontSize',25);

%% Holling type II
% two figures with co-existence and without co-existence equilibrium
% with co-existence equilibrium
pB=0.1;
thetaA=5;
thetaB=5;
pA=1.5*pB;%or pA=pB - for co-existence
pars=[c,pA,pB,thetaA,thetaB,k,omega];
T=20;
%stable co-existence equilibrium
for counter=1:num_traj
   init=[nainit(counter), nbinit(counter)];
   [t,y]=ode45(@(t,y)TwoOpAss(t,y,pars),[0,T], init,opts);
   figure(3);plot(t,y(:,1),'LineWidth',4,'color',colors_p(counter,:));hold on;
end
xlim([0,T]);
xlabel('Time','interpreter','latex');
ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','interpreter','latex');
title(['H2, $$p_{a}/p_{b}=',num2str(pA/pB),'$$'],'interpreter','latex');
set(gca,'FontSize',25);

%% Holling type H3 with k1=1.6
% two figures with one or three co-existence equilibria
% with three co-existence equilibria
pB=0.1;
thetaA=5;
thetaB=5;
k=1.6;
pA=pB;
pars=[c,pA,pB,thetaA,thetaB,k,omega];
T=60;
%stable co-existence equilibrium
for counter=1:num_traj
   init=[nainit(counter), nbinit(counter)];
   [t,y]=ode45(@(t,y)TwoOpAss(t,y,pars),[0,T], init,opts);
   figure(4);plot(t,y(:,1),'LineWidth',4,'color',colors_p(counter,:));hold on;
end
xlim([0,T]);
xlabel('Time','interpreter','latex');
ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','interpreter','latex');
title(['$$k=',num2str(k),'$$, $$\theta_{a}=\theta_{b}=',num2str(thetaA),'$$, $$p_{a}/p_{b}=',num2str(pA/pB),'$$, H3'],'interpreter','latex');
set(gca,'FontSize',25);


% with one st co-existence equilibria
pB=0.1;
thetaA=5;
thetaB=5;
k=1.6;
pA=1.5*pB;
pars=[c,pA,pB,thetaA,thetaB,k,omega];
T=120;
%stable co-existence equilibrium
for counter=1:num_traj
   init=[nainit(counter), nbinit(counter)];
   [t,y]=ode45(@(t,y)TwoOpAss(t,y,pars),[0,T], init,opts);
   figure(5);plot(t,y(:,1),'LineWidth',4,'color',colors_p(counter,:));hold on;
end
xlim([0,T]);
xlabel('Time','interpreter','latex');
ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','interpreter','latex');
title(['$$k=',num2str(k),'$$, $$\theta_{a}=\theta_{b}=',num2str(thetaA),'$$, $$p_{a}/p_{b}=',num2str(pA/pB),'$$, H3'],'interpreter','latex');
set(gca,'FontSize',25);

%% Holling type H3 with k1=2.7
% two figures both with one co-existence equilibria with different levels
% of n_{a}
% with three co-existence equilibria
pB=0.1;
thetaA=5;
thetaB=5;
k=2.7;
pA=pB;
pars=[c,pA,pB,thetaA,thetaB,k,omega];
T=60;
%stable co-existence equilibrium
for counter=1:num_traj
   init=[nainit(counter), nbinit(counter)];
   [t,y]=ode45(@(t,y)TwoOpAss(t,y,pars),[0,T], init,opts);
   figure(6);plot(t,y(:,1),'LineWidth',4,'color',colors_p(counter,:));hold on;
end
xlim([0,T]);
xlabel('Time','interpreter','latex');
ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','interpreter','latex');
title(['$$k=',num2str(k),'$$, $$p_{a}/p_{b}=',num2str(pA/pB),'$$, H3'],'interpreter','latex');
set(gca,'FontSize',25);


% with one st co-existence equilibria
pB=0.1;
thetaA=5;
thetaB=5;
k=2.7;
pA=1.5*pB;
pars=[c,pA,pB,thetaA,thetaB,k,omega];
T=120;
%stable co-existence equilibrium
for counter=1:num_traj
   init=[nainit(counter), nbinit(counter)];
   [t,y]=ode45(@(t,y)TwoOpAss(t,y,pars),[0,T], init,opts);
   figure(7);plot(t,y(:,1),'LineWidth',4,'color',colors_p(counter,:));hold on;
end
xlim([0,T]);
xlabel('Time','interpreter','latex');
ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','interpreter','latex');
title(['$$k=',num2str(k),'$$, $$p_{a}/p_{b}=',num2str(pA/pB),'$$, H3'],'interpreter','latex');
set(gca,'FontSize',25);

