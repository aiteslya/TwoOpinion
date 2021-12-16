% produce Holling I bifurcation diagram
% prepare state space

close all;
clear variables;
clc;
format long;

pB=0.4;
thetaA=5;
thetaB=5;
mesh_size=200;
%find pA for which there is no stable co-existence
pA1=pB/10;
pA2=pB;
pA3=2*pB;
pAarr1=linspace(pA1,pA2,mesh_size/2);
pAarr2=linspace(pA2,pA3,mesh_size/2);
pAarr=linspace(pA1,pA3,mesh_size);

%plot
figure(1);

h1=fill([pAarr1/pB fliplr(pAarr1/pB)],[pAarr1*0 pAarr1*0+1],'c');hold on;
set(h1,'facecolor',[251, 206, 177]/255) 
set(h1,'EdgeColor','none')

h1=fill([pAarr2/pB fliplr(pAarr2/pB)],[pAarr2*0 pAarr2*0+1],'c');hold on;
set(h1,'facecolor',[173, 216, 230]/255);%[208, 240, 192]/255)
set(h1,'EdgeColor','none')

% h1=fill([na3 fliplr(na3)],[na3*0 na3*0+1],'c');hold on;
% set(h1,'facecolor',[173, 216, 230]/255)
% set(h1,'EdgeColor','none')

h1=plot(pAarr1/pB,pAarr1*0,'bo');
set(h1,'markerFacecolor',get(h1,'color'));hold on;
h2=plot(pAarr1/pB,pAarr1*0+1,'ro');
set(h2,'markerFacecolor',get(h2,'color'));
%overplot
h1=plot(pAarr2/pB,pAarr2*0,'ro');
set(h1,'markerFacecolor',get(h1,'color'));hold on;
h2=plot(pAarr2/pB,pAarr2*0+1,'bo');
set(h2,'markerFacecolor',get(h2,'color'));hold on;

xlabel('$$p_{a}/p_{b}$$','interpreter','latex');
ylabel('Density of $$N_{a}$$ population,$$n_{a}$$','interpreter','latex');
title(['$$k=1$$, $$\theta_{a}=\theta_{b}=0$$'],'interpreter','latex');
set(gca,'FontSize',25);
xlim([pA1/pB,pA3/pB]);
