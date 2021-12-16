%this code generates  R0 contour maps for permutation of the following parameters:
% c,omega,betaA: 2 figures in total. For Omega we take only two settings:
% omega=0.1 and omega=0.9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heterogeneous contact rate with assortativity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear variables;
clc;
format long;

%set up integration options
Atol=1e-12;
RelTol=1e-10;
opts = odeset('RelTol',RelTol,'AbsTol',Atol);
calc_err=5e-11;
fig_c=1;
%%
%c-omega
%set up the parameters
num=30;
pA=0.4;
pB=0.4;
k=1.6;
thetaA=5;
thetaB=5;
betaA=0.9;
betaB=1.2;
gammaA=1;
gammaB=1;

%check R0 for no assortativity
omega=0;
c=5;
parsR0=[betaA,betaB,omega,pA,pB,k,thetaA,thetaB,gammaA,gammaB,c];
r0=R0(0.5,parsR0);

%set up of the parameter array
carr=linspace(1,105,num);

%omega=0.1
%set up the parameters
omega=0.1;

%set up of the parameter array
betaAarr=linspace(0.01,betaB,num);

[BetaAarr,Carr]=meshgrid(betaAarr,carr);

%heterogeneous contacts with assortativity
Res=zeros(num,num);
SA=zeros(num,num);
for i1=1:1:num
    for i2=1:1:num
       betaA=BetaAarr(i1,i2);
       c=Carr(i1,i2);
       %need to calculate sa and sb
       parsNa=[k,thetaA,thetaB,pA,pB];
       sa=NA(parsNa);
       sb=1-sa;
       SA(i1,i2)=sa;
       parsR0=[betaA,betaB,omega,pA,pB,k,thetaA,thetaB,gammaA,gammaB,c];
       r0=R0(sa,parsR0);
       Res(i1,i2)=r0;
       %add here investigation for R0a: it will be marked with dark blue
       %line on the screen
    end
end

figure(fig_c);[M,c]=contourf(BetaAarr,Carr,Res,sort([0,[0.6,0.7,0.8 0.9 1.1]]),'ShowText','on');
hold on;

c.LineWidth=1;
c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
clabel(M,c,'FontSize',25,'FontWeight','bold');
%highlight R0=1
v = [1,1];
[M,c]=contour(BetaAarr,Carr,Res,v,'r','ShowText','on');
c.LineWidth=6;
clabel(M,c,'FontSize',25,'FontWeight','bold');
caxis([0.6,1.21]);
colormap(spring);
%colorbar;
xlabel('$$\beta_{a}$$','interpreter','latex','FontWeight','bold');
ylabel('$$c$$','interpreter','latex','FontWeight','bold');
title(['$$\omega=$$',num2str(omega)],'interpreter','latex','FontWeight','bold');
set(gca,'FontSize',25);

fig_c=fig_c+1;

%c-betaA: omega=0.9
%set up the parameters
omega=0.9;
%heterogeneous contacts with assortativity
Res=zeros(num,num);
SA=zeros(num,num);
for i1=1:1:num
    for i2=1:1:num
       betaA=BetaAarr(i1,i2);
       c=Carr(i1,i2);
       %need to calculate sa and sb
       parsNa=[k,thetaA,thetaB,pA,pB];
       sa=NA(parsNa);
       sb=1-sa;
       SA(i1,i2)=sa;
       parsR0=[betaA,betaB,omega,pA,pB,k,thetaA,thetaB,gammaA,gammaB,c];
       r0=R0(sa,parsR0);
       Res(i1,i2)=r0;
       %add here investigation for R0a: it will be marked with dark blue
       %line on the screen
    end
end

figure(fig_c);[M,c]=contourf(BetaAarr,Carr,Res,sort([0,[0.6,0.7,0.8 0.9 1.1]]),'ShowText','on');
hold on;

c.LineWidth=1;
c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
clabel(M,c,'FontSize',25,'FontWeight','bold');
%highlight R0=1
v = [1,1];
[M,c]=contour(BetaAarr,Carr,Res,v,'r','ShowText','on');
c.LineWidth=6;
clabel(M,c,'FontSize',20,'FontWeight','bold');
caxis([min(min(Res)),max(max(Res))]);
colormap(spring);
colorbar;
xlabel('$$\beta_{a}$$','interpreter','latex','FontWeight','bold');
ylabel('$$c$$','interpreter','latex','FontWeight','bold');
title(['$$\omega=$$',num2str(omega)],'interpreter','latex','FontWeight','bold');
set(gca,'FontSize',25);

fig_c=fig_c+1;

