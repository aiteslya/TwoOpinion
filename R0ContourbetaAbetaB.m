%this code generates  R0 for permutation of the following parameters:
% betaA and betaB
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
num=40;
pAarr=[0.3,0.4,0.55];
pB=0.4;
k=1.6;
thetaA=5;
thetaB=5;
betaA=0.8;
betaB=2;
gammaA=1;
gammaB=1;

omega=0;
C=40;

%set up of the parameter array

arraybetaA=[0.01,0.8;0.2,0.6];
arraybetaB=[1,6;5,7];
arraypA=[0.4;0.283];
betaAarr=linspace(0.2,0.6,num);

betaBarr=linspace(5,7,num);

[BetaAarr,BetaBarr]=meshgrid(betaAarr,betaBarr);

for c1=1:2
    betaAarr=linspace(arraybetaA(c1,1),arraybetaA(c1,2),num);

    betaBarr=linspace(arraybetaB(c1,1),arraybetaB(c1,2),num);
    pA=arraypA(c1,1);
    [BetaAarr,BetaBarr]=meshgrid(betaAarr,betaBarr);
    Res=zeros(num,num);
    SA=zeros(num,num);
    for i1=1:1:num
        for i2=1:1:num
           betaA=BetaAarr(i1,i2);
           betaB=BetaBarr(i1,i2);
           %need to calculate sa and sb
           parsNa=[k,thetaA,thetaB,pA,pB];
           sa=NA(parsNa);
           sb=1-sa;
           SA(i1,i2)=sa;
           parsR0=[betaA,betaB,omega,pA,pB,k,thetaA,thetaB,gammaA,gammaB,C];
           r0=R0(sa,parsR0);
           Res(i1,i2)=r0;
           %add here investigation for R0a: it will be marked with dark blue
           %line on the screen
        end
    end

    figure(1);
    subplot(1,2,c1);
    if c1==1
        [M,c]=contourf(BetaAarr,BetaBarr,Res,sort([0,1.4,linspace(0.95*min(min(Res)),1.05*max(max(Res)),5)]),'ShowText','on');
    else
        [M,c]=contourf(BetaAarr,BetaBarr,Res,sort([0,1.4,linspace(0.95*min(min(Res)),1.05*max(max(Res)),7)]),'ShowText','on');
    end
    hold on;

    c.LineWidth=1;
    c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
    clabel(M,c,'FontSize',25,'FontWeight','bold');
    %highlight R0=1
    v = [1,1];
    [M,c]=contour(BetaAarr,BetaBarr,Res,v,'r','ShowText','on');
    c.LineWidth=6;
    clabel(M,c,'FontSize',25,'FontWeight','bold');
    
    colormap(parula);
    caxis([0.1,6]);
    if c1==2
        colorbar
    end
    %colorbar;
    xlabel('Infection rate of $$N_{a}$$ individuals, $$\beta_{a}$$','interpreter','latex','FontWeight','bold');
    ylabel('Infection rate of $$N_{b}$$ individuals, $$\beta_{b}$$','interpreter','latex','FontWeight','bold');
    
    set(gca,'FontSize',25);

end