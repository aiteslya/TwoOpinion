%this code generates  R0 for permutation of the following parameters:
% c,omega, pA0 2 figures in total. For pA0 we take only 3 settings settings:
% 0.3,0.4, 0.55
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

%set up of the parameter array
carr=linspace(1,105,num);

omegaArr=linspace(0.1,0.9,num);

[OmegaAarr,Carr]=meshgrid(omegaArr,carr);

for pA=1.25*pB%pAarr
    Res=zeros(num,num);
    SA=zeros(num,num);
    for i1=1:1:num
        for i2=1:1:num
           omega=OmegaAarr(i1,i2);
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

    figure(fig_c);[M,c]=contourf(OmegaAarr,Carr,Res,sort([0,1.4,linspace(0.95*min(min(Res)),1.05*max(max(Res)),6)]),'ShowText','on');
    hold on;

    c.LineWidth=1;
    c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
    clabel(M,c,'FontSize',25,'FontWeight','bold');
    %highlight R0=1
    v = [1,1];
    [M,c]=contour(OmegaAarr,Carr,Res,v,'r','ShowText','on');
    c.LineWidth=6;
    clabel(M,c,'FontSize',25,'FontWeight','bold');
    %caxis([0.6,1.21]);
    colormap(parula);
    %colorbar;
    xlabel('$$\omega$$','interpreter','latex','FontWeight','bold');
    ylabel('$$c$$','interpreter','latex','FontWeight','bold');
    title(['$$p_{a}(0)=$$',num2str(pA)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);

    fig_c=fig_c+1;
end