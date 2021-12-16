%this code generates contour map of peak n_{a} in different regions of
%betaA-betaB-m space
%for assortativity 0

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
%c-omega
%set up the parameters
pB=0.4;
pA0=0.4;%0.283;
pA1=1;

k=1.6;
thetaA=5;
thetaB=5;
meshsize=100;
betaAarr=linspace(0.01,0.8,meshsize);
betaBarr=linspace(1.1,6,meshsize);
C=40;
marr=[25,50,75];
gammaA=1;
gammaB=1;
omega=0;

[BetaAarr,BetaBarr]=meshgrid(betaAarr,betaBarr);
figc=1;
%define colors:

l=numel(marr);
darkblue = [18, 10, 143]/255;
lightblue =  [0.45,0.76,0.98];
col = [linspace(lightblue(1),darkblue(1),l)', linspace(lightblue(2),darkblue(2),l)', linspace(lightblue(3),darkblue(3),l)'];

count=1;
for m=marr

    NumelPeaks=zeros(meshsize,meshsize);
    SwitchA=zeros(meshsize,meshsize);
    for i1=1:1:meshsize
        for i2=1:1:meshsize
           betaA=BetaAarr(i1,i2);
           betaB=BetaBarr(i1,i2);
           %need to calculate sa and sb
           parsNa=[k,thetaA,thetaB,pA0,pB];
           sa=NA(parsNa);
           sb=1-sa;
           T=500;
           fl=0;
           init=[sa,0,0,sb-6e-8,6e-8,0];
           pars=[C,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
           while ~fl
            [t,y]=ode45(@(t,y)TwoOpAssSIR(t,y,pars),[0,T], init,opts);
            infect=y(:,2)+y(:,5);
            if infect(end)<6e-8
                fl=1;
            else
               T=2*T;
            end
           end
           na=y(:,1)+y(:,2)+y(:,3);
           [pks,locs] = findpeaks(infect,t);
           NumelPeaks(i1,i2)=numel(find(pks>6e-8));
           if abs(na(end)-1)<calc_err
               SwitchA(i1,i2)=1;
           end
        end
    end
    figure(fig_c);[M,c]=contourf(BetaAarr,BetaBarr,NumelPeaks,[0,1,2,3,4,5,6],'ShowText','on');
    hold on;
    c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
    clabel(M,c,'FontSize',20,'FontWeight','bold');
    %see where the flip line falls
    [M,c]=contour(BetaAarr,BetaBarr,SwitchA,[0,1]);
    c.LineColor='r';
    c.LineWidth=4;
    
    caxis([0,6]);
    %colormap(parula);
    colorbar;
    xlabel('Infection rate of $$N_{a}$$ individuals, $$\beta_{a}$$','interpreter','latex');
    ylabel('Infection rate of $$N_{b}$$ individuals, $$\beta_{b}$$','interpreter','latex');
    if m==25
        title(['Low sensitivity, $$m=$$',num2str(m)],'interpreter','latex');
    elseif m==50
        title(['Medium sensitivity, $$m=$$',num2str(m)],'interpreter','latex');
    else
        title(['High sensitivity, $$m=$$',num2str(m)],'interpreter','latex');
    end
    set(gca,'FontSize',25);

    figure(fig_c+1);[M,c]=contourf(BetaAarr,BetaBarr,SwitchA,[0,1],'ShowText','on');
    hold on;
    c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
    clabel(M,c,'FontSize',20,'FontWeight','bold');
    caxis([0,1]);
    %colormap(parula);
    colorbar;
    xlabel('$$\beta_{a}$$','interpreter','latex','FontWeight','bold');
    ylabel('$$\beta_{b}$$','interpreter','latex','FontWeight','bold');
    title(['Switch to opinion $$a$$, $$m=$$',num2str(m)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);
    
    figure(7);
    [M,c]=contour(BetaAarr,BetaBarr,SwitchA,[0,1]);
    c.LineColor=col(count,:);
    c.LineWidth=4;
    hold on;
    count=count+1;
    fig_c=fig_c+2;
end