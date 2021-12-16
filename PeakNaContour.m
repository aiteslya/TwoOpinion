%this code generates contour map of peak n_{a} in different regions of c-m
%for different assortativity degrees

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
calc_err=5e-11;
fig_c=1;
%%
%c-omega
%set up the parameters
pA0=0.4;
pA1=1;
pB=0.4;
k=1.6;
thetaA=5;
thetaB=5;
betaA=0.8;
betaB=2;
gammaA=1;
gammaB=1;

num=30;
%set up of the parameter array
carr=linspace(1,105,num);
marr=linspace(5,100,num);

omegaArr=[0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.9,0.92,0.95,0.97];%linspace(0.05,0.95,10);

[Marr,Carr]=meshgrid(marr,carr);
figc=1;
for omega=omegaArr
    ResnA=zeros(num,num);
    Resinfect=zeros(num,num);
    NumelPeaks=zeros(num,num);
    for i1=1:1:num
        for i2=1:1:num
           m=Marr(i1,i2);
           C=Carr(i1,i2);
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
           peaknA=max(na);
           ResnA(i1,i2)=peaknA;
           peakInfect=max(infect);
           Resinfect(i1,i2)=peakInfect;
           [pks,locs] = findpeaks(infect,t);
           NumelPeaks(i1,i2)=numel(find(pks>6e-8));
        end
    end

    figure(fig_c);[M,c]=contourf(Marr,Carr,ResnA,[0,0.6,0.7,0.8,0.9,1],'ShowText','on');
    hold on;
    c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
    clabel(M,c,'FontSize',25,'FontWeight','bold');
    %highlight R0=1
    v = [1,0.9];
    [M,c]=contour(Marr,Carr,ResnA,v,'r','ShowText','on');
    c.LineWidth=6;
    clabel(M,c,'FontSize',25,'FontWeight','bold');
    caxis([0,1])
    colormap(parula);
    %colorbar;
    xlabel('$$m$$','interpreter','latex','FontWeight','bold');
    ylabel('$$c$$','interpreter','latex','FontWeight','bold');
    title(['Peak $$n_{a}$$, $$\omega=$$',num2str(omega)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);
    
    figure(fig_c+1);[M,c]=contourf(Marr,Carr,Resinfect,[0,0.005,0.01,0.02,0.03],'ShowText','on');
    hold on;
    c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
    clabel(M,c,'FontSize',25,'FontWeight','bold');
    caxis([0,0.08])
    colormap(parula);
    %colorbar;
    xlabel('$$m$$','interpreter','latex','FontWeight','bold');
    ylabel('$$c$$','interpreter','latex','FontWeight','bold');
    title(['Peak $$i$$, $$\omega=$$',num2str(omega)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);
    
    figure(fig_c+2);[M,c]=contourf(Marr,Carr,NumelPeaks,[0,1,2,3,4,5],'ShowText','on');
    hold on;
    c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
    clabel(M,c,'FontSize',25,'FontWeight','bold');
    %caxis([0.6,1.21]);
    colormap(parula);
    %colorbar;
    xlabel('$$m$$','interpreter','latex','FontWeight','bold');
    ylabel('$$c$$','interpreter','latex','FontWeight','bold');
    title(['Number of peaks of $$i$$, $$\omega=$$',num2str(omega)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);

    fig_c=fig_c+3;
end