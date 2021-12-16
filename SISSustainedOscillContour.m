%this script produces contour plots of existence of periodic solutions,
%their amplitudes and period for an SIS coupled with opinion
%competition searching for sustained oscillations in betA-betaB-m space
%prepare space
clc;
clear variables;
close all;
format long;

%set parameters
pB=0.4;
meshsize=20;
pA0=0.283;
pA1=0.6;
C=10;%%16
marr=[45,65];
gammaA=1;
gammaB=1;
k=1.6;
thetaA=5;
thetaB=5;
betaAarr=linspace(0.2,0.6,meshsize);
betaBarr=linspace(5,7,meshsize);
[BetaAarr,BetaBarr]=meshgrid(betaAarr,betaBarr);
omega=0;
% set up the integration settings
Atol=1e-12;
RelTol=1e-10;
opts = odeset('RelTol',RelTol,'AbsTol',Atol);
%tolerance for calculation of the equilibrium
calc_err=1e-7;
T=1e4;
fig_c=1;
trans=2*T/3;
%need to calculate sa and sb to use as init

parsNa=[k,thetaA,thetaB,pA0,pB];
sa=NA(parsNa);
sb=1-sa;
init=[sa,0,sb-6e-8,6e-8];

for m=marr
    %initialize containers for the results
    % period of the solution
    Per=zeros(meshsize,meshsize);
    %amplitude of the solution
    Ampl=zeros(meshsize,meshsize);
    %infected
    MeanInfect=zeros(meshsize,meshsize);
    % switch to pure a solution
    NA=zeros(meshsize,meshsize);
    for i1=1:meshsize
        for i2=1:meshsize
            betaA=BetaAarr(i1,i2);
            betaB=BetaBarr(i1,i2);
            pars=[C,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
            [t,y]=ode45(@(t,y)TwoOpAssSIS(t,y,pars),[0,T], init,opts);
            infect=y(:,2)+y(:,4);
            
            ind=find(t>trans,1);
            ttrim=t(ind:end);
            infecttrim=infect(ind:end);
            [pks,locs] = findpeaks(infecttrim,ttrim);
            NumelPeaks=numel(find(pks>6e-8));
            MeanInfect(i1,i2)=mean(infecttrim);
            if NumelPeaks>0 & max(pks)-min(infecttrim)>calc_err
                %further analys
                avpeaks=mean(pks);
                ampl=(pks-min(infecttrim))/2;
                meanampl=mean(ampl);
                if max(abs(ampl-meanampl))<calc_err % single peak per cycle
                    Ampl(i1,i2)=meanampl;
                    %calculate period
                    Per(i1,i2)=locs(3)-locs(2);
                    MeanInfect(i1,i2)=min(infecttrim)+meanampl;
                else
                    figure(100);
                    plot(ttrim,infecttrim);hold on;
                    disp('Did not integrate long enough or several local maxima per cycle');
                     disp(['pA0=',num2str(pA0),', pA1=',num2str(pA1),', C=',num2str(C)]);
                        disp(['m=',num2str(m),', betaA=',num2str(betaA),' and betaB=',num2str(betaB)]);
                end
            else %check if switched to na population
                na=y(:,1)+y(:,2);
                natrim=na(ind:end);
                if abs((natrim(end)-1))<calc_err
                    NA(i1,i2)=1;
                end
            end
       end
    end 
    %output the results;
    figure(fig_c);
    surf(BetaAarr,BetaBarr,Per);
    hold on;
    view(2);
    %c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
    %clabel(M,c,'FontSize',20,'FontWeight','bold');
    %see where the flip line falls
    [M,c]=contour(BetaAarr,BetaBarr,NA,[0,1]);
    c.LineColor='r';
    c.LineWidth=4;   
    colorbar;
    caxis([0.99*min(min(Per)),1.01*max(max(Per))]);
    xlabel('Infection rate of $$N_{a}$$ individuals, $$\beta_{a}$$','interpreter','latex','FontWeight','bold');
    ylabel('Infection rate of $$N_{b}$$ individuals, $$\beta_{b}$$','interpreter','latex','FontWeight','bold');
    title(['Period of the cycle, $$m=$$',num2str(m)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);
    
    figure(fig_c+1);
    surf(BetaAarr,BetaBarr,Ampl);
    hold on;
    view(2);
    
    %see where the flip line falls
    [M,c]=contour(BetaAarr,BetaBarr,NA,[0,1]);
    c.LineColor='r';
    c.LineWidth=4;    
    
    colorbar;
    caxis([0.99*min(min(Ampl)),1.01*max(max(Ampl))]);
    xlabel('Infection rate of $$N_{a}$$ individuals, $$\beta_{a}$$','interpreter','latex','FontWeight','bold');
    ylabel('Infection rate of $$N_{b}$$ individuals, $$\beta_{b}$$','interpreter','latex','FontWeight','bold');
    title(['Amplitude of the cycle, $$m=$$',num2str(m)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);
    
    figure(fig_c+2);
    surf(BetaAarr,BetaBarr,MeanInfect);
    hold on;
    view(2);
    
    %see where the flip line falls
    [M,c]=contour(BetaAarr,BetaBarr,NA,[0,1]);
    c.LineColor='r';
    c.LineWidth=4;    
       
    colorbar;
    caxis([0.99*min(min(MeanInfect)),1.01*max(max(MeanInfect))]);
    xlabel('Infection rate of $$N_{a}$$ individuals, $$\beta_{a}$$','interpreter','latex','FontWeight','bold');
    ylabel('Infection rate of $$N_{b}$$ individuals, $$\beta_{b}$$','interpreter','latex','FontWeight','bold');
    title(['Average prevalence $$i$$, $$m=$$',num2str(m)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);
    fig_c=fig_c+3;
end