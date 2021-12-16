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
meshsize=50;
pA0=0.283;
pA1=0.6;
C=10;

gammaA=1;
gammaB=1;
k=1.6;
thetaA=5;
thetaB=5;
marr=[50,75];
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
    NA=zeros(meshsize,meshsize);
    for i1=1:meshsize
        for i2=1:meshsize
            betaA=BetaAarr(i1,i2);
            betaB=BetaBarr(i1,i2);
                
            pars=[C,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
            fl=0;
            T=1e4;
            trans=2*T/3;
            while ~fl
               if T>2e6
                   disp('the integration time is too long');
               end
            [t,y]=ode45(@(t,y)TwoOpAssSIS(t,y,pars),[0,T], init,opts);
            infect=y(:,2)+y(:,4);
            na=y(:,1)+y(:,2);
            if t(end)>=T
                ind=find(t>trans,1);
                ttrim=t(ind:end);
                infecttrim=infect(ind:end);
                natrim=na(ind:end);
                %figure(100);plot(ttrim,infecttrim);
                if max(abs(max(infecttrim)-min(infecttrim)),abs(max(natrim)-min(natrim)))<calc_err % settled to an equilibrium
                    if max(infecttrim)<calc_err
                       MeanInfect(i1,i2)=0;
                    else
                       MeanInfect(i1,i2)=max(infecttrim);
                    end
                    if abs(max(natrim)-1)<calc_err
                        NA(i1,i2)=1;
                    else
                        NA(i1,i2)=max(natrim);
                    end
                    fl=1;
                else
                    %check for periodic orbit
                    [pks,locs] = findpeaks(infecttrim,ttrim);
                    indpks=find(pks>calc_err);
                    pks=pks(indpks);
                    locs=locs(indpks);
                    if numel(pks)>2 % 
                        avpeaks=mean(pks);
                        ampl=(pks-min(infecttrim))/2;
                        meanampl=mean(ampl);
                        if max(abs(ampl-meanampl))/meanampl<1e-4% relative error<calc_err % located a periodic orbit
                            fl=1;
                            %find period
                            per=locs(3)-locs(2);
                            Per(i1,i2)=per;
                            ind1=find(ttrim==locs(2),1);
                            ind2=find(ttrim==locs(3),1);
                            %find average over the period
                            MeanInfect(i1,i2)=trapz(ttrim(ind1:ind2),infecttrim(ind1:ind2))/per;
                            if abs(max(natrim)-1)<calc_err
                                NA(i1,i2)=1;
                            else
                                NA(i1,i2)=trapz(ttrim(ind1:ind2),natrim(ind1:ind2))/per;
                            end
                            %find amplitude
                            Ampl(i1,i2)=max(infecttrim)-MeanInfect(i1,i2);
                        else % did not finish integrating
                            T=2*T;
                            trans=2*T/3;
                            init=y(end,:);
                        end
                    else
                        T=2*T;
                        trans=2*T/3;
                        init=y(end,:);
                    end
                end
            else
                if abs(1-natrim(end))<calc_err | abs(y(end,1)-1)<calc_err
                    MeanInfect(i1,i2)=0;
                    NA(i1,i2)=1;
                    fl=1;
                else
                    figure(100);plot(t,infect);hold on;
                    figure(101);plot(t,na); hold on;
                    disp('Did not integrate until the end');
                end
            end
           end
   end
end 
    %output the results;
    figure(1);
    subplot(1,numel(marr),fig_c);
    s=surf(BetaAarr,BetaBarr,Per);
    s.EdgeColor = 'none';
    hold on;
    view(2);
    
    %see where the flip line falls
    [M,c]=contour(BetaAarr,BetaBarr,NA,[0,1]);
    c.LineColor='r';
    c.LineWidth=4;   
    if fig_c==2
        colorbar;
        %caxis([0.99*min(min(Per)),1.01*max(max(Per))]);
    end
    
    xlabel({'Infection rate of $$N_{a}$$ individuals, $$\beta_{a}$$'},'interpreter','latex','FontWeight','bold');
    ylabel({'Infection rate of $$N_{b}$$ individuals, $$\beta_{b}$$'},'interpreter','latex','FontWeight','bold');
    title(['Sensitivity of reaction $$m=',num2str(m),'$$'],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);
    
    figure(2);
    subplot(1,numel(marr),fig_c);
    s=surf(BetaAarr,BetaBarr,Ampl);
    s.EdgeColor = 'none';
    hold on;
    view(2);
   
    %see where the flip line falls
    [M,c]=contour(BetaAarr,BetaBarr,NA,[0,1]);
    c.LineColor='r';
    c.LineWidth=4;    
    
    if fig_c==2
        colorbar;
        
    end
    xlabel({'Infection rate of $$N_{a}$$ individuals, $$\beta_{a}$$'},'interpreter','latex','FontWeight','bold');
    ylabel({'Infection rate of $$N_{b}$$ individuals, $$\beta_{b}$$'},'interpreter','latex','FontWeight','bold');
    title(['Sensitivity of reaction $$m=',num2str(m),'$$'],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);
    
    figure(3);
    subplot(1,numel(marr),fig_c)
    s=surf(BetaAarr,BetaBarr,MeanInfect);
    s.EdgeColor = 'none';
    hold on;
    view(2);
    
    %see where the flip line falls
    [M,c]=contour(BetaAarr,BetaBarr,NA,[0,1]);
    c.LineColor='r';
    c.LineWidth=4;    
    
    if fig_c==2
        colorbar;
        %axis([0.99*min(min(MeanInfect)),1.01*max(max(MeanInfect))]);
    end
    
    xlabel({'Infection rate of $$N_{a}$$ individuals, $$\beta_{a}$$'},'interpreter','latex','FontWeight','bold');
    ylabel({'Infection rate of $$N_{b}$$ individuals, $$\beta_{b}$$'},'interpreter','latex','FontWeight','bold');
    title(['Sensitivity of reaction $$m=',num2str(m),'$$'],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);

    figure(4);
    subplot(1,numel(marr),fig_c)
    s=surf(BetaAarr,BetaBarr,NA);
    s.EdgeColor = 'none';
    hold on;
    view(2);
    
    %see where the flip line falls
    [M,c]=contour(BetaAarr,BetaBarr,NA,[0,1]);
    c.LineColor='r';
    c.LineWidth=4;    
    
    if fig_c==2
        colorbar;
        
    end
    
    xlabel({'Infection rate of $$N_{a}$$ individuals, $$\beta_{a}$$'},'interpreter','latex','FontWeight','bold');
    ylabel({'Infection rate of $$N_{b}$$ individuals, $$\beta_{b}$$'},'interpreter','latex','FontWeight','bold');
    title(['Sensitivity of reaction $$m=',num2str(m),'$$'],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);

    fig_c=fig_c+1;
end