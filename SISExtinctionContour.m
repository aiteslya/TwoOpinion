%this script produces contour maps which map  regions of pa(1) and m for
%which the infection goes extinct when R0>1, for the disease which follows
%SIS dynamics coupled with the competition of two opinions, where the
%population reacts to the growth of prevalence

%prepare space
clc;
clear variables;
close all;
format long;

%set parameters
pB=0.4;
pA0=pB;
meshsize=60;
marr=linspace(1,100,meshsize);
pA1arr=linspace(pB,1,meshsize);
gammaA=1;
gammaB=1;
k=1.6;
thetaA=5;
thetaB=5;
betaA=0.8;
betaB=1.5;
C=40;

[PA1arr,Marr]=meshgrid(pA1arr,marr);
omegaarr=[0.7,0.8,0.9];
% set up the integration settings
Atol=1e-12;
RelTol=1e-10;
opts = odeset('RelTol',RelTol,'AbsTol',Atol);
%tolerance for calculation of the equilibrium
calc_err=1e-8;
%need to calculate sa and sb
parsNa=[k,thetaA,thetaB,pA0,pB];
sa=NA(parsNa);
sb=1-sa;
T=4000;
init=[sa,0,sb-6e-8,6e-8];
fig_c=1;
for omega=omegaarr
    %allocate result containers:
    AvernA=nan(meshsize,meshsize);
    AverPrev=nan(meshsize,meshsize);
    for i1=1:meshsize
        for i2=1:meshsize
           pA1=PA1arr(i1,i2);
           m=Marr(i1,i2);
           pars=[C,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
           T=4000;
           trans=2*T/3;
           fl=0;
           init=[sa,0,sb-6e-8,6e-8];
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
                if max(abs(max(infecttrim)-min(infecttrim)),abs(max(natrim)-min(natrim)))<calc_err % settled to an equilibrium
                    if max(infecttrim)<calc_err
                        AverPrev(i1,i2)=0;
                    else
                        AverPrev(i1,i2)=max(infecttrim);
                    end
                    AvernA(i1,i2)=max(natrim);
                    fl=1;
                else
                    %check for periodic orbit
                    [pks,locs] = findpeaks(infecttrim,ttrim);
                    indpks=find(pks>calc_err);
                    pks=pks(indpks);
                    locs=locs(indpks);
                    if numel(pks)>2
                        avpeaks=mean(pks);
                        ampl=(pks-min(infecttrim))/2;
                        meanampl=mean(ampl);
                        if max(abs(ampl-meanampl))/meanampl<1e-5% relative error<calc_err
                            fl=1;
                            %find period
                            per=locs(3)-locs(2);
                            ind1=find(ttrim==locs(2),1);
                            ind2=find(ttrim==locs(3),1);
                            %find average over the period
                            AverPrev(i1,i2)=trapz(ttrim(ind1:ind2),infecttrim(ind1:ind2))/per;
                            AvernA(i1,i2)=trapz(ttrim(ind1:ind2),natrim(ind1:ind2))/per;
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
                    AverPrev(i1,i2)=0;
                    AvernA(i1,i2)=1;
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
    figure(1);subplot(1,3,fig_c);
    surf(PA1arr,Marr,AvernA);
    view(2);
    caxis([0.5,1]);
    colormap(flipud(parula))
    %colormap(parula);
    if fig_c==2
        colB=colorbar('southoutside');
        colB.Label.String = 'Average density of $$N_{a}$$ population';
        colB.Label.Interpreter='latex';
    end
    xlabel('$$p_{a}(1)$$','interpreter','latex','FontWeight','bold');
    ylabel('Sensitivity of reaction, $$m$$','interpreter','latex','FontWeight','bold');
    title(['Assortativity degree $$\omega=$$',num2str(omega)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);

    figure(2);subplot(1,3,fig_c);
    surf(PA1arr,Marr,AverPrev);
    view(2);

    if fig_c==2
        colB=colorbar('southoutside');
        colB.Label.String = 'Average prevalence of infected';
        colB.Label.Interpreter='latex';
    end
    xlabel('$$p_{a}(1)$$','interpreter','latex','FontWeight','bold');
    ylabel('Sensitivity of reaction, $$m$$','interpreter','latex','FontWeight','bold');
    title(['Assortativity degree $$\omega=$$',num2str(omega)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',25);

    fig_c=fig_c+1;
end
