%This script generates time series for two opinion system coupled with SIR
%infectious disease. These time series are outputed as ease as a phase
%curves
clc;
clear variables;
close all;
format long;
%set up the parameters

thetaA=5;
thetaB=5;
betaA=0.8;
betaB=2;
gammaA=1;
gammaB=1;
k=1.6;
theta=5;
pB=0.4;
omega=0;

pA0arr=[0.35,0.4,0.5];
pA1arr=[0.55,0.8,1];
Carr=[5,10,20,40,80,100];

%generate time series
%pre-set colors
mnum=20;
marr=linspace(10,100,mnum);
% l=numel(marr);
% darkblue = [18, 10, 143]/255;%[0.01, 0.4, 0.76];
% lightblue =  [0.45,0.76,0.98];
% medblue=(lightblue+darkblue)/2;
% noadap=[0.9,0.13,0.13];
% yearcol=[0.17,0.17,0.17];
% colors_p = [linspace(lightblue(1),darkblue(1),l)', linspace(lightblue(2),darkblue(2),l)', linspace(lightblue(3),darkblue(3),l)'];
% colors_p=[noadap;colors_p];
% marr=[0,marr];

%set up integration options
Atol=1e-12;
RelTol=1e-10;
opts = odeset('RelTol',RelTol,'AbsTol',Atol);
calc_err=5e-11;
%%

%set up the parameters

figc=1;
for pA1=pA1arr
    for pA0=pA0arr
        ccounter=1;
        Res=zeros(numel(Carr),numel(marr));
        Resna=zeros(numel(Carr),numel(marr));
        for C=Carr
            mcounter=1;
            for m=marr
                %calculate initial conditions - stable co-existence of
                %opinions
               parsNa=[k,thetaA,thetaB,pA0,pB];
               sa=NA(parsNa);
               sb=1-sa;
               init=[sa,0,0,sb-6e-8,6e-8,0];
               if m==0
                   pars=[C,pA0,pA0,1,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
               else
                   pars=[C,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
               end
               T=500;
               fl=0;
               while ~fl
                [t,y]=ode45(@(t,y)TwoOpAssSIR(t,y,pars),[0,T], init,opts);
                infect=y(:,2)+y(:,5);
                if infect(end)<6e-8
                    fl=1;
                    cumul=y(:,2)+y(:,3)+y(:,5)+y(:,6);
                    na=y(:,1)+y(:,2)+y(:,3);
                    if abs(na(end)-1)<6e-8
                        Resna(ccounter,mcounter)=1;
                    else
                        Resna(ccounter,mcounter)=sa;
                    end
                    Res(ccounter,mcounter)=max(infect);
                else
                    T=T*2;
                end
               end
               mcounter=mcounter+1;
            end
            figure(figc);
            h1(ccounter)=plot(marr,Res(ccounter,:),'*');hold on;
            figure(figc+1);
            h2(ccounter)=plot(marr,Resna(ccounter,:),'*');hold on;
            ccounter=ccounter+1;
        end
        figure(figc);
        title(['$$p_{a}(0)/p_{b}=',num2str(pA0/pB),'$$ and $$p_{a}(1)/p_{b}=',num2str(pA1/pB),'$$'],'interpreter','latex');
        xlabel('$$m$$','interpreter','latex');
        ylabel('Size of the maximum peak of $$i$$','interpreter','latex');
        legend(h1,['$$c=',num2str(Carr(1)),'$$'],['$$c=',num2str(Carr(2)),'$$'],['$$c=',num2str(Carr(3)),'$$'],['$$c=',num2str(Carr(4)),'$$'],['$$c=',num2str(Carr(5)),'$$'],['$$c=',num2str(Carr(6)),'$$'],'interpreter','latex');
        set(gca,'FontSize',25);
        figure(figc+1);
        title(['$$p_{a}(0)/p_{b}=',num2str(pA0/pB),'$$ and $$p_{a}(1)/p_{b}=',num2str(pA1/pB),'$$'],'interpreter','latex');
        xlabel('$$m$$','interpreter','latex');
        ylabel('Density of $$N_{a}$$ individuals after the epidemic','interpreter','latex');
        legend(h2,['$$c=',num2str(Carr(1)),'$$'],['$$c=',num2str(Carr(2)),'$$'],['$$c=',num2str(Carr(3)),'$$'],['$$c=',num2str(Carr(4)),'$$'],['$$c=',num2str(Carr(5)),'$$'],['$$c=',num2str(Carr(6)),'$$'],'interpreter','latex');
        set(gca,'FontSize',25);
        figc=figc+2;
    end
end