%This script generates time series for two opinion system coupled with SIR
%infectious disease. These time series are outputed as ease as a phase
%curves
clc;
clear variables;
close all;
format long;
%generate bifurcation diagram on which things will be overimposed.
k=1.6;
theta=5;
pB=0.4;
pAnum=400;
pAarr=sort([pB,linspace(0.1,1,pAnum)]);
Res=[];
ResSt=[];
for pA=pAarr
    ga=@(x)(pA*x^(k-1))/(1+theta*x^k);
    gb=@(x)(pB*x^(k-1))/(1+theta*x^k);
    
    num=1e3;
    a=linspace(1e-5,1-1e-5,num);
    b=1-a;
    fl=0;
    while ~fl
       eqn=-arrayfun(gb,b)+arrayfun(ga,a);
       nulleqn=eqn(1:(num-1)).*eqn(2:num);
       ind=find(nulleqn<0);
       if numel(ind)==0
           error('NA: check the parameters provided, no co-existence eq found');
       elseif numel(ind)==2 & num<1e6
           num=2*num;
           a=linspace(1e-5,1-1e-5,num);
           b=1-a;
       elseif numel(ind)==1 
           Res=[Res;pA a(ind(1))];
           fl=1;
       elseif numel(ind)==3
           Res=[Res;pA a(ind(1))];
           Res=[Res;pA a(ind(3))];
           ResSt=[ResSt; pA a(ind(2))];
           fl=1;
       else
          error('NA: could not obtain correct number of co-existence equilibria');
       end           
           
    end
end
figure(1);plot(Res(:,1)./pB,Res(:,2),'ro');hold on
plot(ResSt(:,1)./pB,ResSt(:,2),'bo');hold on
xlabel('$$p_{a}(0)/p_{b}$$','interpreter','latex');
ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','interpreter','latex');
set(gca,'FontSize',25);

IND=find(Res(:,1)==pB);
ind=IND(2)

%generate time series
%pre-set colors
marr=[20,25,50,75,80];
l=numel(marr);
darkblue = [18, 10, 143]/255;%[0.01, 0.4, 0.76];
lightblue =  [0.45,0.76,0.98];
medblue=(lightblue+darkblue)/2;
noadap=[0.9,0.13,0.13];
yearcol=[0.17,0.17,0.17];
colors_p = [linspace(lightblue(1),darkblue(1),l)', linspace(lightblue(2),darkblue(2),l)', linspace(lightblue(3),darkblue(3),l)'];
colors_p=[noadap;colors_p];
marr=[0,marr];
% one set omega=0.1 and pA=pB %double peak - maybe we can find better
%another set omega=0.1 and pA=pB=1.25
pA1=1;
C=105;
%set up integration options
Atol=1e-12;
RelTol=1e-10;
opts = odeset('RelTol',RelTol,'AbsTol',Atol);
calc_err=5e-11;
%%
%set up the parameters
pA0=pB;
omega=0.9;
k=1.6;
thetaA=5;
thetaB=5;
betaA=0.8;
betaB=2;
c=105;
gammaA=1;
gammaB=1;

%set up the parameters
T=500;
mcounter=1;
for m=marr
   parsNa=[k,thetaA,thetaB,pA0,pB];
   sa=NA(parsNa);
   sb=1-sa;
   init=[sa,0,0,sb-6e-8,6e-8,0];
   if m==0
       pars=[C,pA0,pA0,1,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
   else
       pars=[C,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
   end
   [t,y]=ode45(@(t,y)TwoOpAssSIR(t,y,pars),[0,T], init,opts);
   infect=y(:,2)+y(:,5);
   na=y(:,1)+y(:,2)+y(:,3);
   if m>0
    pA=pA0+(pA0-pA1)*(m+1)*(1./(1+m*infect)-1)/m;
   end
   figure(1);
   plot(pA/pB,na,'color',colors_p(mcounter,:));
   figure(2);
   plot(t,infect,'color',colors_p(mcounter,:),'LineWidth',4);hold on;
   figure(3);
   plot(t,na,'color',colors_p(mcounter,:),'LineWidth',4);hold on;
   figure(4);
   plot(t,pA/pB,'color',colors_p(mcounter,:),'LineWidth',4);hold on;
   mcounter=mcounter+1;
end
figure(3);plot(t,t*0+Res(ind,2),'g');