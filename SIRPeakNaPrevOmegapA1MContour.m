%this code generates contour map of peak n_{a} and peak prevalence in
% different regions of omega-c space for SIR disease

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
calc_err=5e-10;
fig_c=1;
%%
%c-omega
%set up the parameters
pA0=0.4;

pB=0.4;
k=1.6;
thetaA=5;
thetaB=5;
betaA=0.8;
betaB=2;
gammaA=1;
gammaB=1;
C=40;

% find all co-existence equilibria
ga=@(x)(pA0*x^(k-1))/(1+thetaA*x^k);
gb=@(x)(pB*x^(k-1))/(1+thetaB*x^k);

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
   elseif numel(ind)==1 | numel(ind)==3
       a(ind)
       fl=1;
   else
      error('NA: could not obtain correct number of co-existence equilibria');
   end           

end

num=45;
%set up of the parameter array
pA1arr=linspace(0.4,1,num);
marr=linspace(1,100,num);

omegaArr=[0,0.9];

%[Omegaarr,Carr]=meshgrid(omegaArr,carr);
[PA1arr,Marr]=meshgrid(pA1arr,marr);
figc=1;
for omega=omegaArr
    ResnA=zeros(num,num);
    Resinfect=zeros(num,num);
    for i1=1:1:num
        for i2=1:1:num
%            omega=Omegaarr(i1,i2);
%            C=Carr(i1,i2);
            pA1=PA1arr(i1,i2);
            m=Marr(i1,i2);
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
            na=y(:,1)+y(:,2)+y(:,3);
            peaknA=max(na);
            ResnA(i1,i2)=peaknA;
            peakInfect=max(infect);
            Resinfect(i1,i2)=peakInfect;
           end
        end
    end
    figure(1);
    subplot(1,2,fig_c)
    surf(PA1arr/pB,Marr,ResnA);
    view(2);
    xlim([1,1/pB])
    ylim([1,100])
    grid off

    colormap(parula);
    caxis([0.5,1]);
    colorbar;
    xlabel('$$p_{a}(1)/p_{b}$$','interpreter','latex','FontWeight','bold');
    ylabel('Sensitivity of reaction, $$m$$','interpreter','latex','FontWeight','bold');
    title(['Assortativity degree $$\omega=$$',num2str(omega,2)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',30);
    
    figure(2);
    subplot(1,2,fig_c)
    surf(PA1arr/pB,Marr,Resinfect);
    xlim([1,1/pB])
    ylim([1,100])
    view(2);
    grid off
    colormap(parula);
    caxis([0.004,0.08]);
    colorbar;
    xlabel('$$p_{a}(1)/p_{b}$$','interpreter','latex','FontWeight','bold');
    ylabel('Sensitivity of reaction, $$m$$','interpreter','latex','FontWeight','bold');
    title(['Assortativity degree $$\omega=$$',num2str(omega,2)],'interpreter','latex','FontWeight','bold');

    set(gca,'FontSize',30);
  
    fig_c=fig_c+1;
 
end