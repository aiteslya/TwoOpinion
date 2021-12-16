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
pA1=1;
pB=0.4;
k=1.6;
thetaA=5;
thetaB=5;
betaA=0.8;
betaB=2;
gammaA=1;
gammaB=1;

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
carr=linspace(1,105,num);
marr=[25,50,75];

omegaArr=linspace(0.05,0.95,num);

[Omegaarr,Carr]=meshgrid(omegaArr,carr);
figc=1;
for m=marr
    ResnA=zeros(num,num);
    Resinfect=zeros(num,num);
    for i1=1:1:num
        for i2=1:1:num
           omega=Omegaarr(i1,i2);
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
            na=y(:,1)+y(:,2)+y(:,3);
            peaknA=max(na);
            ResnA(i1,i2)=peaknA;
            peakInfect=max(infect);
            Resinfect(i1,i2)=peakInfect;
           end
        end
    end
    figure(fig_c);
    surf(Omegaarr,Carr,ResnA);
    view(2);
    xlim([0.05,0.95])
    ylim([1,105])
    grid off

    colormap(flipud(parula))
    caxis([0.5,1]);
    colorbar;
    xlabel('Degree of assortativity, $$\omega$$','interpreter','latex','FontWeight','bold');
    ylabel('Social contact rate, $$c$$','interpreter','latex','FontWeight','bold');
    title(['Sensitivity of reaction $$m=$$',num2str(m,2)],'interpreter','latex','FontWeight','bold');
    set(gca,'FontSize',30);
    
     figure(fig_c+1);
     surf(Omegaarr,Carr,Resinfect);
     xlim([0.5,0.95])
     ylim([1,105])
     view(2);
    grid off
    colormap(parula);
    caxis([0,0.04]);
    colorbar;
    xlabel('Degree of assortativity, $$\omega$$','interpreter','latex','FontWeight','bold');
    ylabel('Social contact rate, $$c$$','interpreter','latex','FontWeight','bold');
    title(['Sensitivity of reaction $$m=$$',num2str(m,2)],'interpreter','latex','FontWeight','bold');

    set(gca,'FontSize',30);
  
    fig_c=fig_c+2;
 
end