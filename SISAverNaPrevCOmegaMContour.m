%this code generates contour map of average n_{a} and average prevalence in
% different regions of omega-c space

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
betaB=1.5;%2;
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

num=80;%45;
%set up of the parameter array
carr=linspace(1,105,num);
marr=[25,50,75];

omegaArr=linspace(0.05,0.95,num);

[Omegaarr,Carr]=meshgrid(omegaArr,carr);
figc=1;
for m=marr
    AvernA=zeros(num,num);
    AverPrev=zeros(num,num);
    for i1=1:1:num
        for i2=1:1:num
           omega=Omegaarr(i1,i2);
           C=Carr(i1,i2);
           %need to calculate sa and sb
           parsNa=[k,thetaA,thetaB,pA0,pB];
           sa=NA(parsNa);
           sb=1-sa;
           T=1000;
           trans=2*T/3;
           fl=0;
           init=[sa,0,sb-6e-8,6e-8];
           
           pars=[C,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
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

    figure(fig_c);
    surf(Omegaarr,Carr,AvernA);
    view(2);
    xlim([0.5,0.95])
    ylim([1,105])
    grid off

    colormap(parula);
    caxis([0.6,1]);
    colorbar;
    xlabel('Degree of assortativity, $$\omega$$','interpreter','latex','FontWeight','bold');
    ylabel('Social contact rate, $$c$$','interpreter','latex','FontWeight','bold');
    title(['Sensitivity of reaction, $$m=$$',num2str(m,2)],'interpreter','latex','FontWeight','bold');

    set(gca,'FontSize',30);
    
     figure(fig_c+1);
     surf(Omegaarr,Carr,AverPrev);
     xlim([0.5,0.95])
     ylim([1,105])
     view(2);

    grid off
    colormap(parula);
    caxis([0,0.017]);
    colorbar;
    %colorbar;
    xlabel('Degree of assortativity, $$\omega$$','interpreter','latex','FontWeight','bold');
    ylabel('Social contact rate, $$c$$','interpreter','latex','FontWeight','bold');
    title(['Sensitivity of reaction, $$m=$$',num2str(m,2)],'interpreter','latex','FontWeight','bold');

    set(gca,'FontSize',30);
  
    fig_c=fig_c+2;
 
end