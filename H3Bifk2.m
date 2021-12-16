% produce Holling III bifurcation diagram
% prepare state space

close all;
clear variables;
clc;
format long;

thetaA=5;
thetaB=5;
mesh_size=200;

k=2.7;
theta=5;
pB=0.4;
pAarr=linspace(0.1,1,mesh_size);
Res=[];
Resst=[];
Resunst1=[];
Resunst2=[];
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
           Res=[Res;pA a(ind)];
           fl=1;
       elseif numel(ind)==3
           Resst=[Resst; pA a(ind(2))];
           Resunst1=[Resunst1; pA a(ind(1))];
           Resunst2=[Resunst2; pA a(ind(3))];
           fl=1;
       else
          error('NA: could not obtain correct number of co-existence equilibria');
       end           
           
    end
end

figure(2);plot(Res(:,1)/pB,Res(:,2),'o');hold on;

green=[208, 240, 192]/255;
melon=[251, 206, 177]/255;
lblue=[173, 216, 230]/255;

%plot
figure(1);

h1=fill([pAarr fliplr(pAarr)]/pB,[pAarr*0 pAarr*0+1],'c');hold on;
set(h1,'facecolor',melon) %blue
set(h1,'EdgeColor','none')

h1=fill([pAarr fliplr(pAarr)]/pB,[(Res(:,2))' pAarr*0+1],'c');hold on;
set(h1,'facecolor',lblue) %blue
set(h1,'EdgeColor','none')


h1=plot(pAarr/pB,pAarr*0,'bo');
set(h1,'markerFacecolor',get(h1,'color'));hold on;
h2=plot(pAarr/pB,pAarr*0+1,'bo');
set(h2,'markerFacecolor',get(h1,'color'));
%overplot
h1=plot(pAarr/pB,pAarr*0,'bo');
set(h1,'markerFacecolor',get(h1,'color'));hold on;
h2=plot(pAarr/pB,pAarr*0+1,'bo');
set(h2,'markerFacecolor',get(h1,'color'));hold on;
h1=plot(Res(:,1)/pB,Res(:,2),'ro');
set(h1,'markerFacecolor',get(h1,'color'));hold on;

xlabel('$$p_{a}/p_{b}$$','interpreter','latex');
ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','interpreter','latex');
title(['$$k=',num2str(k),'$$, $$\theta_{a}=\theta_{b}=5$$, H3'],'interpreter','latex');
set(gca,'FontSize',25);
xlim([min(pAarr)/pB max(pAarr)/pB]);
