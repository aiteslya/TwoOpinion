% produce Holling III bifurcation diagram
% prepare state space

close all;
clear variables;
clc;
format long;

thetaA=5;
thetaB=5;
mesh_size=200;

k=1.6;
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
figure(2);plot(Resunst1(:,1)/pB,Resunst1(:,2),'o');hold on;
figure(2);plot(Resst(:,1)/pB,Resst(:,2),'o');hold on;
figure(2);plot(Resunst2(:,1)/pB,Resunst2(:,2),'o');hold on;

indpA1=find(pAarr==min(Resst(:,1)),1);
indpA2=find(pAarr==max(Resst(:,1)),1);
indpA2Res=find(Res(:,1)>max(Resst(:,1)),1);

na1=pAarr(1,1:indpA1)./pB;
na2=pAarr(1,indpA1:indpA2,1)./pB;
na3=pAarr(1,indpA2:end)./pB;

green=[208, 240, 192]/255;
melon=[251, 206, 177]/255;
lblue=[173, 216, 230]/255;

%plot
figure(1);

h1=fill([na1 fliplr(na1)],[na1*0 na1*0+1],'c');hold on;
set(h1,'facecolor',melon) %blue
set(h1,'EdgeColor','none')

h1=fill([na1(1:(numel(na1)-1)) min(na2) min(na2) fliplr(na1(1:(numel(na1)-1)))],[na1(1:(numel(na1)-1))*0+1 1 Resunst2(1,2) fliplr((Res(1:indpA1-1,2))')],'c');hold on;
%h1=fill([na1 min(na2) min(na2) fliplr(na1)],[na1*0+1 1 Resunst2(1,2) fliplr((Res(1:indpA1-1,2))')],'c');hold on;
set(h1,'facecolor',lblue)
set(h1,'EdgeColor','none')

h1=fill([na2 fliplr(na2)],[na2*0 na2*0+1],'c');hold on;
set(h1,'facecolor',green) % green
set(h1,'EdgeColor','none')

h1=fill([na2 fliplr(na2)],[na2*0+1 fliplr(Resunst2(:,2)')],'c');hold on;
set(h1,'facecolor',lblue)
set(h1,'EdgeColor','none')

h1=fill([na2 fliplr(na2)],[na2*0 fliplr(Resunst1(:,2)')],'c');hold on;
set(h1,'facecolor',melon)
set(h1,'EdgeColor','none')

h1=fill([na3 fliplr(na3)],[na3*0 na3*0+1],'c');hold on;
set(h1,'facecolor',lblue) %blue
set(h1,'EdgeColor','none')

h1=fill([Resunst1(end,1) na3(2:end) max(na2) max(na2) fliplr(na3(2:end)) Resunst1(end,1)],[Resunst1(end,2) (Res(indpA2Res:end,2))' 0 0 na3(2:end)*0 0],'c');hold on;
set(h1,'facecolor',melon)
set(h1,'EdgeColor','none')

h1=plot(pAarr/pB,pAarr*0,'bo');
set(h1,'markerFacecolor',get(h1,'color'));hold on;
h2=plot(pAarr/pB,pAarr*0+1,'bo');
set(h2,'markerFacecolor',get(h1,'color'));
%overplot
h1=plot(na1,na1*0,'bo');
set(h1,'markerFacecolor',get(h1,'color'));hold on;
h2=plot(na2,na2*0+1,'bo');
set(h2,'markerFacecolor',get(h1,'color'));hold on;
h1=plot(Res(:,1)/pB,Res(:,2),'ro');
set(h1,'markerFacecolor',get(h1,'color'));hold on;

h1=plot(Resunst1(:,1)/pB,Resunst1(:,2),'ro');
set(h1,'markerFacecolor',get(h1,'color'));hold on;
h1=plot(Resunst2(:,1)/pB,Resunst2(:,2),'ro');
set(h1,'markerFacecolor',get(h1,'color'));hold on;
h1=plot(Resst(:,1)/pB,Resst(:,2),'bo');
set(h1,'markerFacecolor',get(h1,'color'));hold on;

xlabel('$$p_{a}/p_{b}$$','interpreter','latex');
ylabel('Density of $$N_{a}$$ population,$$n_{a}$$','interpreter','latex');
title('$$k=1.6$$, $$\theta_{a}=\theta_{b}=5$$, H3','interpreter','latex');
set(gca,'FontSize',25);
xlim([min(pAarr)/pB max(pAarr)/pB]);
