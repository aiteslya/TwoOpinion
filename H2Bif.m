% produce Holling II bifurcation diagram
% prepare state space
 
clear variables;
clc;
format long;

pB=0.1;
thetaA=5;
thetaB=5;
mesh_size=200;
%find pA for which there is no stable co-existence
pA1=pB/(1+thetaB);
pA2=pB*(1+thetaA);
pAarr=linspace(0.01*pA1,1.10*pA2,mesh_size);

num=1e3;
na=linspace(1e-6,1-1e-6,num);
nb=1-na;

Res=[];
for pA=pAarr
   eqn=-pB./(1+thetaB*nb)+pA./(1+thetaA*na);
   nulleqn=eqn(1:(num-1)).*eqn(2:num);
   ind=find(nulleqn<0);
   if numel(ind)>0
       Res=[Res;pA na(ind)];
   end
end



na1=linspace(0.01*pA1/pB,pA1/pB,10);
na2=linspace(pA2/pB,1.10*pA2/pB,20);
na3=linspace(pA1/pB,pA2/pB,20);

%plot
figure(1);

h1=fill([na1 fliplr(na1)],[na1*0 na1*0+1],'c');hold on;
set(h1,'facecolor',[251, 206, 177]/255)
set(h1,'EdgeColor','none')

h1=fill([na2 fliplr(na2)],[na2*0 na2*0+1],'c');hold on;
set(h1,'facecolor',[208, 240, 192]/255)
set(h1,'EdgeColor','none')

h1=fill([na3 fliplr(na3)],[na3*0 na3*0+1],'c');hold on;
set(h1,'facecolor',[173, 216, 230]/255)
set(h1,'EdgeColor','none')



h1=plot(pAarr/pB,pAarr*0,'ro');
set(h1,'markerFacecolor',get(h1,'color'));hold on;
h2=plot(pAarr/pB,pAarr*0+1,'ro');
set(h2,'markerFacecolor',get(h1,'color'));
%overplot
h1=plot(na1,na1*0,'bo');
set(h1,'markerFacecolor',get(h1,'color'));hold on;
h2=plot(na2,na2*0+1,'bo');
set(h2,'markerFacecolor',get(h1,'color'));hold on;
h1=plot(Res(:,1)/pB,Res(:,2),'bo');
set(h1,'markerFacecolor',get(h1,'color'));hold on;
xlabel('$$p_{a}/p_{b}$$','interpreter','latex');
ylabel('Density of $$N_{a}$$ population,$$n_{a}$$','interpreter','latex');
title(['$$k=1$$, $$\theta_{a}=\theta_{b}=5$$, H2'],'interpreter','latex');
set(gca,'FontSize',25);
xlim([0.01*pA1/pB,1.10*pA2/pB]);
