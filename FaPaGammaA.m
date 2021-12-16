% script that produces fa as a function of na, pa and gamma a as a function
% of number of infected

%prepare space
close all;
clc;
clear variables; format long;

fig_c=1;
pA=0.5;
%plot fa for 2 different pairs of k and theta
thetakarr=[5 5; 1.6 2.7];
col_n=size(thetakarr,2);
n_points=50;
Na=linspace(1e-6,1-1e-6,n_points);
figure(fig_c);plot(Na,pA*Na,'-+');hold on;
for j=1:col_n
   theta=thetakarr(1,j);
   k=thetakarr(2,j);
   if j==1 
       fa=@(x)pA*x/(1+theta*x);
       figure(fig_c);h2=plot(Na,arrayfun(fa,Na),'-o');hold on;
   end
   fa=@(x)pA*x^(k)/(1+theta*x^k);
   figure(fig_c);h=plot(Na,arrayfun(fa,Na),'-*');hold on;
end
xlabel('Density of $$N_{a}$$ population, $$n_{a}$$','Interpreter','latex');
ylabel('Switch rate function $$f_{a}(n_{a})$$','Interpreter','latex');
legend('$$\theta=0$$ and $$k=1$$, H1',['$$\theta=',num2str(thetakarr(1,1)),'$$ and $$k=1$$, H2'],['$$\theta=',num2str(thetakarr(1,1)),'$$ and $$k=',num2str(thetakarr(2,1)),'$$, H3'],['$$\theta=',num2str(thetakarr(1,2)),'$$ and $$k=',num2str(thetakarr(2,2)),'$$, H3'],'interpreter','latex');
set(gca,'FontSize',25);
fig_c=fig_c+1;

%% map pa(i)

%preset colors
k_arr=[1,25,50,75];
l=numel(k_arr);
darkblue = [18, 10, 143]/255;%[0.01, 0.4, 0.76];
lightblue =  [0.45,0.76,0.98];
medblue=(lightblue+darkblue)/2;
novacc=[0.9,0.13,0.13];
yearcol=[0.17,0.17,0.17];
colors_p = [linspace(lightblue(1),darkblue(1),l)', linspace(lightblue(2),darkblue(2),l)', linspace(lightblue(3),darkblue(3),l)'];

infect=linspace(0,1,n_points);

pA0=0.5;
pA1=0.9;
k_count=1;
for k=k_arr
    B=(pA1-pA0)*(1+1/k);
    pA=@(x)pA0+B*(1-1/(k*x+1));
    figure(fig_c);plot(infect,arrayfun(pA,infect),'LineWidth',3,'color',colors_p(k_count,:));hold on;
    k_count=k_count+1;
end
xlabel('Density of infectected individuals, $$i$$','interpreter','latex');
ylabel({'Probability of switching to';'opinion $$a$$ per contact, $$p_{a}$$'},'interpreter','latex');
legend(['$$m=',num2str(k_arr(1)),'$$'],['$$m=',num2str(k_arr(2)),'$$'],['$$m=',num2str(k_arr(3)),'$$'],['$$m=',num2str(k_arr(4)),'$$'],'interpreter','latex');
set(gca,'FontSize',25);
fig_c=fig_c+1;

