% This script generates bifurcation diagram tracing na component of the the
% co-existence equilibrium as pa changes
clc;
clear variables;
close all;
format long;

k=1.6;
theta=5;
pB=0.4;
pAnum=80;
pAarr=linspace(0.1,1,pAnum);
Res=[];
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
       elseif numel(ind)==1 | numel(ind)==3
           for j=1:numel(ind)
               Res=[Res;pA a(ind(j))];
           end
           fl=1;
       else
          error('NA: could not obtain correct number of co-existence equilibria');
       end           
           
    end
end

figure(1);plot(Res(:,1),Res(:,2),'o')