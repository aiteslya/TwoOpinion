%this code generates time series for a system which couples an SIR disease
%with competition of two opinions and shows mechasnism through which the
%jump to the pure a opinion equilibrium happens, for two different
%assortativity levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear variables;
clc;
format long;

%set up integration options
Atol=1e-12;
RelTol=1e-10;
opts = odeset('RelTol',RelTol,'AbsTol',Atol);
calc_err=5e-11;
fig_c=1;
%%
%set up the parameters
pA0=0.4;
pA1=1;
c=40;
omegaArr=[0.2,0.6];
marr=[25,50,75];
l=numel(marr);
darkgreen = [53,94,59]/255;
lightgreen =  [152, 251, 152]/255;
noadap=[0.9,0.13,0.13];
yearcol=[0.17,0.17,0.17];
colors_p = [linspace(lightgreen(1),darkgreen(1),l)', linspace(lightgreen(2),darkgreen(2),l)', linspace(lightgreen(3),darkgreen(3),l)'];

pB=0.4;
k=1.6;
thetaA=5;
thetaB=5;
theta=5;
betaA=0.8;
betaB=2;
gammaA=1;
gammaB=1;

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
figure(1);
subplot(2,3,3);plot(Res(:,1)./pB,Res(:,2),'ro');hold on
plot(ResSt(:,1)./pB,ResSt(:,2),'bo');hold on
tempArr=linspace(0.1,1,pAnum);
plot(tempArr./pB,tempArr*0,'bo');
plot(tempArr./pB,tempArr*0+1,'bo');
xlim([0.1/pB,1/pB]);
xlabel('$$p_{a}/p_{b}$$','interpreter','latex');
ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','interpreter','latex');
set(gca,'FontSize',22);

subplot(2,3,6);plot(Res(:,1)./pB,Res(:,2),'ro');hold on
plot(ResSt(:,1)./pB,ResSt(:,2),'bo');hold on
tempArr=linspace(0.1,1,pAnum);
plot(tempArr./pB,tempArr*0,'bo');
plot(tempArr./pB,tempArr*0+1,'bo');
xlim([0.1/pB,1/pB]);
xlabel('$$p_{a}/p_{b}$$','interpreter','latex');
ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','interpreter','latex');
set(gca,'FontSize',22);


%omega=0.1
%set up the parameters
T=100;
%find the initial sa
parsNa=[k,thetaA,thetaB,pA0,pB];
sa=NA(parsNa);
sb=1-sa;
init=[sa,0,0,sb-6e-8,6e-8,0];

counter=1;
for omega=omegaArr
    mcounter=1;
    for m=marr
        pars=[c,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
        [t,y]=ode45(@(t,y)TwoOpAssSIR(t,y,pars),[0,T], init,opts);
        infect=y(:,2)+y(:,5);
        na=y(:,1)+y(:,2)+y(:,3);
        pA=pA0+(pA0-pA1)*(m+1)*(1./(1+m*infect)-1)/m;
       
        figure(1);
        subplot(2,3,(counter-1)*3+3);
        plot(pA/pB,na,'color',colors_p(mcounter,:),'LineWidth',3);hold on;
        %figure(2);
        subplot(2,3,(counter-1)*3+1);
        plot(t,na,'LineWidth',4,'color',colors_p(mcounter,:));hold on;

        if mcounter==1
            %add unstable upper equilibrium
            ga=@(x)(pA0*x^(k-1))/(1+thetaA*x^k);
            gb=@(x)(pB*x^(k-1))/(1+thetaB*x^k);
    
            num=1e3;
            a=linspace(1e-5,1-1e-5,num);
            b=1-a;
    
            eqn=-arrayfun(gb,b)+arrayfun(ga,a);
            nulleqn=eqn(1:(num-1)).*eqn(2:num);
            ind=find(nulleqn<0);
            if numel(ind)==3
                figure(1);subplot(2,3,(counter-1)*3+1);plot(t,t*0+a(max(ind)),'or');
            else
                error('a problem');
            end
        end
        %figure(3);
        subplot(2,3,(counter-1)*3+2);
        %collect the legends only for the second row
        if counter==2
            h(mcounter)=plot(t,infect,'LineWidth',4,'color',colors_p(mcounter,:));hold on;
        else
            plot(t,infect,'LineWidth',4,'color',colors_p(mcounter,:));hold on;
        end 
        mcounter=mcounter+1;
    end
    counter=counter+1;
end    
%set up figure settings

%figure(2);
subplot(2,3,1);
xlabel('Time','Interpreter','latex');
ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','Interpreter','latex');
set(gca,'FontSize',22);
subplot(2,3,4);
xlabel('Time','Interpreter','latex');
ylabel('Density of $$N_{a}$$ population, $$n_{a}$$','Interpreter','latex');
set(gca,'FontSize',22);

%figure(3);
subplot(2,3,2);
xlabel('Time','Interpreter','latex');
ylabel({'Density of infected individuals, $$i$$'},'Interpreter','latex');
set(gca,'FontSize',22);

subplot(2,3,5);
xlabel('Time','Interpreter','latex');
ylabel({'Density of infected individuals, $$i$$'},'Interpreter','latex');
legend(h,['$$m=',num2str(marr(1)),'$$, low sensitivity'],['$$m=',num2str(marr(2)),'$$, medium sensitivity'],['$$m=',num2str(marr(3)),'$$, high sensitivity'],'interpreter','latex','location','southoutside','NumColumns',3);
set(gca,'FontSize',22);

%write annotations

dim1 = [0.01,0.99,0,0];
str1 = 'a';
annotation('textbox',dim1,'String',str1,'FitBoxToText','on','FontSize',25,'FontWeight','bold','EdgeColor','none');

dim2 = [0.1,0.99,0,0];
str2 = 'b';
annotation('textbox',dim2,'String',str2,'FitBoxToText','on','FontSize',25,'FontWeight','bold','EdgeColor','none');

dim3 = [0.2,0.99,0,0];
str3 = 'c';
annotation('textbox',dim3,'String',str3,'FitBoxToText','on','FontSize',25,'FontWeight','bold','EdgeColor','none');

dim4 = [0.01,0.49,0,0];
str4 = 'd';
annotation('textbox',dim4,'String',str4,'FitBoxToText','on','FontSize',25,'FontWeight','bold','EdgeColor','none');

dim5 = [0.1,0.49,0,0];
str5 = 'e';
annotation('textbox',dim5,'String',str5,'FitBoxToText','on','FontSize',25,'FontWeight','bold','EdgeColor','none');

dim6 = [0.2,0.49,0,0];
str6 = 'f';
annotation('textbox',dim6,'String',str6,'FitBoxToText','on','FontSize',25,'FontWeight','bold','EdgeColor','none');

%assortativity legends for the rows
dim7 = [0.01,0.79,0,0];
str7 = 'Degree of assortativity $$\omega=0.2$$';
annotation('textbox',dim7,'String',str7,'FitBoxToText','on','FontSize',25,'FontWeight','bold','EdgeColor','none','Interpreter','latex');

dim8 = [0.01,0.29,0,0];
str8 = 'Degree of assortativity $$\omega=0.6$$';
annotation('textbox',dim8,'String',str8,'FitBoxToText','on','FontSize',25,'FontWeight','bold','EdgeColor','none','Interpreter','latex');

