%Pure opinion model, homogeneous contact rates
%investigation of the number of equilibria as p_{a} changes and switching
%rates assume different shapes
close all; 
clear variables;
clc;
format long;
%set up parameters
pB=0.1;
figCount=1;
panum=200;
% Holling 1

countRes=1;
countBi=1;
Res=[];
Resbi=[];
Resmin=[];
Resmax=[];
ResMax=[];
ResMin=[];
countMM=1;
pAarr1=linspace(0.1*pB,pB,panum/2);
pAarr2=linspace(pB,2*pB,panum/2);

figure(figCount);
h1=fill([pAarr1 fliplr(pAarr1)]/pB,[pAarr1*0+1 pAarr1*0],'c');hold on
h2=fill([pAarr2 fliplr(pAarr2)]/pB,[pAarr2*0+1 pAarr2*0],'m');

set(h1,'EdgeColor','none')
set(h2,'EdgeColor','none')

plot(pAarr1/pB,pAarr1*0+1,'ro','MarkerSize',6,'MarkerFaceColor','r');hold on;
plot(pAarr1/pB,pAarr1*0,'bo','MarkerSize',6,'MarkerFaceColor','b');hold on;
plot(pAarr2/pB,pAarr2*0+1,'bo','MarkerSize',6,'MarkerFaceColor','b');hold on;
plot(pAarr2/pB,pAarr2*0,'ro','MarkerSize',6,'MarkerFaceColor','r');hold on;
vline=linspace(0,1,100);
xvline=vline*0+1;
plot(xvline,vline,'bo','MarkerSize',6,'MarkerFaceColor','b');hold on;
xlim([min(pAarr1)/pB,max(pAarr2)/pB]);

xlabel('p_{a}/p_{b}','FontWeight','bold');
ylabel('n_{a}','FontWeight','bold');
title('k=1, $$\theta=0$$, H1','Interpreter','latex');
%title(['n=',num2str(n)],'interpreter','tex','FontWeight','bold');
box on;
%grid on;
set(gca,'fontsize',25);
clear b;
clear expr;
clear nullexpr;
clear nullind;
figCount=figCount+1;

% Holling 2
n=1;
theta=5;
thetaA=theta;
thetaB=theta;

countRes=1;
countBi=1;
Res=[];
Resbi=[];
Resmin=[];
Resmax=[];
ResMax=[];
ResMin=[];
countMM=1;

pA1=pB/(1+thetaB);
pA2=pB*(1+thetaA);
pAarr1=linspace(0.1*pB,pA1*1.001,50);
pAarr=linspace(pA1*1.001,1,100);

for pA=pAarr
    fl=0;
    anum=2e6;
    a=linspace(1e-6,0.999999,anum);
    b=1-a;
    expr=-pB.*(1+thetaA*a.^n)+pA.*(1+thetaB*b.^n);
    nullexpr=expr(1:(anum-1)).*expr(2:anum);
    nullind=setdiff(find(nullexpr<0),[1,anum]);
    if numel(nullind)>0
        ResMax(1,countMM)=max(a(nullind));
        ResMin(1,countMM)=min(a(nullind));
        countMM=countMM+1;
        for count=1:1:numel(nullind)
            Res(1,countRes)=pA;
            Res(2,countRes)=a(nullind(count));
            countRes=countRes+1;
        end
        if numel(nullind)>2
            Resbi(1,countBi)=pA;
            Resbi(2,countBi)=a(nullind(2));
            Resmin(1,countBi)=pA;
            Resmin(2,countBi)=a(nullind(1));
            Resmax(1,countBi)=pA;
            Resmax(2,countBi)=a(nullind(3));
            countBi=countBi+1;
        elseif numel(nullind)==2
            error('Two equilibria found')
        end
    end
end
if numel(Resbi)>0
    figure(figCount);
    A1=find(pAarr==Resbi(1,1),1);
    A2=find(pAarr==Resbi(1,end),1);
    A2Res=find(Res(1,:)==Resbi(1,end),1);
%            
    l=size(Res,2);
    h1=fill([pAarr(1:A2) pAarr(A2) pAarr(A2+1:end) fliplr(pAarr)]/pB,[ResMax(1:A2) ResMin(A2) ResMax(A2+1:end) pAarr*0+1],'c');hold on
    h2=fill([pAarr(1:A1) pAarr(A1) pAarr(A1:end) fliplr(pAarr)]/pB,[ResMax(1:A1) ResMax(A1) ResMin(A1:end) pAarr*0],'m');
    h3=fill([pAarr(A1:A2) fliplr(pAarr(A1:A2))]/pB,[ResMax(A1:A2) fliplr(ResMin(A1:A2))],'g');
    set(h1,'EdgeColor','none')
    set(h2,'EdgeColor','none')
    set(h3,'EdgeColor','none')
    plot(Res(1,:)/pB,Res(2,:),'ro','MarkerSize',6,'MarkerFaceColor','r');hold on;
    plot(pAarr/pB,pAarr*0,'bo','MarkerSize',6,'MarkerFaceColor','b');hold on;
    plot(pAarr/pB,pAarr*0+1,'bo','MarkerSize',6,'MarkerFaceColor','b');hold on;
    xlim([min(pAarr)/pB,max(pAarr)/pB]);
    plot(Resbi(1,:)/pB,Resbi(2,:),'bo','MarkerSize',6,'MarkerFaceColor','b');hold on;
else
    figure(figCount);
    h1=fill([pAarr1 fliplr(pAarr1)]/pB,[pAarr1*0+1 pAarr1*0],'c');hold on
    h3=fill([pAarr fliplr(pAarr)]/pB,[pAarr*0+1 pAarr*0],'g');hold on
    set(h1,'EdgeColor','none')
    set(h3,'EdgeColor','none')
    
    plot(pAarr1/pB,pAarr1*0+1,'ro','MarkerSize',6,'MarkerFaceColor','r');hold on;
    plot(pAarr/pB,pAarr*0,'ro','MarkerSize',6,'MarkerFaceColor','r');hold on;
    plot(pAarr/pB,pAarr*0+1,'ro','MarkerSize',6,'MarkerFaceColor','r');hold on;
    plot(Res(1,:)/pB,Res(2,:),'bo','MarkerSize',6,'MarkerFaceColor','b');hold on;
    plot(pAarr1/pB,pAarr1*0,'bo','MarkerSize',6,'MarkerFaceColor','b');hold on;
    xlim([min(pAarr1)/pB,max(pAarr)/pB]);
end
xlabel('p_{a}/p_{b}','FontWeight','bold');
ylabel('n_{a}','FontWeight','bold');
title('k=1, $$\theta_{a}=\theta_{b}=5$$, H2','Interpreter','latex');
box on;
%grid on;
set(gca,'fontsize',25);
clear b;
clear expr;
clear nullexpr;
clear nullind;
figCount=figCount+1;
%Holling 3

narr=linspace(1.2,2,4);%[1.5,2.5,4];

pBarr=[0.2,0.5,0.8];

for pB=0.4%pBarr
    for n=1.6%narr

        countRes=1;
        countBi=1;
        Res=[];
        Resbi=[];
        Resmin=[];
        Resmax=[];
        ResMax=[];
        ResMin=[];
        countMM=1;
        
        pAarr=linspace(0.1*pB,2*pB,panum);
        for pA=pAarr
            fl=0;
            anum=2e6;
            a=linspace(1e-9,0.999999999,anum);
            b=1-a;
            expr=-pB*((b).^(n-1)).*(1+thetaA*a.^n)+pA*((a).^(n-1)).*(1+thetaB*b.^n);
            nullexpr=expr(1:(anum-1)).*expr(2:anum);
            nullind=setdiff(find(nullexpr<0),[1,anum]);
            ResMax(1,countMM)=max(a(nullind));
            ResMin(1,countMM)=min(a(nullind));
            countMM=countMM+1;
            for count=1:1:numel(nullind)
                Res(1,countRes)=pA;
                Res(2,countRes)=a(nullind(count));
                countRes=countRes+1;
            end
            if numel(nullind)>2
                Resbi(1,countBi)=pA;
                Resbi(2,countBi)=a(nullind(2));
                Resmin(1,countBi)=pA;
                Resmin(2,countBi)=a(nullind(1));
                Resmax(1,countBi)=pA;
                Resmax(2,countBi)=a(nullind(3));
                countBi=countBi+1;
            elseif numel(nullind)==2
                error('Two equilibria found')
            end
        end
        if numel(Resbi)>0
            figure(figCount);
            A1=find(pAarr==Resbi(1,1),1);
            A2=find(pAarr==Resbi(1,end),1);
            A2Res=find(Res(1,:)==Resbi(1,end),1);
%            
            l=size(Res,2);
            h1=fill([pAarr(1:A2) pAarr(A2) pAarr(A2+1:end) fliplr(pAarr)]/pB,[ResMax(1:A2) ResMin(A2) ResMax(A2+1:end) pAarr*0+1],'c');hold on
            h2=fill([pAarr(1:A1) pAarr(A1) pAarr(A1:end) fliplr(pAarr)]/pB,[ResMax(1:A1) ResMax(A1) ResMin(A1:end) pAarr*0],'m');
            h3=fill([pAarr(A1:A2) fliplr(pAarr(A1:A2))]/pB,[ResMax(A1:A2) fliplr(ResMin(A1:A2))],'g');
            set(h1,'EdgeColor','none')
            set(h2,'EdgeColor','none')
            set(h3,'EdgeColor','none')
            plot(Res(1,:)/pB,Res(2,:),'ro','MarkerSize',6,'MarkerFaceColor','r');hold on;
            plot(pAarr/pB,pAarr*0,'bo','MarkerSize',6,'MarkerFaceColor','b');hold on;
            plot(pAarr/pB,pAarr*0+1,'bo','MarkerSize',6,'MarkerFaceColor','b');hold on;
            xlim([min(pAarr)/pB,max(pAarr)/pB]);
            plot(Resbi(1,:)/pB,Resbi(2,:),'bo','MarkerSize',6,'MarkerFaceColor','b');hold on;
        else
            figure(figCount);
            h1=fill([pAarr fliplr(pAarr)]/pB,[Res(2,:)*0+1 fliplr(Res(2,:))],'c');hold on
            h2=fill([pAarr fliplr(pAarr)]/pB,[Res(2,:)*0 fliplr(Res(2,:))],'m');
            
            set(h1,'EdgeColor','none')
            set(h2,'EdgeColor','none')
                        
            plot(Res(1,:)/pB,Res(2,:),'ro','MarkerSize',6,'MarkerFaceColor','r');hold on;
            plot(pAarr/pB,pAarr*0,'bo','MarkerSize',6,'MarkerFaceColor','b');hold on;
            plot(pAarr/pB,pAarr*0+1,'bo','MarkerSize',6,'MarkerFaceColor','b');hold on;
            xlim([min(pAarr)/pB,max(pAarr)/pB]);
        end
        xlabel('p_{a}/p_{b}','FontWeight','bold');
        ylabel('n_{a}','FontWeight','bold');
        title(['k=',num2str(n)],'Interpreter','latex');
        %title(['n=',num2str(n)],'interpreter','tex','FontWeight','bold');
        box on;
        %grid on;
        set(gca,'fontsize',25);
        clear b;
        clear expr;
        clear nullexpr;
        clear nullind;
        figCount=figCount+1;
    end
end