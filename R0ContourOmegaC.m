%this code generates  R0 for permutation of the following parameters:
% c and omega: 1 figure
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
%c-omega
%set up the parameters
num=60;
pB=0.4;
pAarr=[0.8,1,1.25]*pB;

k=1.6;
thetaA=5;
thetaB=5;
betaA=0.8;
betaB=1.5;
gammaA=1;
gammaB=1;

carr=linspace(1,105,num);
omegaArr=linspace(0.05,0.95,num);%linspace(0.05,0.95,num);

[Omegaarr,Carr]=meshgrid(omegaArr,carr);
str='abc';
for pA=pAarr
    %heterogeneous contacts with assortativity
    Res=zeros(num,num);
    Rescheck=zeros(num,num);
    for i1=1:1:num
        for i2=1:1:num
           omega=Omegaarr(i1,i2);
           C=Carr(i1,i2);
           %need to calculate sa and sb
           parsNa=[k,thetaA,thetaB,pA,pB];
           sa=NA(parsNa);
           sb=1-sa;
           parsR0=[betaA,betaB,omega,pA,pB,k,thetaA,thetaB,gammaA,gammaB,C];
           r0=R0(sa,parsR0);
           fa=@(x)pA*x^k/(1+thetaA*x^k);
           fb=@(x)pB*x^k/(1+thetaB*x^k);
           numer=betaA*sa*(1-omega)*(gammaB+C*fa(sa)*(1-omega)+C*fb(sb)*(1-omega))+betaB*sb*(1-omega)*(gammaA+C*fb(sb)*(1-omega)+C*fa(sa)*(1-omega))-betaA*betaB*(1-omega)*omega;
           den=(gammaA-betaA*omega)*(gammaB-betaB*omega)+C*(1-omega)*(fa(sa)*(gammaA-betaA*omega)+fb(sb)*(gammaB-betaB*omega));
           R0check=numer/den;
           Res(i1,i2)=r0;
           Rescheck(i1,i2)=R0check;
           %add here investigation for R0a: it will be marked with dark blue
           %line on the screen
        end
    end

    figure(fig_c);[M,c]=contourf(Omegaarr,Carr,Res,linspace(0.95*min(min(Res)),1.05*max(max(Res)),10),'ShowText','on');
    %figure(fig_c);[M,c]=contourf(Omegaarr,Carr,Res);
    hold on;

    c.LineWidth=1;
    c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
    clabel(M,c,'FontSize',35,'FontWeight','bold');
    %highlight R0=1
    v = [1,1];
    [M,c]=contour(Omegaarr,Carr,Res,v,'r','ShowText','on');
    c.LineWidth=6;
    clabel(M,c,'FontSize',35,'FontWeight','bold');
    caxis([0.9,1.5]);
    %colormap(spring);
    if fig_c==3
        colorbar;
    end
    xlabel('Degree of assortativity, $$\omega$$','interpreter','latex','FontWeight','bold');
    ylabel('Social contact rate, $$c$$','interpreter','latex','FontWeight','bold');
    title(['$$p_{a}/p_{b}=',num2str(pA/pB),'$$'],'interpreter','latex');
    figure(fig_c);annotation('textbox', [0.03, 0.999, 0, 0], 'string', str((fig_c+1)/2),'FontWeight','bold','FontSize',35)
    set(gca,'FontSize',35);
    disp(['Max R0 is ',num2str(max(max(Res))),'and Min R0 is ',num2str(min(min(Res)))]);

    fig_c=fig_c+2;
end