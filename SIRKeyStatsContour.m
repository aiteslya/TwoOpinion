%this code generates contour maps for key statistiscs for a system which couples an SIR disease
%with competition of two opinions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear variables;
clc;
format long;

%pre-set colors
meshsize=20;
marr=linspace(20,120,meshsize);

%set up integration options
Atol=1e-12;
RelTol=1e-10;
opts = odeset('RelTol',RelTol,'AbsTol',Atol);
calc_err=5e-11;
%%
%set up the parameters
pAarr=linspace(0.3,0.55,meshsize);
omegaArr=[0.1,0.9];
pB=0.4;
pA1=1;%1;
k=1.6;
thetaA=5;
thetaB=5;
betaA=0.8;
betaB=2;
Carr=[25,50,105];
gammaA=1;
gammaB=1;

%set up the parameters
T=500;

[PA0arr,Marr]=meshgrid(pAarr,marr);
figc=1;

for C=Carr
    for omega=omegaArr
         % initialize the containers
         NAs=zeros(meshsize,meshsize);
         MaxPeak=nan(meshsize,meshsize);
         MaxPeakTime=nan(meshsize,meshsize);
         NumelPeak=zeros(meshsize,meshsize);
         for i1=1:meshsize
            for i2=1:meshsize
               pA0=PA0arr(i1,i2);
               m=Marr(i1,i2);
               %find the initial sa
               parsNa=[k,thetaA,thetaB,pA0,pB];
               sa=NA(parsNa);
               sb=1-sa;
               init=[sa,0,0,sb-6e-8,6e-8,0];
               pars=[C,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];

               [t,y]=ode45(@(t,y)TwoOpAssSIR(t,y,pars),[0,T], init,opts);
               infect=y(:,2)+y(:,5);
               na=y(:,1)+y(:,2)+y(:,3);
               % see if the population switched to zero
               if abs(na(end)-1)<calc_err
                   NAs(i1,i2)=1;
               end
               [pks,locs] = findpeaks(infect,t);
               maxpeak=max(pks);
               if numel(maxpeak)>0
                ind=find(infect==maxpeak,1);
                lockmaxpeak=t(ind);
                MaxPeak(i1,i2)=maxpeak;
                MaxPeakTime(i1,i2)=lockmaxpeak;
               end
               NumelPeak(i1,i2)=numel(find(pks>1e-8));
            end
         end
        figure(figc);
        [M,c]=contourf(PA0arr,Marr,NAs,[0,1,2],'ShowText','off');
        hold on;
        %c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
        %clabel(M,c,'FontSize',25,'FontWeight','bold');
        colormap(summer);
        xlabel('$$p_{a}(0)$$','interpreter','latex','FontWeight','bold');
        ylabel('$$m$$','interpreter','latex','FontWeight','bold');
        title(['Switch to $$a$$ with $$\omega=$$',num2str(omega)],'interpreter','latex','FontWeight','bold');
        set(gca,'FontSize',25);

        figure(figc+1);
        [M,c]=contourf(PA0arr,Marr,NumelPeak,[0,1,2,3,4],'ShowText','off');
        hold on;
        %c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
        colormap(lines);
        %clabel(M,c,'FontSize',25,'FontWeight','bold');
        xlabel('$$p_{a}(0)$$','interpreter','latex','FontWeight','bold');
        ylabel('$$m$$','interpreter','latex','FontWeight','bold');
        title(['Number of peaks with $$\omega=$$',num2str(omega)],'interpreter','latex','FontWeight','bold');
        set(gca,'FontSize',25);

        figure(figc+2);
        [M,c]=contourf(PA0arr,Marr,MaxPeak,linspace(0.9*min(min(MaxPeak)),1.2*max(max(MaxPeak)),7),'ShowText','on');
        hold on;
        c.LevelList=round(c.LevelList,2);  %rounds levels to 3rd decimal place
        clabel(M,c,'FontSize',20,'FontWeight','bold');
        xlabel('$$p_{a}(0)$$','interpreter','latex','FontWeight','bold');
        ylabel('$$m$$','interpreter','latex','FontWeight','bold');
        title(['Size of the maximum peak with $$\omega=$$',num2str(omega)],'interpreter','latex','FontWeight','bold');
        set(gca,'FontSize',25);

        figure(figc+3);
        [M,c]=contourf(PA0arr,Marr,MaxPeakTime,linspace(0.9*min(min(MaxPeakTime)),1.1*max(max(MaxPeakTime)),7),'ShowText','on');
        hold on;
        c.LevelList=round(c.LevelList,1);  %rounds levels to 3rd decimal place
        clabel(M,c,'FontSize',20,'FontWeight','bold');
        xlabel('$$p_{a}(0)$$','interpreter','latex','FontWeight','bold');
        ylabel('$$m$$','interpreter','latex','FontWeight','bold');
        title(['Time of the maximum peak with $$\omega=$$',num2str(omega)],'interpreter','latex','FontWeight','bold');
        set(gca,'FontSize',25);

        figc=figc+4;
    end
end

%reset y-axes

