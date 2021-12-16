%this script produces parameter sweeps for an SIS coupled with opinion
%competition searching for sustained oscillations in betA-betaB-m space
%prepare space
clc;
clear variables;
close all;
format long;

%set parameters
pB=0.4;
meshsize=10;
pA0arr=linspace(0.283,0.56,meshsize);
pA1arr=linspace(0.6,1,meshsize);
marr=10*2.^(linspace(0,6.7,meshsize));
gammaA=1;
gammaB=1;
k=1.6;
thetaA=5;
thetaB=5;
betaAarr=linspace(1e-3,0.6,meshsize);
betaBarr=linspace(2,7,meshsize);
Carr=10*2.^(linspace(0,6.7,meshsize));

omega=0;
% set up the integration settings
Atol=1e-12;
RelTol=1e-10;
opts = odeset('RelTol',RelTol,'AbsTol',Atol);
%tolerance for calculation of the equilibrium
calc_err=1e-8;
T=3e3;
fig_c=1;
trans=2*T/3;
%need to calculate sa and sb
for pA0=pA0arr
    parsNa=[k,thetaA,thetaB,pA0,pB];
    sa=NA(parsNa);
    sb=1-sa;
    init=[sa,0,sb-6e-8,6e-8];
    for pA1=pA1arr
        for C=Carr
            for m=marr
                for betaA=betaAarr
                   for betaB=betaBarr
                    pars=[C,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
                    [t,y]=ode45(@(t,y)TwoOpAssSIS(t,y,pars),[0,T], init,opts);
                    infect=y(:,2)+y(:,4);
                    ind=find(t>trans,1);
                    ttrim=t(ind:end);
                    infecttrim=infect(ind:end);
                    [pks,locs] = findpeaks(infecttrim,ttrim);
                    NumelPeaks=numel(find(pks>6e-8));

                    if NumelPeaks>0 & max(pks)-min(infecttrim)>calc_err
                        disp(['Potential sustained oscillations for pA0=',num2str(pA0),', pA1=',num2str(pA1),', C=',num2str(C)]);
                        disp(['m=',num2str(m),', betaA=',num2str(betaA),' and betaB=',num2str(betaB)]);
                        disp(['Figure ',num2str(fig_c)]);
                        if mod(fig_c,20)==1 & fig_c>1
                            close(fig_c-1:-1:fig_c-20);
                        end
                        figure(fig_c);
                        plot(ttrim,infecttrim);
                        fig_c=fig_c+1;
                    end
                   end
                end   
            end
   
        end
    end
end
