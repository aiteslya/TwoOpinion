function r0=R0(na,pars)
% function takes initial distribution of opinions and respective parameter
% set up and returns the basic reproductive number
% pars=[betaA,betaB,omega,pA,pB,k,thetaA,thetaB,gammaA,gammaB,c];
    betaA=pars(1);
    betaB=pars(2);
    omega=pars(3);
    pA=pars(4);
    pB=pars(5);
    k=pars(6);
    thetaA=pars(7);
    thetaB=pars(8);
    gammaA=pars(9);
    gammaB=pars(10);
    c=pars(11);
    nb=1-na;
    F(1,1)=na*betaA*(omega/na+1-omega);
    F(1,2)=na*betaA*(1-omega);
    F(2,1)=nb*betaB*(1-omega);
    F(2,2)=nb*betaB*(omega/nb+1-omega);

    fa=@(x)pA*x^k/(1+thetaA*x^k);
    fb=@(x)pB*x^k/(1+thetaB*x^k);
    V(1,1)=c*(1-omega)*fb(nb)+gammaA;
    V(1,2)=-c*(1-omega)*fa(na);
    V(2,1)=-c*(1-omega)*fb(nb);
    V(2,2)=c*(1-omega)*fa(na)+gammaB;
    Vinv=inv(V);
    FV=F*Vinv;
    e=eig(FV);
    r0=max(e);
end