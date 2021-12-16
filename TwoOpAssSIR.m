function dydt=TwoOpAssSIR(t,y,pars)
    % RHS of the system that keeps track of the opinion at the time of the
    % infection
    dydt=zeros(1,6);
    sa=y(1);
    ia=y(2);
    ra=y(3);
    sb=y(4);
    ib=y(5);
    rb=y(6);
    na=sa+ia+ra;
    nb=sb+ib+rb;
    infec=ia+ib;
    %pars=[c,pA0,pA1,m,pB,thetaA,thetaB,k,betaA,betaB,gammaA,gammaB,omega];
    c=pars(1);
    pA0=pars(2);
    pA1=pars(3);
    m=pars(4);
    pB=pars(5);
    thetaA=pars(6);
    thetaB=pars(7);
    k=pars(8);
    betaA=pars(9);
    betaB=pars(10);
    gammaA=pars(11);
    gammaB=pars(12);
    omega=pars(13);
    pA=pA0+(pA0-pA1)*(m+1)*(1./(1+m*infec)-1)/m;
    fa=@(x)pA*x^k/(1+thetaA*x^k);
    fb=@(x)pB*x^k/(1+thetaB*x^k);
    %na compartment
    if nb>1e-13
        dydt(1)=-c*sa*(1-omega)*fb(nb)+c*sb*(1-omega)*fa(na)-sa*betaA*(omega*ia/na+(1-omega)*(ia+ib));
        dydt(2)=-c*ia*(1-omega)*fb(nb)+c*ib*(1-omega)*fa(na)+sa*betaA*(omega*ia/na+(1-omega)*(ia+ib))-gammaA*ia;
        dydt(3)=-c*ra*(1-omega)*fb(nb)+c*rb*(1-omega)*fa(na)+gammaA*ia;
        dydt(4)=c*sa*(1-omega)*fb(nb)-c*sb*(1-omega)*fa(na)-sb*betaB*(omega*ib/nb+(1-omega)*(ia+ib));
        dydt(5)=c*ia*(1-omega)*fb(nb)-c*ib*(1-omega)*fa(na)+sb*betaB*(omega*ib/nb+(1-omega)*(ia+ib))-gammaB*ib;
        dydt(6)=c*ra*(1-omega)*fb(nb)-c*rb*(1-omega)*fa(na)+gammaB*ib;
    else
        dydt(1)=-sa*betaA*ia;
        dydt(2)=sa*betaA*ia-gammaA*ia;
        dydt(3)=gammaA*ia;
        dydt(4)=0;
        dydt(5)=0;
        dydt(6)=0;
    end
    
    %nb compartment
    
    dydt=dydt';
    
end