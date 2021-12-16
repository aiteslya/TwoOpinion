function dydt=TwoOpAssSIRE(t,y,pars)
    % RHS of the system that keeps track of the opinion at the time of the
    % infection
    dydt=zeros(1,10);
    sa=y(1);
    ia=y(2);
    iba=y(3);%had b opinion but is currently in a
    ra=y(4);
    rba=y(5);
    Ia=ia+iba;
    sb=y(6);
    ib=y(7);
    iab=y(8);% had a opinion but is currently in b
    rb=y(9);
    rab=y(10);
    Ib=ib+iab;
    na=sa+ia+iba+ra+rba;
    nb=sb+ib+iab+rb+rab;
    infec=ia+ib+iba+iab;
    %pars=[ca0,cA1,C,cb,pa0,pa1,k1,pb,thetaA, thetaB,n, Pa,Pb, gammaA, gammaB,omega];
    cA0=pars(1);
    cA1=pars(2);
    C=pars(3);
    cB=pars(4);
    pA0=pars(5);
    pA1=pars(6);
    k1=pars(7);
    pB=pars(8);
    thetaA=pars(9);
    thetaB=pars(10);
    n=pars(11);
    Pa=pars(12);
    Pb=pars(13);
    gammaA=pars(14);
    gammaB=pars(15);
    omega=pars(16);
    cA=cA0+(cA0-cA1)*(C+1)*(1./(1+C*infec)-1)/C;
    pA=pA0+(pA0-pA1)*(k1+1)*(1./(1+k1*infec)-1)/k1;
    fa=@(x)pA*x^n/(1+thetaA*x^n);
    fb=@(x)pB*x^n/(1+thetaB*x^n);
    %na compartment
    dydt(1)=-cA*sa*(1-omega)*fb(cB*nb/(cA*na+cB*nb))+cB*sb*(1-omega)*fa(cA*na/(cA*na+cB*nb))-sa*cA*Pa*(omega*Ia/na+(1-omega)*((cA*Ia+cB*Ib)/(cA*na+cB*nb)));
    dydt(2)=-cA*ia*(1-omega)*fb(cB*nb/(cA*na+cB*nb))+cB*iab*(1-omega)*fa(cA*na/(cA*na+cB*nb))+sa*cA*Pa*(omega*Ia/na+(1-omega)*((cA*Ia+cB*Ib)/(cA*na+cB*nb)))-gammaA*ia;
    dydt(3)=cB*ib*(1-omega)*fa(cA*na/(cA*na+cB*nb))-cA*iba*(1-omega)*fb(cB*nb/(cA*na+cB*nb))-gammaA*iba;%iba
    dydt(4)=-cA*ra*(1-omega)*fb(cB*nb/(cA*na+cB*nb))+cB*rab*(1-omega)*fa(cA*na/(cA*na+cB*nb))+gammaA*ia;
    dydt(5)=cB*rb*(1-omega)*fa(cA*na/(cA*na+cB*nb))-cA*rba*(1-omega)*fb(cB*nb/(cA*na+cB*nb))+gammaA*iba;%rba
    %nb compartment
    dydt(6)=cA*sa*(1-omega)*fb(cB*nb/(cA*na+cB*nb))-cB*sb*(1-omega)*fa(cA*na/(cA*na+cB*nb))-sb*cB*Pb*(omega*Ib/nb+(1-omega)*((cA*Ia+cB*Ib)/(cA*na+cB*nb)));
    dydt(7)=cA*iba*(1-omega)*fb(cB*nb/(cA*na+cB*nb))-cB*ib*(1-omega)*fa(cA*na/(cA*na+cB*nb))+sb*cB*Pb*(omega*Ib/nb+(1-omega)*((cA*Ia+cB*Ib)/(cA*na+cB*nb)))-gammaB*ib;
    dydt(8)=-cB*iab*(1-omega)*fa(cA*na/(cA*na+cB*nb))+cA*ia*(1-omega)*fb(cB*nb/(cA*na+cB*nb))-gammaB*iab;%iab
    dydt(9)=cA*rba*(1-omega)*fb(cB*nb/(cA*na+cB*nb))-cB*rb*(1-omega)*fa(cA*na/(cA*na+cB*nb))+gammaB*ib;
    dydt(10)=-cB*rab*(1-omega)*fa(cA*na/(cA*na+cB*nb))+cA*ra*(1-omega)*fb(cB*nb/(cA*na+cB*nb))+gammaB*iab;%rab
    dydt=dydt';
    
end