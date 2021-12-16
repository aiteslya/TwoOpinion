function dydt=TwoOpAss(t,y,pars)
    % RHS of the system that keeps track of competition between two
    % opposing health opinions
    dydt=zeros(1,2);
    na=y(1);
    nb=y(2);
    %pars=[c,pA,pB,thetaA,thetaB,k,omega];
    c=pars(1);
    pA=pars(2);
    pB=pars(3);
    thetaA=pars(4);
    thetaB=pars(5);
    k=pars(6);
    omega=pars(7);
    
    fa=@(x)pA*x^k/(1+thetaA*x^k);
    fb=@(x)pB*x^k/(1+thetaB*x^k);
    %na compartment
    dydt(1)=-na*c*(1-omega)*fb(nb)+nb*c*(1-omega)*fa(na);    
    %nb compartment
    dydt(2)=na*c*(1-omega)*fb(nb)-nb*c*(1-omega)*fa(na);
    dydt=dydt';
    
end