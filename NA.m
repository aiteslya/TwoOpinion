function na=NA(pars)
    %this function returns stable co-existence equilibrium given that the
    %switch rate function is Holling type III
    %parsNa=[k,thetaA,thetaB,pA,pB];
    k=pars(1);
    thetaA=pars(2);
    thetaB=pars(3);
    pA=pars(4);
    pB=pars(5);
    if k<1 | thetaA<=0 |thetaB<=0
        error('NA: switch rate function is not Holling type III');
    end
    ga=@(x)(pA*x^(k-1))/(1+thetaA*x^k);
    gb=@(x)(pB*x^(k-1))/(1+thetaB*x^k);
    
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
       elseif numel(ind)==3
           na=a(ind(2));
           fl=1;
       elseif numel(ind)>3
           error('NA: more than three equilibria found');
       elseif num<1e6
           num=2*num;
           a=linspace(1e-5,1-1e-5,num);
           b=1-a;
       else
          error('NA: could not obtain three co-existence equilibria');
       end           
           
    end