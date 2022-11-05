function F=derivepid1(t,X)

global alpha a b c Betamax tauvar


if t<=tauvar             %Beta varie linéairement à partir de tauvar de 1 à Betamax à partir de t>tauvar
    Beta=1;
elseif t<=tauvar+50
    Beta=((Betamax-1)/50) * (t-tauvar) +1;
else
    Beta=Betamax;
end


F=[-a*Beta * X(1) * alpha * X(2); 
    a*Beta * X(1) * alpha *X(2) - (c+b) * X(2);
    b * X(2); 
    c * X(2)];