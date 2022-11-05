function F=derivepid2(t,X)

global alpha a b c Betamax tauvar h hr hc nlits;


if t<=tauvar             %Beta varie linéairement à partir de tauvar de 1 à Betamax à partir de t>tauvar
    Beta=1;
elseif t<=tauvar+50
    Beta=((Betamax-1)/50) * (t-tauvar) +1;
else
    Beta=Betamax;
end

F=[-a*Beta * X(1) * alpha * X(2); 
    a*Beta * X(1) * alpha *X(2) - (c+b) * X(2)-min(h * X(2),max(0,nlits-X(5))); 
    b * X(2) + hr * X(5); 
    c * X(2) + hc * X(5);
    -(hc+hr) * X(5) + min(h * X(2),max(0,nlits-X(5)))];