function F=derivepid3(t,X)

global alpha a b c Betamax tauvar M

F=zeros(16,1);


if t<=tauvar             %Beta varie linéairement à partir de tauvar de 1 à Betamax à partir de t>tauvar
    Beta=1;
elseif t<=tauvar+50
    Beta=((Betamax-1)/50) * (t-tauvar) +1;
else
    Beta=Betamax;
end


for k=1 : 4
    F(4*(k-1) + 1)= -a * alpha * Beta * X(4*(k-1)+1) * X(4*(k-1)+2);
    F(4*(k-1)+2) = a*alpha*Beta*X(4*(k-1)+1) * X(4*(k-1)+2) - (c+b) * X(4*(k-1)+2);
    F(4*(k-1)+3) = b * X(4*(k-1)+2);
    F(4*(k-1)+4) = c * X(4*(k-1)+2);
end


T=zeros(4,4);
for i=1 : 4
    for j=1 : 4
        T(i,j) = X(4*(i-1) + j);
    end
end

R = M*T;

S=zeros(1,4);

for j = 1 : 4
    s=0;
    for i = 1 : 4
        s=s+M(i,j);
    end
    S(j)=s;
end

for k = 1 : 4
    for p = 1 : 4
        F(4*(k-1)+p) = F(4*(k-1)+p)+R(k,p)-S(k)*X(4*(k-1)+p);
    end
end
