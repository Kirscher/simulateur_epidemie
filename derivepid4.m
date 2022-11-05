function F=derivepid4(t,X)

global alpha a b c Betamax tauvar M h hr hc nlits taufermeture Mfermeture

%20 puisqu'il y a 4 pays et 5 types de population.
F=zeros(20,1);

if t<=tauvar             %Beta varie linéairement à partir de tauvar de 1 à Betamax à partir de t>tauvar
    Beta=1;
elseif t<=tauvar+50
    Beta=((Betamax-1)/50) * (t-tauvar) +1;
else
    Beta=Betamax;
end

%C'est ici qu'on change la matrice M
if t>taufermeture
    M=Mfermeture;  % le pays 4 change ses frontières à partir de taufermeture
end


%Ici on applique les quations de transmissions de la maladie avant
%d'utiliser les déplacements de population. 
for k=1 : 4
    F(5*(k-1) + 1)= -a * alpha * Beta * X(5*(k-1)+1) * X(5*(k-1)+2);
    F(5*(k-1)+2) = a*alpha*Beta*X(5*(k-1)+1) * X(5*(k-1)+2) - (c+b) * X(5*(k-1)+2) - min(h * X(5*(k-1)+2),max(0,nlits-X(5*(k-1)+5)));
    F(5*(k-1)+3) = b * X(5*(k-1)+2) + hr * X(5*(k-1)+5);
    F(5*(k-1)+4) = c * X(5*(k-1)+2) + hc * X(5*(k-1)+5);
    F(5*(k-1)+5) = -(hc+hr) * X(5*(k-1)+5) + min(h * X(5*(k-1)+2),max(0,nlits-X(5*(k-1)+5)));
end

%T est une matrice qui va nous aider à faire plus rapidement un calcul. En
%effet les données au temps i sont renseignés comme une matrice colonne,
%mais il est plus simple d'avoir un format avec une ligen par pays et une
%colonne par type de population.
T=zeros(4,5);
for i=1 : 4
    for j=1 : 5
        T(i,j) = X(5*(i-1) + j);
    end
end

%R est ici la matrice de la même taille que T qui permet de voir le nombre
%de personnes arrivant dans un pays.
R = M*T;


%On veut calculer ici le pourcentage de population quittant les pays (il ne
%faut pas l'oublier pour éviter que la population totale des pays
%augmentent). S est une matrice ligne renseignant le pourcentage totale de
%population à quitter le pays.
S=zeros(1,5);


for j = 1 : 4
    s=0;
    for i = 1 : 4
        s=s+M(i,j); %on parcourt la matrice M par colonne.
    end
    S(j)=s;
end

%Ici on rajoute à la matrice colonne F existante les différentes arrivées
%et départs, doù on a cette equation.
for k = 1 : 4
    for p = 1 : 5
        F(5*(k-1)+p) = F(5*(k-1)+p)+R(k,p)-S(k)*X(5*(k-1)+p);
    end
end