% Essai sur l'équation y'=2ty  sur [0,1], dont on connait bien la solution
global alpha a b c Betamax tauvar d h hr hc nlits M taufermeture Mfermeture;
%%%%%%%%%%%%%%%%%MODELISATION 4 CASES SEULEMENT
% coef mortalité
a=0.06; % coef infection sains
b=6e-3; % coef guerison
c=2e-2;

Betamax=1;  %Betamax est la valeur finale du coefficient Beta (cf dans les fonctions)
            % par défaut il est égal à 1 pour qu'il n'impacte pas le
            % programme
tauvar=200;

yinit1 = [0.995,0.005,0,0]; 

LTspan = [0,365]; % a modifier sous forme de matrice si on veut rajouter des etapes de confinement, 
                   %sous la forme [0,30;30,60;60,100;100,150;150,365] si
                   %l'on veut faire comme dans l'énoncé.
                   %[0,365] est le cas le plus simple ou alpha ne change pas
Alpha = [1];   %Alpha est la liste des alphas qui diffèrent selon les intervalles 
                %de temps de LTspan. Il doit être de la forme
                %[1,0.1,1,0.3,0.7]
                %[1] est donc le cas par défaut où alpha ne change pas et
                %où sa valeur n'affecte pas le coefficient a

hEuler = 1;   % Pas pour les méthodes d'Euler
paraODE = hEuler;
options = odeset('maxstep', paraODE); % Options pour ode45, pour avoir un pas de même ordre (max) 

[n,m]=size(LTspan); 




x1=[];
y1=[];

for i =1 : n; %on parcourt toute la liste des intervalles
    alpha=Alpha(i);
    
    %xp et yp (provisoires) sont les solutions que nous fournit le solveur pour un
    %intervalle de temps choisit et un yinitial donné
    [xp,yp] = Solveur(@derivepid1, LTspan(i,:), yinit1, hEuler, 2); 
    
    %les deux lignes concatènes les deux listes précédentes avec les listes
    %finales.
    x1 = [x1;xp];
    y1 = [y1;yp];
    
    %les trois lignes visent à récupérer les valeurs y1 à la fin pour
    %qu'elle puisse être utilisées en temps que condition initiale dans le
    %cas où la simulation n'est pas terminée et qu'on a pas fini de
    %parcourir la liste des intervalles de temps
    [k,l]=size(y1);
    yinit1 = y1(k,:);
    yinit1 = yinit1';
end


%affichage d'une figure contenant l'évolution des 4 cas en fonction du
%temps
figure()
plot(x1,y1(:,1),"r",x1,y1(:,2),"g",x1,y1(:,3),"b",x1,y1(:,4),"c",'LineWidth',3)
legend('sains','inféctés','rétablis','décédés')
grid on
xlabel('Temps (jours)')
ylabel('Populations %')
title(["Evolution des 4 cas en fonction du temps","Betamax = "+num2str(Betamax),"Tauvar = "+num2str(tauvar),"Alpha = "+num2str(Alpha),"LTspan = "+num2str(LTspan)])


%%%%%%%%%%%%%%%%%%%%%%MODELISATION 5 CASES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% coef mortalité
alpha = 1;
a=1.e-1; % coef infection sains
b=6e-3; % coef guerison
c=2e-2;

Betamax=1;
tauvar=200;

%Les 4 variables suivantes sont nouvelles car elles prennent en compte un
%nouveau type de la population.
h=0.05; % pourcentage des cas graves à hospitaliser
hr=0.1;
hc=0.1;
nlits=0.05;

%yinit doit aussi prendre en compte qu'il n'y a pas d'hospitalisés au début
yinit2 = [0.995,0.005,0,0,0];


LTspan = [0,365]; % a modifier sous forme de matrice si on veut rajouter des etapes de confinement
Alpha = [1];

hEuler = 1;   % Pas pour les méthodes d'Euler
paraODE = hEuler;
options = odeset('maxstep', paraODE); % Options pour ode45, pour avoir un pas de même ordre (max) 


[n,m]=size(LTspan);




x2=[];
y2=[];

for i =1 : n;
    alpha=Alpha(i);
    
    [xp,yp] = Solveur(@derivepid2, LTspan(i,:), yinit2, hEuler, 2);
    x2 = [x2;xp];
    y2 = [y2;yp];
    
    [k,l]=size(y2);
    yinit2 = y2(k,:);
    yinit2 = yinit2';
end


%On affiche comme précédement une figure contenant 5 courbes en focntion du
%temps. Une pour chaque cas.
figure()
plot(x2,y2(:,1),"r",x2,y2(:,2),"g",x2,y2(:,3),"b",x2,y2(:,4),"c",x2,y2(:,5),'LineWidth',3)
legend('sains','inféctés','rétablis','décédés','hospitalisés')
grid on
xlabel('Temps (jours)')
ylabel('Populations %')
title(["Evolution des 5 cas en fonction du temps","Betamax = "+num2str(Betamax),"Tauvar = "+num2str(tauvar),"Alpha = "+num2str(Alpha),"LTspan = "+num2str(LTspan)])


%%%%%%%%%%%%%%%%%AVEC LES VILLES


% coef mortalité
alpha = 1;
a=1.e-1; % coef infection sains
b=6e-3; % coef guerison
c=2e-2;
Betamax=1;
tauvar=200;
%On se retouve dans le cas où l'on a que 4 types de population

%M est la patrice de déplacements entre les deux pays, ici il y a 4 pays.
%On a M(i,j) le pourcentage de population voyageant du pays j au pays i par
%pas de temps. On considère que tous les types de population voyagent de
%manière equivalente.
M = [0,1,1,1;1,0,1,1;0,1,0,1;0,0,1,0]/100;

yinit3 = [0.995,0.005,0,0];

%la matrice y finale est construite de manière à ce que la hauteur soit le
%temps et que la longueur soit définie comme suit : on renseigne chaqye
%type de population à la suite (chaque colonne) et on met côte à côte
%chaque pays. D'où le fait que les coditions initiales soient une
%répétition. On peut par exemple changer cela en indiquant que tout le
%monde est sain dans chaque pays et que la proportion initiale de malade 
%est positive dans un seul pays 
yinit3 = [yinit3,yinit3,yinit3,yinit3];

LTspan = [0,365]; % a modifier sous forme de matrice si on veut rajouter des etapes de confinement
Alpha = [1];

hEuler = 1;   % Pas pour les méthodes d'Euler
paraODE = hEuler;
options = odeset('maxstep', paraODE); % Options pour ode45, pour avoir un pas de même ordre (max) 
[n,m]=size(LTspan);




x3=[];
y3=[];

for i =1 : n;
    alpha=Alpha(i);
    
    [xp,yp] = Solveur(@derivepid3, LTspan(i,:), yinit3, hEuler, 2);
    x3 = [x3;xp];
    y3 = [y3;yp];
    
    [k,l]=size(y3);
    yinit3 = y3(k,:);
    yinit3 = yinit3';
end

%la matrice ynhabtot renseigne le pourcentage total de la population en
%fonction du temps (cela permet de bien voir les déplacements de
%population)
ynhabtot=zeros(size(y3,1),4);
ynhabtot(:,1)=y3(:,1)+y3(:,2)+y3(:,3)+y3(:,4);
ynhabtot(:,2)=y3(:,5)+y3(:,6)+y3(:,7)+y3(:,8);
ynhabtot(:,3)=y3(:,9)+y3(:,10)+y3(:,11)+y3(:,12);
ynhabtot(:,4)=y3(:,13)+y3(:,14)+y3(:,15)+y3(:,16);

%On affiche ici 4 figures différentes (une pour chaque pays) avec 4 courbes
%correnspondant aux 4 types de population.
figure()
plot(x3,y3(:,1),"r",x3,y3(:,2),"g",x3,y3(:,3),"b",x3,y3(:,4),"c",x3,ynhabtot(:,1),"k",'LineWidth',3)
legend('sains','inféctés','rétablis','décédés','total')
grid on
xlabel('Temps (jours)')
ylabel('Populations %')
title(["Pays n° 1 Evolution des 4 cas en fonction du temps","Betamax = "+num2str(Betamax),"Tauvar = "+num2str(tauvar),"Alpha = "+num2str(Alpha),"LTspan = "+num2str(LTspan)])
figure()
plot(x3,y3(:,5),"r",x3,y3(:,6),"g",x3,y3(:,7),"b",x3,y3(:,8),"c",x3,ynhabtot(:,2),"k",'LineWidth',3)
legend('sains','inféctés','rétablis','décédés','total')
grid on
xlabel('Temps (jours)')
ylabel('Populations %')
title(["Pays n° 2 Evolution des 4 cas en fonction du temps","Betamax = "+num2str(Betamax),"Tauvar = "+num2str(tauvar),"Alpha = "+num2str(Alpha),"LTspan = "+num2str(LTspan)])
figure()
plot(x3,y3(:,9),"r",x3,y3(:,10),"g",x3,y3(:,11),"b",x3,y3(:,12),"c",x3,ynhabtot(:,3),"k",'LineWidth',3)
legend('sains','inféctés','rétablis','décédés','total')
grid on
xlabel('Temps (jours)')
ylabel('Populations %')
title(["Pays n° 3 Evolution des 4 cas en fonction du temps","Betamax = "+num2str(Betamax),"Tauvar = "+num2str(tauvar),"Alpha = "+num2str(Alpha),"LTspan = "+num2str(LTspan)])
figure()
plot(x3,y3(:,13),"r",x3,y3(:,14),"g",x3,y3(:,15),"b",x3,y3(:,16),"c",x3,ynhabtot(:,4),"k",'LineWidth',3)
legend('sains','inféctés','rétablis','décédés','total')
grid on
xlabel('Temps (jours)')
ylabel('Populations %')
title(["Pays n° 4 Evolution des 4 cas en fonction du temps","Betamax = "+num2str(Betamax),"Tauvar = "+num2str(tauvar),"Alpha = "+num2str(Alpha),"LTspan = "+num2str(LTspan)])



%%%%%%%%%%%%%%%%%AVEC LES VILLES ET 5 CASES


% coef mortalité
alpha = 1;
a=1.e-1; % coef infection sains
b=6e-3; % coef guerison
c=2e-2;

Betamax=1;
tauvar=200;

M = [0,1,1,1;1,0,1,1;0,1,0,1;0,0,1,0]/100;


taufermeture = 100; %temps à partir duquel les déplacemnts de population changent

%La matrice Mfermeture est la nouvelle matrice de déplacements, ici le pays
%4 a totalement fermé ses frontières. M change dans la fonctions
%@derivepid4
Mfermeture = [0,1,1,0;1,0,1,0;0,1,0,0;0,0,0,0]/100;

d=0.; % coef resucceptibilité
h=0.05; % pourcentage des cas graves à hospitaliser
hr=0.1;
hc=0.1;
nlits=0.05;

yinit4 = [0.995,0.005,0,0,0];
yinit4 = [yinit4,yinit4,yinit4,yinit4];

LTspan = [0,365]; % a modifier sous forme de matrice si on veut rajouter des etapes de confinement
Alpha = [1];

hEuler = 1;   % Pas pour les méthodes d'Euler
paraODE = hEuler;
options = odeset('maxstep', paraODE); % Options pour ode45, pour avoir un pas de même ordre (max) 

[n,m]=size(LTspan);

x4=[];
y4=[];

for i =1 : n;
    alpha=Alpha(i);
    
    [xp,yp] = Solveur(@derivepid4, LTspan(i,:), yinit4, hEuler, 2);
    x4 = [x4;xp];
    y4 = [y4;yp];
    
    [k,l]=size(y4);
    yinit4 = y4(k,:);
    yinit4 = yinit4';
end

%Meme fonctionnement que précédement sauf qu'il faut faire attention au
%format de y4 qui contient un type de population par ville supplémentaire
size(y4)
ynhabtot=zeros(size(y4,1),4);
ynhabtot(:,1)=y4(:,1)+y4(:,2)+y4(:,3)+y4(:,4)+y4(:,5);
ynhabtot(:,2)=y4(:,6)+y4(:,7)+y4(:,8)+y4(:,9)+y4(:,10);
ynhabtot(:,3)=y4(:,11)+y4(:,12)+y4(:,13)+y4(:,14)+y4(:,15);
ynhabtot(:,4)=y4(:,16)+y4(:,17)+y4(:,18)+y4(:,19)+y4(:,20);

%Affichage des 5 courbes pour les 4 pays.
figure()
plot(x4,y4(:,1),"r",x4,y4(:,2),"g",x4,y4(:,3),"b",x4,y4(:,4),"c",x4,y4(:,5),"m",x4,ynhabtot(:,1),"k",'LineWidth',3)
legend('sains','inféctés','rétablis','décédés','hospitalisés')
grid on
xlabel('Temps (jours)')
ylabel('Populations %')
title(["Pays n° 1 Evolution des 5 cas en fonction du temps","Betamax = "+num2str(Betamax),"Tauvar = "+num2str(tauvar),"Alpha = "+num2str(Alpha),"LTspan = "+num2str(LTspan),"M = ",num2str(M(1,:)),num2str(M(2,:)),num2str(M(3,:)),num2str(M(4,:)), "Mfermeture = ",num2str(Mfermeture(1,:)),num2str(Mfermeture(2,:)),num2str(Mfermeture(3,:)),num2str(Mfermeture(4,:))])
figure()
plot(x4,y4(:,6),"r",x4,y4(:,7),"g",x4,y4(:,8),"b",x4,y4(:,9),"c",x4,y4(:,10),"m",x4,ynhabtot(:,2),"k",'LineWidth',3)
legend('sains','inféctés','rétablis','décédés','hospitalisés')
grid on
xlabel('Temps (jours)')
ylabel('Populations %')
title(["Pays n° 2 Evolution des 5 cas en fonction du temps","Betamax = "+num2str(Betamax),"Tauvar = "+num2str(tauvar),"Alpha = "+num2str(Alpha),"LTspan = "+num2str(LTspan)])
figure()
plot(x4,y4(:,11),"r",x4,y4(:,12),"g",x4,y4(:,13),"b",x4,y4(:,14),"c",x4,y4(:,15),"m",x4,ynhabtot(:,3),"k",'LineWidth',3)
legend('sains','inféctés','rétablis','décédés','hospitalisés')
grid on
xlabel('Temps (jours)')
ylabel('Populations %')
title(["Pays n° 3 Evolution des 5 cas en fonction du temps","Betamax = "+num2str(Betamax),"Tauvar = "+num2str(tauvar),"Alpha = "+num2str(Alpha),"LTspan = "+num2str(LTspan)])
figure()
plot(x4,y4(:,16),"r",x4,y4(:,17),"g",x4,y4(:,18),"b",x4,y4(:,19),"c",x4,y4(:,20),"m",x4,ynhabtot(:,4),"k",'LineWidth',3)
legend('sains','inféctés','rétablis','décédés','hospitalisés')
grid on
xlabel('Temps (jours)')
ylabel('Populations %')
title(["Pays n° 4 Evolution des 5 cas en fonction du temps","Betamax = "+num2str(Betamax),"Tauvar = "+num2str(tauvar),"Alpha = "+num2str(Alpha),"LTspan = "+num2str(LTspan)])