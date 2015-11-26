clear all;

TableJanaf;


global R 
global Mi 
global Alow Ahigh 
global Xu Xb 
global ufu0 ufb0
global f


%Données globales
R = 8.314;
T=293.15; %en K
p=1e5; %en Pa
Mi= [1 2 16 32 17 18 28 44 14 28 30 114]*1e-3;
f = 0.95;

x = 8; y = 18;
Xu=[0 0 0 ((x+y/4)/(1+(x+y/4)*(1+3.76))) 0 0 0 0 0 (3.76*(x+y/4)/(1+(x+y/4)*(1+3.76))) 0 (1/(1+(x+y/4)*(1+3.76)))]';
Xb=[0 0 0 0 0 (y/4)/(1+(x+y/4)*(1+3.76)) 0 (x)/(1+(x+y/4)*(1+3.76)) 0 (3.76*(x+y/4)/(1+(x+y/4)*(1+3.76))) 0 0]';


[ufu0, ufb0] = fct_uformation();

%Exemple Mélange d'air
Xi= [ 0 0 0 0.21 0 0 0 0 0 0.79 0 0]';

% %Exemple Mélange de gaz frais à la stoechiométrique 
% x=8; y=18; 

% %Exemple Mélange de gaz brulés
% x=8; y=18; 


 out =fct_thermo(Xi,T,'u')
 % prim_Cv = fct_prim_Cv(Xi,Alow,T)
 