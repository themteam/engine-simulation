function [ Xu,Xb ] = fct_composition(x,y,phi)

%Calcul la compositon du mélange gaz frais (GF), gas brulés (GB) pour un 
%carburant CxHy à une richesse phi.
%
%[Xu,Xb]=fct_composition(x,y,phi)
%
%Xu : vecteur fraction molaire des des gaz frais [-]
%Xb : vecteur fraction molaire des des gaz brûlés [-]
%Selon la structure Xi(i=u ou i=b)=[H H2 O O2 OH H2O CO CO2 N N2 NO CxHy]
%
%x   : Nombre d'élément Carbone dans le carburant CxHy [-]
%y   : Nombre d'élément Hydrogène dans le carburant CxHy [-]
%phi : Richesse [-]

    %Proportion O2/N2 dans l'air
    psi = 3.76;

    %Mélange de gaz frais
    ntgf=phi+(x+y/4)*(1+psi);       %nombre de moles de gas frais
    Xu=[0 0 0 ((x+y/4)/ntgf) 0 0 0 0 0 (psi*(x+y/4)/ntgf) 0 (phi/ntgf)]';

    %Mélange de gaz brûlés
    ntgb=x*phi+y/2*phi+(1-phi)*(x+y/4)+(x+y/4)*psi;     %nombre de moles de gaz brulés
    Xb=[0 0 0 (1-phi)*(x+y/4)/(ntgb) 0 (y/2*phi)/(ntgb) 0 (x*phi)/(ntgb)...
        0 psi*(x+y/4)/(ntgb) 0 0]';

end

