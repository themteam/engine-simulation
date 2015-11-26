function out = systemeFunction1(theta,y,flag) 

%Fonction exprimant le système My=F. A resoudre avec ODE45. Fait appel à
%MatF1, fonction représentant F (vecteur du second membre) et MatM1, 
%fonction représentatnt M (matrice masse).
%
%theta : angle du vilebrequin représentant le temps
%y     : vecteur des grandeurs physiques inconnues dans la chambre de combustion
%y(1)  : P pression
%y(2)  : T température
%Y(3)  : m masse totale des gaz
%y(4)  : mu masse totale des gaz frais
%y(5)  : mb masse totale des gaz brûlés
%y(6)  : f=mb/m fraction de gaz brûlés
%y(7)  : mcapa masse de gaz brûlés dans la capacité remplie par backflow

    switch flag
        % calcul de la matrice excitation F
        case '' 
            out=MatF1(theta,y);

        %calcul de la matrice masse
        case 'mass' 
            out=MatM1(theta,y);
    end
end

