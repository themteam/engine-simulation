function out = systemeFunction(theta,y,flag) 
    switch flag
        % calcul de la matrice masse
        case '' 
            out=matF(theta,y);

        %calcul de la fonction excitation F
        case 'mass' 
            out=matMass(theta,y);
    end
end