function CD=fct_CD(lv)
%Fonction qui calcule le coefficient de débit en fonction de la levée de la
%soupape.
%
%CD=fct_CD(lv)
%
%CD : coefficient de débit [-]
%lv : levée de la soupape [m]

    CD = 0.05*10^3*lv+0.65;

end
