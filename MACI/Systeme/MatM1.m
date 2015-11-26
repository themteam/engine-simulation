function M = MatM1( theta, y )

%Fonction représentant la matrice masse M du système à résoudre My=F.
%Elle est appelée uniquement par la fonction systemeFunction1 et ne peut 
%être utilisée de façon indépendante.
%
%M=matM1(theta,y)
%
%theta : angle du vilebrequin représentant le temps
%y     : vecteur des inconnus [p,T,m,mu,mb,f,mcapa]
%M     : matrice masse
   
    global R Mu Mb ufb0 ufu0 Xu Xb
    
    %X : vecteur fractions molaires du mélange
    X = ((y(4)/Mu)*Xu+(y(5)/Mb)*Xb)./(y(4)/(Mu)+y(5)/(Mb));
    Cv=fct_thermo(X,y(2),'cv');
    u=fct_thermo(X,y(2),'u',y(6));
    r=fct_thermo(X,y(2),'r');
    [V]=fct_volume(theta);
    
    M=zeros(7,7);
    
    M(1,:)=[V -y(3)*r -r*y(2) 0 0 -y(3)*y(2)*(R/Mb-R/Mu) 0];
    M(2,:)=[0 y(3)*Cv u 0 0 y(3)*(ufb0-ufu0) 0];
    M(3,:)=[0 0 1 -1 -1 0 0];
    M(4,:)=[0 0 0 1 0 0 0];
    M(5,:)=[0 0 0 0 1 0 0];
    M(6,:)=[0 0 y(6) 0 -1 y(3) 0];
    M(7,:)=[0 0 0 0 0 0 1];
end

