function out = MatF(theta,y,flag)
global Cv r 
[V,dv]=fct_volume(theta); 
switch flag
    
       % calcul de la matrice masse
    case 'mass' 
        out=zeros(3,3);
        out(1,1)=1;
        out(2,2)=y(1)*Cv;
        out(3,2)=-y(1)*r;
        out(3,3)=V;
        
        %calcul de la fonction excitation F
    case '' 
        out=zeros(3,1);
        out(2)=-y(3)*dv;
        out(3)=-y(3)*dv;
end

end

