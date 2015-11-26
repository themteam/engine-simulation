function F = matF(theta,y)
    %Defines the right side F of the system My'=F
    
    [~,dv]=fct_volume(theta); 
    F=zeros(3,1);
    F(2)=-y(3)*dv;
    F(3)=-y(3)*dv;
end
