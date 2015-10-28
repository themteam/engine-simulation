function M=MatMass(theta,y)
global Cv r

    [V]=fct_volume(theta);
    M=zeros(3,3);
    M(1,1)=1;
    M(2,2)=y(1)*Cv;
    M(3,2)=-y(1)*r;
    M(3,3)=V;

end 