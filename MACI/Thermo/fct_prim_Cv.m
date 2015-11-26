function prim_Cv = fct_prim_Cv(Xi,A,T)

    vectT=[1 T T^2 T^3 T^4 T^5]';
    vectFraction=[1 1/2 1/3 1/4 1/5]';
    r = fct_thermo(Xi,T,'r');
    prim_Cv = (r * ((A(:,1)-1)*T + A(:,2:5)* (vectFraction(2:5).*vectT(3:6))))'*Xi;

end