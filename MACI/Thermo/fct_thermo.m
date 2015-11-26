function out = fct_thermo(X,T,type,f)

    global R Mi Alow Ahigh ufu0 ufb0
    TableJanaf

    if (T<=1000 & T>=250)
        A=Alow;
    else (T<=3000 & T>1000);
        A=Ahigh;
    end

    switch type 
        case 'r'
            %Masse Molaire du melange
            M=Mi*X;

            %Constante des gaz parfaits
            out=R/M;

        case 'cp'
            vectT=[1 T T^2 T^3 T^4 T^5]';
            %Masse Molaire du melange
            M=Mi*X;
            %Chaleur spécifique à pression constante
            Cp=R*A(:,1:5)*vectT(1:5);
            out=(Cp'*X)/M;

        case 'cv'
            vectT=[1 T T^2 T^3 T^4 T^5]';
            %Masse Molaire du melange
            M=Mi*X;

            %Constante des gaz parfaits!!!
            r=R/M;

            %Chaleur spécifique à pression constante
            Cp=R*A(:,1:5)*vectT(1:5);
            cp=(Cp'*X)/M;

            %Chaleur spécifique à volume constant
            out=cp-r;

        case 'gamma'        
            vectT=[1 T T^2 T^3 T^4 T^5]';
            %Masse Molaire du melange
            M=Mi*X;
            %Constante des gaz parfaits!!!
            r=R/M; 
            %Chaleur spécifique à pression constante
            Cp=R*A(:,1:5)*vectT(1:5);
            cp=(Cp'*X)/M;
            %Chaleur spécifique à volume constant
            cv=cp-r;

            %Coefficient de Laplace!!!
            out=cp/cv;

        case 'h'
            vectT=[1 T T^2 T^3 T^4 T^5]';
            vectFraction=[1 1/2 1/3 1/4 1/5]';
            %Masse Molaire du melange
            M=Mi*X;

            %Enthalpie
            Hi=R*(A(:,1:5)*(vectFraction.*vectT(2:6))+A(:,6));
            out=(Hi'*X)/M;


         case 's'       
            vectT=[1 T T^2 T^3 T^4 T^5]';
            vectFraction=[1 1/2 1/3 1/4 1/5]';
            %Masse Molaire du melange
            M=Mi*X;
            %Entropie
            Si=R*(A(:,1)*log(vectT(2))+A(:,2:5)*(vectFraction(1:4).*vectT(1:4))+A(:,7));
            out=(Si'*X)/M;         

        case 'u'
            % Calcul de l'intégral 
            if (T>=250 && T<=1000) 
                T0 = 273.15;
                int_Cv = fct_prim_Cv(X,Alow,T) - fct_prim_Cv(X,Alow,T0);
            else
                Tlim = 1000; T0 = 273.15;
                int_Cv = fct_prim_Cv(X,Ahigh,T) - fct_prim_Cv(X,Ahigh,Tlim) + ...
                    fct_prim_Cv(X,Alow,Tlim) - fct_prim_Cv(X,Alow,T0);
            end
            out = (1-f)*ufu0 + f*ufb0 + int_Cv ;

    end

end