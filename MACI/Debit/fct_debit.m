function [dm_ae,dm_ae_bf]=fct_debit(theta,y,type)
%Fonction qui calcule le débit d'air entrant ou sortant de la chambre de
%combustion par les soupapes.
%
%[dm_ae,dm_ae_bf]=fct_debit(theta,y,type)
%
%theta  : angle du vilebrequin représentant le temps
%y      : vecteur des inconnus [p,T,m,mu,mb,f,mcapa]
%type   : 'adm' pour le cas admission ou 'ech' pour le cas de l'échappement
%dm_ae  : débit de masse d'air entrant (resp. sortant) de la chambre dans
%         le cas admission (resp. echappement)
%dm_ae_bf: débit de masse d'air sortant (resp. entrant) de la chambre en 
%          backflow dans le cas admission(reps. sortant)

    global padm Tadm pech Tech Xb Xu Dadm Dech Mu Mb N

    switch type
        case 'adm'
            lv=fLevee(theta,'adm');
            S= pi*Dadm*lv ;
            CD=fct_CD(lv);
            sct=S*CD;
            if padm > y(1)
                flagadm = 1;
                pamont = padm;
                Tamont=Tadm;
                paval=y(1);
                if y(7)>0
                    X=Xb;
                else
                    X=Xu;
                end
            else
                flagadm=0;
                pamont=y(1);
                Tamont=y(2);
                paval=padm;
                X = ((y(4)/Mu)*Xu+(y(5)/Mb)*Xb)/(y(4)/(Mu)+y(5)/(Mb));
            end
            gamma_am=fct_thermo(X,Tamont,'gamma');
            r_am=fct_thermo(X,Tamont,'r');

            Xcrit=(2/(gamma_am+1))^(gamma_am/(gamma_am-1));
            if (paval/pamont > Xcrit)
                XP=paval/pamont;
            else
                XP=Xcrit;
            end

            dm = 1/(6*N)*sct*pamont*sqrt(2*gamma_am/((gamma_am-1)*r_am*Tamont)*(XP^(2/gamma_am)-XP^((gamma_am+1)/gamma_am)));

            if flagadm
                dm_ae= dm;
                dm_ae_bf = 0;
            else
                dm_ae= 0;
                dm_ae_bf = -dm;
            end

        case 'ech'
            lv=fLevee(theta,'ech');
            S= pi*Dech*lv ;
            CD=fct_CD(lv);
            sct=S*CD;
            if pech > y(1)
                flagadm = 1;
                pamont = pech;
                Tamont=Tech;
                paval=y(1);
                X=Xb;
            else
                flagadm=0;
                pamont=y(1);
                Tamont=y(2);
                paval=pech;
                X = ((y(4)/Mu)*Xu+(y(5)/Mb)*Xb)/(y(4)/(Mu)+y(5)/(Mb));
            end
            gamma_am=fct_thermo(X,Tamont,'gamma');
            r_am=fct_thermo(X,Tamont,'r');

            Xcrit=(2/(gamma_am+1))^(gamma_am/(gamma_am-1));
            if (paval/pamont > Xcrit)
                XP=paval/pamont;
            else
                XP=Xcrit;
            end

            dm = 1/(6*N)*sct*pamont*sqrt(2*gamma_am/((gamma_am-1)*r_am*Tamont)*(XP^(2/gamma_am)-XP^((gamma_am+1)/gamma_am)));

            if flagadm
                dm_ae = 0;
                dm_ae_bf = dm;
            else
                dm_ae = -dm;
                dm_ae_bf = 0;
            end
    end
end