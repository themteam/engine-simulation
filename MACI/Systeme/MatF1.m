function F = MatF1(theta,y)

%Fonction représentant le second membre F du système à résoudre My=F.
%Elle est appelée uniquement par la fonction systemeFunction1 et ne peut 
%être utilisée de façon indépendante.
%
%F=MatF1(theta,y)
%
%theta : angle du vilebrequin représentant le temps
%y     : vecteur des inconnus [p,T,m,mu,mb,f,mcapa]
%F     : vecteur du second membre du système My=F

%-------------------------------------------------------------------------
    
%-------------------------------------------------------------------------
    global Xu Xb Mu Mb Tadm Tech P0struct T0struct
    
    %X : vecteur fractions molaires du mélange
    X = ((y(4)/Mu)*Xu+(y(5)/Mb)*Xb)./(y(4)/(Mu)+y(5)/(Mb));
    F=zeros(7,1);
    if isempty(P0struct)||isempty(T0struct)
        disp('dq_parois nul');
        dq_parois = 0; %hypothèse sans transfert de chaleur
    else
        h_c = fct_coef_echange(y, theta);
        dq_parois = fct_echange_chaleur( h_c, y, theta )
    end
    if y(7)<=0
        %Cas mcapa nulle
        h_adm=fct_thermo(Xu,Tadm,'h');
        h_cyl=fct_thermo(X,y(2),'h');
        h_ech=fct_thermo(Xb,Tech,'h'); 
        % hypothèse : uniquement des gaz brûlés à l'échappement ?
        [~,dv]=fct_volume(theta);
        [dm_adm, dm_adm_bf]=fct_debit(theta,y,'adm');
        [dm_ech,dm_ech_bf]=fct_debit(theta,y,'ech');
        F(1)=-y(1)*dv;
        F(2)=-y(1)*dv+dq_parois+dm_adm*h_adm+dm_adm_bf*h_cyl+...
              dm_ech*h_cyl+dm_ech_bf*h_ech; 
        F(3)=0;    
        F(4)=dm_adm+(1-y(6))*dm_adm_bf+(1-y(6))*dm_ech;       
        F(5)=y(6)*dm_ech+y(6)*dm_adm_bf+dm_ech_bf;
        F(6)=0; 
        F(7)=-y(6)*dm_adm_bf;
        
    else
        %Cas mcapa positif
        h_adm=fct_thermo(Xb,Tadm,'h'); % pourquoi T=Tadm ==> hyp brassage 
        % température immédiate ? 
        h_cyl=fct_thermo(X,y(2),'h');
        h_ech=fct_thermo(Xb,Tech,'h'); % hypothèse : uniquement des gaz
        %brûlés à l'échappement ? (à écrire dans la fonction débit pour les
        %conditions sur la pression d'échappement)
        [~,dv]=fct_volume(theta);
        [dm_adm, dm_adm_bf]=fct_debit(theta,y,'adm');
        [dm_ech,dm_ech_bf]=fct_debit(theta,y,'ech');
        
        F(1)=-y(1)*dv;
        F(2)=-y(1)*dv+dq_parois+dm_adm*h_adm+dm_adm_bf*h_cyl+dm_ech*h_cyl+dm_ech_bf*h_ech;
        F(3)=0; 
        F(4)=(1-y(6))*dm_adm_bf+(1-y(6))*dm_ech;       
        F(5)=dm_adm+y(6)*dm_ech+y(6)*dm_adm_bf+dm_ech_bf;
        F(6)=0; 
        F(7)=-dm_adm-y(6)*dm_adm_bf;
    end
end

