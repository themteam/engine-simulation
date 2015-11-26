function [] = valODE45Pression( theta,P,V,Po )

    %Validation  par PV^(gamma)=constante
    %Theta : angle villebrequin [degrees]
    %P     : pression en Pa issue de la fonction ode45 [Pa]
    %V     : volume [m^3]
    %Po    : conditions initiales en pression [pa]

    gamma=1.4;      %air
    res=linspace(0,0,length(P));
    for i=1:length(P)
        res(i)=P(i)*V(i)^(gamma);
    end
    ref=Po*V(1)^(gamma);
    figure;
    plot(theta,res,'o',theta,ref,'r*-');
    legend('PV ^{\gamma} au cours du temps','P_0V_0^{\gamma} aux conditions initales');
    xlabel('Angle du vilequin \theta (degrés)');
    ylabel('PV ^{\gamma}');
end


