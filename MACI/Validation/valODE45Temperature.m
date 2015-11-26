function [] = valODE45Temperature( theta,T,V,To )

    %Validation  par TV^(gamma-1)=constante
    %Theta : angle villebrequin [degrees]
    %T     : température en K issue de la fonction ode45 [K]
    %V     : volume
    %To    : conditions initiales en température [K]

    gamma=1.4;      %air
    res=linspace(0,0,length(T));
    for i=1:length(T)
        res(i)=T(i)*V(i)^(gamma-1);
    end
    ref=To*V(1)^(gamma-1);
    figure;
    plot(theta,res,'o',theta,ref,'r*-');
    legend('TV ^{\gamma-1} au cours du temps','T_0V_0^{\gamma-1} aux conditions initales');
    xlabel('Angle du vilequin \theta (degrés)');
    ylabel('TV ^{\gamma-1}');
end

