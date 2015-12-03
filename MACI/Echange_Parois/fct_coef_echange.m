function h_c = fct_coef_echange(y, theta)
%ign = angle de début de combustion
%tps_comb = durée de combustion
%alésage = 70 mm
%course = 62 mm
%   Detailed explanation goes here
    global Cu d_alesage N course ign tps_comb P0struct T0struct
    
    Cm = course*N/30;
    C0 = 120;
    V1 = fct_volume(180);
    P1 = interp1(P0struct.theta,P0struct.pression,180);
    T1 = interp1(T0struct.theta,T0struct.temperature,180);
    P0actu=interp1(P0struct.theta,P0struct.pression,theta);
    
    if (fLevee(theta,'adm')==0 && fLevee(theta,'ech')==0)
        C1=2.28;
    else
        C1=6.18;
    end

    if (theta>=ign && theta<=(ign+tps_comb))
        if  y(1)-P0actu>0
            C2=6.22e-3;
        else
            disp('P-P0_<0 error');
            return
        end
    else
        C2=0;
    end

    %Conversion en bar
    P1=P1;
    P=y(1);
    P0actu=P0actu;
    disp(theta);
    h_c=C0*(d_alesage^(-0.2)*P^0.8*(C1*Cm+C2*Cu*T1/(P1*V1)*(P-P0actu))^0.8*y(2)^(-0.53))
end

