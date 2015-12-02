function [ h_c ] = fct_coef_echange(y, theta, ign, tps_comb)
%ign = angle de début de combustion
%tps_comb = durée de combustion
%alésage = 70 mm
%course = 62 mm
%   Detailed explanation goes here
    global Cu D_alesage N course ths Xb
    
    %Conditions initiales
    cond_init.P0 = 3e5;                           %Pa
    cond_init.T0 = 400;                           %K
    cond_init.r = fct_thermo(Xb, cond_init.T0, 'r');
    cond_init.m0 = cond_init.P0 * fct_volume(0) / (cond_init.r * cond_init.T0); %kg
    cond_init.f0 = 0.95;
    cond_init.mu0 = (1 - cond_init.f0) * cond_init.m0;                %kg
    cond_init.mb0 = cond_init.f0 * cond_init.m0;                      %kg
    cond_init.mcapa0 = 0;                         %kg
    
    % Détermination de la fonction P0.
    [theta, M] = fct_P0(cond_init, ths);
    P0 = M(:, 1);
    T0 = M(:, 2);

    theta_ = 0:0.01:720;
    P0_ = interp1(theta, P0, theta_);
    T0_ = interp1(theta, T0, theta_);
    
    Cm = course*N/30;
    C0 = 120;
    V1 = fct_volume(180);
    P1 = P0_(180);
    T1 = T0_(180);
        
    if  y(1)-P0_>0
        if (fLevee(theta,'adm')==0 && fLevee(theta,'ech')==0)
            C1=2.28;
        else
            C1=6.18;
        end

        if (theta>=ign && theta<=(ign+tps_comb))
            C2=6.22e-3;
        else
            C2=0;
        end
    
        h_c=C0*(D_alesage^(-0.2)*y(1)^0.8*(C1*Cm+C2*Cu*T1/(P1*V1)*(y(1)-P0_))^0.8*y(2)^(-0.53));
   
    else 
        disp('P-P0_<0 error'); 
        return
    end
end

