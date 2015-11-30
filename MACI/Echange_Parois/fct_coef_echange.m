function [ h_c ] = fct_coef_echange(y, theta, ign, tps_comb)
%ign = angle de début de combustion
%tps_comb = durée de combustion
%alésage = 70 mm
%course = 62 mm
%   Detailed explanation goes here
global Cu D p1 t1 P0 N course
[V,dV]=fct_volume(theta)

Cm=course*N/30
C0=120;
V1=fct_volume(180);
if (theta <= 180 && theta > 0)
    p1=y(1);
    t1=y(2);
else 
    P1=p1;
    T1=t1;
end
if theta < 360 && V>Vpmh
     
else
end
        
if  y(1)-P0>0
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
    
    h_c=C0*(D^(-0.2)*y(1)^0.8*(C1*Cm+C2*Cu*T1/(P1*V1)*(y(1)-P0))^0.8*y(2)^(-0.53));
   
else 
    disp('P-P0<0 error'); 
    return
end

