%--------------------------------------------------------------------------
%                           CONSTANTES
%--------------------------------------------------------------------------
global Vm                   %Volume au point mort bas (PMB) [m^3]
global lambda               %Rapport longueur de bielle/rayon manivelle [-]
global Cu                   %Cylindrée [m^3]
global EpsCompression       %Rapport de compression Vpmb/Vpmh [-]
global rmanivel             %Rayon de la manivelle [m]
global lbielle              %Longueur de bielle [m]


global Cv
global r 

lambda=3.15;
EpsCompression=9.4;
Cu=954*10^-6;
Vm=Cu/(EpsCompression-1);

Cv=718;
M_air=29e-3;
R=8.314;
r=R/M_air;
teta=180:1:720;               %Vecteur angle vilebrequin


%--------------------------------------------------------------------------

V=linspace(0,0,length(teta));
dvdta=linspace(0,0,length(teta));
for i=1:length(teta)
    [V(i),dvdta(i)]=fct_volume(toRadians(teta(i)));
end
%plot(teta,V,'*-');
%plot(teta,dvdta,'*-');

%conditions initiales
P0=1e5; %Pa
T0=300; %K
m0=P0*fct_volume(180)/(r*T0); %kg
y0=[m0,T0,P0]; 

%options=odeset('Mass','MatM','RelTol',1e-3,'AbsTol',[1e-7, 1e-1, 1e2]);

options=odeset('Mass','m(t,y)','RelTol',1e-3,'AbsTol',[1e-7, 1e-1, 1e2]);
[theta,y]=ode45('MatF',[180:0.1:540],y0,options);

figure;plot(theta,y(:,1))



