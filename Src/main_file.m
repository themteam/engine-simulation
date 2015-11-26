clc;
clear all;
%--------------------------------------------------------------------------
%                           CONSTANTES
%--------------------------------------------------------------------------
global Vm                   %Volume au point mort bas (PMB) [m^3]
global lambda               %Rapport longueur de bielle/rayon manivelle [-]
global Cu                   %Cylindrée [m^3]
global EpsCompression       %Rapport de compression Vpmb/Vpmh [-]
global rmanivel             %Rayon de la manivelle [m]
global lbielle              %Longueur de bielle [m]
global Xu                   %fraction de gaz frais
global Xb                   %fraction de gaz brûlés 
global Mi                   %vecteur de masse molaire des composants
global R 
global Mu
global Mb
global N 
global Dadm
global ufu0
global ufb0
lambda=3.15;
EpsCompression=9.4;
Cu=954*10^-6;
Vm=Cu/(EpsCompression-1);
M_air=29e-3;
R=8.314;
teta=0:0.01:360;               %Vecteur angle vilebrequin
[Xu, Xb]= fct_composition(8,18,1);
Mi= [1 2 16 32 17 18 28 44 14 28 30 114]*1e-3;
Mu = Mi *Xu;
Mb = Mi *Xb;
N = 3000;  %tr/min
Dadm = 0.0368;
[ufu0, ufb0] = fct_uformation();
%--------------------------------------------------------------------------

V=linspace(0,0,length(teta));
dvdta=linspace(0,0,length(teta));
for i=1:length(teta)
    [V(i),dvdta(i)]=fct_volume(teta(i));
end

% %Verification de la dérivée
% for i=1:length(teta)-1
%     res(i)=V(i+1)-dvdta(i)*(teta(i+1)-teta(i));
% end
% subplot(3,1,1);
% plot(teta,V);
% grid on
% subplot(3,1,2)
% plot(teta(1:end-1),res);
% grid on
% subplot(3,1,3);
% plot(teta,dvdta);
% grid on

% %Conditions initiales
% P0=1e5; %Pa
% T0=300; %K
% m0=P0*fct_volume(180)/(r*T0); %kg
% y0=[m0,T0,P0]; 
% 
% options=odeset('Mass','M(t,y)','RelTol',1e-3,'AbsTol',[1e-7, 1e-1, 1e2]);
% [theta,y]=ode45('systemeFunction',180:0.01:540,y0,options);


% figure;plot(theta,y(:,3))
% 
% valODE45_Pression(theta,y(:,3),V,y0(3))
% 
% dP=y(:,3);
% P=zeros(length(dP),1);
% P(1)=P0;
% 
% for i=2:length(theta)
%     P(i)=P(i-1)+dP(i-1)*(theta(i)-theta(i-1));
% end
% plot(P)




%Conditions initiales
P0=1e5;                       %Pa
T0=300;                       %K
r=fct_thermo(Xb,T0,'r');
m0=P0*fct_volume(0)/(r*T0);   %kg
f0=0.1;
mu0=(1-f0)*m0;                %kg
mb0=f0*m0;                    %kg
mcapa0=0;                     %kg
y0=[P0,T0,m0,mu0,mb0,f0,mcapa0]; 

options=odeset('Mass','M(t,y)','RelTol',1e-3,'AbsTol',[1e2,1e-1,1e-7,1e-7,1e-7,1e-2,1e-7]);
[theta,y]=ode45('systemeFunction1',0:0.01:360,y0,options);
