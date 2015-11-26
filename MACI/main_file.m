clear all;
clc;
addpath(fullfile(pwd,'Volume'));
addpath(fullfile(pwd,'Thermo'));
addpath(fullfile(pwd,'Debit'));
addpath(fullfile(pwd,'Validation'));

%--------------------------------------------------------------------------
%                           CONSTANTES
%--------------------------------------------------------------------------
%Vm             : volume au point mort bas (PMB) [m^3]
%lamba          : rapport longueur de bielle/rayon manivelle [-]
%Cu             : cylindrée [m^3]
%EpsCompression : rapport de compression Vpmb/Vpmh [-]
%rmanivel       : rayon de la manivelle [m]
%lbielle        : longueur de bielle [m]
%
%Xu             : vecteur fraction molaire des gaz frais
%Xb             : vecteur fraction molaire des gaz brûlés
%Mi             : vecteur masse molaire des éléments composants les gaz
%                 [H H2 O O2 OH H2O CO CO2 N N2 NO CxHy]
%Mu             : masse d'une mole de gaz frais
%Mb             : mase d'une mole de gaz brûlé
%
%R              : constante universelle des gaz parfaits
%ufu0           :
%ufb0           : 
%--------------------------------------------------------------------------
global Vm lambda Cu EpsCompression rmanivel lbielle 
global Xu Xb Mi Mu Mb

global R ufu0 ufb0
global Dadm Tadm padm Dech Tech pech

global N 

lambda=3.15;
EpsCompression=9.4;
Cu=954*10^-6;
Vm=Cu/(EpsCompression-1);
M_air=29e-3;
R=8.314;
[Xu, Xb]= fct_composition(8,18,1);
Mi= [1 2 16 32 17 18 28 44 14 28 30 114]*1e-3;
Mu = Mi *Xu;
Mb = Mi *Xb;
N = 3000;  %tr/min
Dadm = 0.0368;
Dech = 0.0368;
[ufu0, ufb0] = fct_uformation();
Tadm=300;
padm=1e5;
Tech=400;
pech=1e5;
%--------------------------------------------------------------------------

% V=linspace(0,0,length(teta));
% dvdta=linspace(0,0,length(teta));
% for i=1:length(teta)
%     [V(i),dvdta(i)]=fct_volume(teta(i));
% end

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
P0=3e5;                       %Pa
T0=400;                       %K
r=fct_thermo(Xb,T0,'r');
m0=P0*fct_volume(0)/(r*T0);   %kg
f0=0.1;
mu0=(1-f0)*m0;                %kg
mb0=f0*m0;                    %kg
mcapa0=0;                     %kg
y0=[P0,T0,m0,mu0,mb0,f0,mcapa0]; 


        %[dm_adm, dm_adm_bf]=fct_debit(0,y0,'adm');
        %[dm_ech,dm_ech_bf]=fct_debit(0,y0,'ech');


options=odeset('Mass','M(t,y)','RelTol',1e-3,'AbsTol',[1e2,1e-1,1e-7,1e-7,1e-7,1e-2,1e-7]);
[theta,y]=ode45('systemeFunction1',0:0.5:720,y0,options);


lanceLev
tita=0:0.5:720;
leve_adm=linspace(0,0,length(tita));
leve_ech=linspace(0,0,length(tita));
for i=1:length(tita)
    leve_adm(i)=fLevee(tita(i),'adm');
    leve_ech(i)=fLevee(tita(i),'ech');
end
figure(1);
subplot(3,1,1);
plot(tita,leve_adm,tita,leve_ech);
subplot(3,1,2);
plot(theta,y(:,1));
ylabel('p');
subplot(3,1,3);
plot(theta,y(:,2));
ylabel('T');

figure(2);
subplot(6,1,1);
plot(tita,leve_adm,tita,leve_ech);
subplot(6,1,2);
plot(theta,y(:,3));
ylabel('m');
subplot(6,1,3);
plot(theta,y(:,4));
ylabel('mu');
subplot(6,1,4);
plot(theta,y(:,5));
ylabel('mb');
subplot(6,1,5);
plot(theta,y(:,6));
ylabel('f');
subplot(6,1,6);
plot(theta,y(:,7));
ylabel('mcapa');
