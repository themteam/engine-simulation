clear all;
clc;
addpath(fullfile(pwd,'Volume'));
addpath(fullfile(pwd,'Thermo'));
addpath(fullfile(pwd,'Debit'));
addpath(fullfile(pwd,'Validation'));
addpath(fullfile(pwd,'Systeme'));

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
%Xu             : vecteur fraction molaire des gaz frais [-]
%Xb             : vecteur fraction molaire des gaz brûlés [-]
%Mi             : vecteur masse molaire des éléments composants les gaz
%                 [H H2 O O2 OH H2O CO CO2 N N2 NO CxHy] [kg/mol]
%Mu             : masse d'une mole de gaz frais [kg]
%Mb             : mase d'une mole de gaz brûlé [kg]
%
%R              : constante universelle des gaz parfaits [J/(mol.K)]
%ufu0           : energie interne massique de formation des gaz frais[J/kg]
%ufb0           : energie interne massique de formation des gaz brûlés[J/kg]
%
%N              : régime moteur critique [rpm]
%Dadm, Dech     : diamètre de la soupape d'admission et d'échappement [m]
%Tadm, padm     : température et pression à l'admission [K] [Pa]
%Tech, pech     : température et pression à l'échappement [K] [Pa]

%NbSoupAdm      : nombre de soupape d'admission [-]
%NbSoupEch      : nombre de soupape d'échappement [-]
%LmaxEch        : levée max de la soupape d'échappement [m]
%LmaxAdm        : levée max de la soupape d'admission [m]
%RFA            : retard à la fermeture de la soupape d'admission défini à 
%                 partir de 0 [°]
%AOA            : avance à l'ouverture de la soupape d'échappement définie 
%                 à partir de 0 [°]
%RFE            : retard à la fermeture de la soupape d'admission défini 
%                 à partir de 0 [°]
%AOE            : retard à l'ouverture de la soupape d'admission défini à 
%                 partir de 0 [°]
%cycle          : valeur d'angle représentant le cycle complet du moteur[°]
%--------------------------------------------------------------------------
global Vm lambda Cu EpsCompression rmanivel lbielle 
global Xu Xb Mi Mu Mb

global R ufu0 ufb0
global N Dadm Tadm padm Dech Tech pech

global AOA RFA AOE RFE LmaxAdm LmaxEch NbSoupAdm NbSoupEch cycle

lambda=3.15;
EpsCompression=9.4;
Cu=954*10^-6;
Vm=Cu/(EpsCompression-1);

[Xu, Xb]= fct_composition(8,18,1);
Mi= [1 2 16 32 17 18 28 44 14 28 30 114]*1e-3;
Mu = Mi *Xu;
Mb = Mi *Xb;

R=8.314;
[ufu0, ufb0] = fct_uformation();

N = 3000;
Dadm = 0.0368;
Dech = 0.029;
Tadm=300;
padm=1e5;
Tech=400;    
pech=1e5;

NbSoupAdm=1;
NbSoupEch=1;
LmaxEch=9.196e-3; 
LmaxAdm=9.196e-3;
RFA=262;
AOA=683; %700
RFE=38;%40
AOE=459;
cycle=720;
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

% Precision threshold (seuil)
ths = 1e-7;

%Conditions initiales
cond_init.P0 = 3e5;                           %Pa
cond_init.T0 = 400;                           %K
cond_init.r = fct_thermo(Xb, cond_init.T0, 'r');
cond_init.m0 = cond_init.P0 * fct_volume(0) / (r * cond_init.T0); %kg
cond_init.f0 = 0.95;
cond_init.mu0 = (1 - cond_init.f0) * m0;                %kg
cond_init.mb0 = cond_init.f0 * cond_init.m0;                      %kg
cond_init.mcapa0 = 0;                         %kg

% Détermination de la fonction P0.
[theta, M] = fct_P0(cond_init, ths);
P0 = M(:, 1);

%theta_ = ???
%P0_ = interp1(theta, P0, theta_)

options=odeset('Mass','M(t,y)','RelTol',1e-3,'AbsTol',[1e2,1e-1,1e-7,1e-7,1e-7,1e-2,1e-7]);
[theta,y]=ode45('systemeFunction1',0:0.5:720,y0,options);

tita=0:0.5:720;
leve_adm=linspace(0,0,length(tita));
leve_ech=linspace(0,0,length(tita));
for i=1:length(tita)
    leve_adm(i)=fLevee(tita(i),'adm');
    leve_ech(i)=fLevee(tita(i),'ech');
end
figure(1);
% subplot(3,1,1);
% plot(tita,leve_adm,tita,leve_ech);
% subplot(3,1,2);
plot(theta,y(:,1));
ylabel('p');
% subplot(3,1,3);
% plot(theta,y(:,2));
% ylabel('T');

figure(2);
% subplot(6,1,1);
% plot(tita,leve_adm,tita,leve_ech);
% subplot(6,1,2);
% plot(theta,y(:,3));
% ylabel('m');
% subplot(6,1,3);
% plot(theta,y(:,4));
% ylabel('mu');
% subplot(6,1,4);
% plot(theta,y(:,5));
% ylabel('mb');
% subplot(6,1,5);
plot(theta,y(:,6));
ylabel('f');
% subplot(6,1,6);
% plot(theta,y(:,7));
% ylabel('mcapa');



