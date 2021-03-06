clear all;
clc;
addpath(fullfile(pwd,'Volume'));
addpath(fullfile(pwd,'Thermo'));
addpath(fullfile(pwd,'Debit'));
addpath(fullfile(pwd,'Validation'));
addpath(fullfile(pwd,'Systeme'));
addpath(fullfile(pwd,'Echange_Parois'));

%--------------------------------------------------------------------------
%                           CONSTANTES
%--------------------------------------------------------------------------
% Vm             : volume au point mort bas (PMB) [m^3]
% lamba          : rapport longueur de bielle/rayon manivelle [-]
% Cu             : cylindr�e [m^3]
% EpsCompression : rapport de compression Vpmb/Vpmh [-]
% rmanivel       : rayon de la manivelle [m]
% lbielle        : longueur de bielle [m]
% d_alesage      : al�sage [m]
% d_piston       : diam�tre du piston [m]
% course         : course du piston [m]
%
% Xu             : vecteur fraction molaire des gaz frais [-]
% Xb             : vecteur fraction molaire des gaz br�l�s [-]
% Mi             : vecteur masse molaire des �l�ments composants les gaz
%                  [H H2 O O2 OH H2O CO CO2 N N2 NO CxHy] [kg/mol]
% Mu             : masse d'une mole de gaz frais [kg]
% Mb             : mase d'une mole de gaz br�l� [kg]
%
% R              : constante universelle des gaz parfaits [J/(mol.K)]
% ufu0           : energie interne massique de formation des gaz frais[J/kg]
% ufb0           : energie interne massique de formation des gaz br�l�s[J/kg]
%
% N              : r�gime moteur critique [rpm]
% Dadm, Dech     : diam�tre de la soupape d'admission et d'�chappement [m]
% Tadm, padm     : temp�rature et pression � l'admission [K] [Pa]
% Tech, pech     : temp�rature et pression � l'�chappement [K] [Pa]
%
% NbSoupAdm      : nombre de soupape d'admission [-]
% NbSoupEch      : nombre de soupape d'�chappement [-]
% LmaxEch        : lev�e max de la soupape d'�chappement [m]
% LmaxAdm        : lev�e max de la soupape d'admission [m]
% RFA            : retard � la fermeture de la soupape d'admission d�fini � 
%                  partir de 0 [�]
% AOA            : avance � l'ouverture de la soupape d'�chappement d�finie 
%                  � partir de 0 [�]
% RFE            : retard � la fermeture de la soupape d'admission d�fini 
%                  � partir de 0 [�]
% AOE            : retard � l'ouverture de la soupape d'admission d�fini � 
%                  partir de 0 [�]
% cycle          : valeur d'angle repr�sentant le cycle complet du moteur[�]
% ign            : angle d'allumage
% tps_comb       : 
% ths            : Precision threshold (seuil) pour la convergence du
%                  syst�me sans combustion (d�termination de P0_)
%--------------------------------------------------------------------------
global Vm lambda Cu EpsCompression rmanivel lbielle 
global d_alesage d_piston course
 
global Xu Xb Mi Mu Mb

global R ufu0 ufb0
global N Dadm Tadm padm Dech Tech pech

global AOA RFA AOE RFE LmaxAdm LmaxEch NbSoupAdm NbSoupEch cycle
global ign tps_comb P0struct T0struct
global ths

lbielle = 126.8e-3;
rmanivel = 38.5e-3;
lambda = lbielle / rmanivel;
EpsCompression=9.4;
Cu=954*10^-6;
Vm=Cu/(EpsCompression-1);
d_piston = 69.95e-3;
d_alesage = 70e-3;
course = 62e-3;

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

ign = 340;
tps_comb = 40;

ths = 1e-5;

%--------------------------------------------------------------------------

% V=linspace(0,0,length(teta));
% dvdta=linspace(0,0,length(teta));
% for i=1:length(teta)
%     [V(i),dvdta(i)]=fct_volume(teta(i));
% end

%Conditions initiales
cond_init.P0 = 3e5;                           %Pa
cond_init.T0 = 500;                           %K
cond_init.r = fct_thermo(Xb, cond_init.T0, 'r');
cond_init.m0 = cond_init.P0 * fct_volume(0) / (cond_init.r * cond_init.T0); %kg
cond_init.f0 = 0.90;
cond_init.mu0 = (1 - cond_init.f0) * cond_init.m0;                %kg
cond_init.mb0 = cond_init.f0 * cond_init.m0;                      %kg
cond_init.mcapa0 = 1e-6;                         %kg

% D�termination de la fonction P0.
[P0struct, T0struct] = fct_P0(cond_init, ths);

 
%%

y0 = [cond_init.P0,  cond_init.T0, cond_init.m0, ...
      cond_init.mu0, cond_init.mb0, cond_init.f0, ...
      cond_init.mcapa0];
  
options=odeset('Mass','M(t,y)','RelTol',1e-3,'AbsTol',[1e2,1e-1,1e-7,1e-7,1e-7,1e-2,1e-7]);
[theta,y]=ode45('systemeFunction1',0:1:380,y0,options);

% tita=0:0.5:720;
% leve_adm=linspace(0,0,length(tita));
% leve_ech=linspace(0,0,length(tita));
% for i=1:length(tita)
%     leve_adm(i)=fLevee(tita(i),'adm');
%     leve_ech(i)=fLevee(tita(i),'ech');
% end

figure(2);
% subplot(3,1,1);
% plot(tita,leve_adm,tita,leve_ech);
subplot(3,1,2);
plot(theta,y(:,1));
ylabel('p');
subplot(3,1,3);
plot(theta,y(:,2));
ylabel('T');

figure(3);
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



