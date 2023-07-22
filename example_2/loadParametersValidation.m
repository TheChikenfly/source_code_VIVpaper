% This file contains the parameters fo the hydro_frame_force example
% For Example 5
accDir = pwd ;
%addpath( genpath( [ accDir '/../../../ONSAS.m/src'] ) );
addpath( genpath( [ accDir '/../shared'] ) );
addpath( genpath( [ accDir '/Zsolutions'] ) );
%
% Time 
%
global finalTime; global vwindMax; global AMBool ; global St; global massratio; global fluidFlowBool;
AMBool = false; fluidFlowBool = true; 

% First try with 'windStepbyStep' mesh analysis
E = 5e10; dt = 0.0028 ; finalTime = 50; rhoFluid = 1020; rho = rhoFluid; massratio = 1; numElements = 100; vwindMax = 0.04; fluidFlowBool = true; load('Zwinstep_va=0.04_Nelem=100_FT50_Nstep=20_coefacc=0.2.mat');
E = 5e10; dt = 0.0028 ; finalTime = 50; rhoFluid = 1020; rho = rhoFluid; massratio = 1; numElements = 50; vwindMax = 0.04; fluidFlowBool = true; load('Zwinstep_va=0.04_Nelem=50_FT50_Nstep=20_coefacc=0.2.mat');
E = 5e10; dt = 0.0028 ; finalTime = 50; rhoFluid = 1020; rho = rhoFluid; massratio = 1; numElements = 30; vwindMax = 0.04; fluidFlowBool = true; load('Zwinstep_va=0.04_Nelem=30_FT50_Nstep=20_coefacc=0.2.mat');
E = 5e10; dt = 0.0028 ; finalTime = 50; rhoFluid = 1020; rho = rhoFluid; massratio = 1; numElements = 15; vwindMax = 0.04; fluidFlowBool = true; load('Zwinstep_va=0.04_Nelem=15_FT50_Nstep=20_coefacc=0.2.mat');
% First try with fluidFlowBool = true;
%E = 5e10; dt = 0.0028 ; finalTime = 50; rhoFluid = 1020; rho = rhoFluid; massratio = 1; numElements = 30; vwindMax = 0.030; load('Zsolutions/ZLinear_va=0.03_Nelem=30_FT50.mat');%'windVelLinear'
% Classic to compare with:
%E = 5e10; dt = 0.0028 ; finalTime = 100; vwindMax = 0.2; rhoFluid = 1020; rho = rhoFluid; numElements = 20; epsilon=0.3; A = 12; load('ZsolV4Nelem=20_FT=100_E=5e+10_dt=0.0028_vmax=0.2_eps=0.3_A=12.mat');
%--------
% Fluid properties
%
%rhoFluid = 1020; nuFluid = 1e-6; 
nuFluid = 1e-6; 
% % absolute mean velocity of water at x=[0 0 0] and time 0 
nameFuncVel = 'windVelUniform'; 
%nameFuncVel = 'windUniformFaded'; 
%nameFuncVel = 'windVelLinear'; 
%nameFuncVel = 'windStepbyStep'; 
va_vector = feval(nameFuncVel, [0, 0, 0], 0);
va        = vwindMax; %
%va = 0.124
% extract lift and drag which are constant (betaRel and Re doesen´t matter)
nameLiftFunc = 'liftCoefWOM'; % 0.3
nameDragFunc = 'dragCoefAmplified'; % 2.0
cL0 = feval(nameLiftFunc, 0, 1); cD0 = feval(nameDragFunc, 0, 1); 
%--------   
%
% Material and geometric properties
% material
nu = .3; 
% % geometry in meters (coral)
% %
% l = 0.1 ; d = 0.002; I = pi * d^4 / 64 ;
% geometry Vortex-induced vibrations of cylinders bent by the flow Tristan Leclercq *, Emmanuel de Langre
%
l = 1 ; d = 0.001; I = pi * d^4 / 64 ;% d = 0.001; I = pi * d^4 / 64 ;
%pretension_strain = 0;
%
St = 0.2; 
fw0 = St*va/d ;%natural shedding frequency
Re = d*va/nuFluid;
Cy = rhoFluid*(va^2)*d*(l^3)/(2*E*I);
CyLeclercq = rhoFluid*cD0*(va^2)*d*(l^3)/(2*E*I);
Gamma = l/d ;
ms = rho * pi * (d/2)^2; % Structural mass
mf = rhoFluid * pi * (d/2)^2; % Fluid added mass
Omegas = sqrt(3*E*I/l^3) ;
Omegaf = 2*pi*St*va/d ; % 2 pi fw0
FA_Damping = 2.0/(4*pi*St) ; % gamma = 0.8 in Facchinetti et al 2004

% UrFaccinetti2 = 2*pi*va/(Omegas*d)
mu = (ms+mf)/(rhoFluid*d^2);
%M = 0.05/mu;
% From Leclercq DeLangre 2018
ma = mf; % added mass
beta = ma / (ma + ms); 
%va =0.0150 
% va=0.0005:0.001:0.0150;
UrLeclercq = (St*(l^2)*va/d)*sqrt((ms+ma)/(E*I))
% Ur linear in vpr: Ur = vpr*(6.46-0.58)/(0.04-0.0036)
%UrLeclercq = (St*(l^2)*va/d)*sqrt((ms)/(E*I)) % No AM 
f1 = sqrt((E*I)/(ms+ma))/(St*(l^2)); % First frequency when vibrating in the fluid! f1 = Structural Mode 2!!
% CyLeclercqback = (2*cD0/pi) * UrLeclercq^2 * beta/(Gamma*St^2); % = CyLeclercq
% UrLeclercqBack = sqrt((pi*Gamma*CyLeclercq*St^2)/(cD0* beta*2)); % =UrLeclercq
% Numerical parameters
% structure damping
cu = 0;
% van der pol equation parameters
%
epsilon = 0.3; A = 12;
B = 2 * pi * St * va / d;
cq = epsilon * B; kq = B^2; 
C2 = A/d; 
%--------
%
% Initial Conditions
q0 = 0.001; dq0 = 0; ddq0 = 0; 
X0 = 0; dX0 = 0; ddX0 = 0;
Y0 = 0; dY0 = 0; ddY0 = 0;
Z0 = 0; dZ0 = 0; ddZ0 = 0;
%--------
% Plot parameters
%
spanPlot = 20 ; lw = 2.0 ; plotfontsize = 20 ; %ms = 11 ; 

% ODE option
options = odeset('RelTol',1e-6);
%
%
% TPE young modulus computation Should be Etpe = 6MPa according to website
% dtpe = 0.005; ltpe = 0.15; mtpe = 0.002;%kg
% rhotpe = mtpe/(ltpe*2*pi*(dtpe/2)^2);
% Itpe = pi * dtpe^4 / 64 ;
% Ttpe = 0.25; %10 periods of vibrations in 2,5 sec
% beta1 = 1.87510; beta2=4.69409; beta1test = 1.87510/ltpe; beta2test=4.69409/ltpe;
% Etpe = (2*pi/beta2^2*Ttpe)^2*(rhotpe*2*pi*(dtpe/2)^2/Itpe); % weird that it doesn't depend on the length
% Etpetest = (2*pi/beta2test^2*Ttpe)^2*(rhotpe/Itpe); % omegai = betai^2 sqrt(EI/rho)