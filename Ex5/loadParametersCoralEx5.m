global massratio; global velocities;
E = 13.6e6; 
E = 14e6; 
nu = 0.3; rho = 792;  % See bottom for exp measures
%cu = 0.00024; % From decay test and for 100 elements 
l = 0.14; d = 0.005; I = pi * d^4 / 64 ;
ms = rho * pi * (d/2)^2;

numElements = 50; Nnodes = numElements+1;
%if (numElements==100); cu = 0.00024; elseif (numElements==50); cu = 0.00006;end; 
cu=0;
%velocities = 0.06:0.02:1.4; % Ur from 10 to 146, 48 points
%velocities = 0.01:0.01:1.4;
%velocities = 1.4;
velocities = [0.01:0.01:0.22 0.24:0.02:1.4]; % better dt and try to reach 0.22
vwindMax = velocities(end);
NR = length(velocities);
%finalTime_Dynamic = 2; finalTime_Static = NR; dt = 0.0005;%dt = 0.005 ; %dt = 0.01 ;
%finalTime_Dynamic = 0.5; finalTime_Static = NR; dt = 0.0005;
%finalTime_Dynamic = 2; finalTime_Static = NR; dt = 0.0002;%dt = 0.02;
finalTime_Dynamic = 2; finalTime_Static = NR; dt = 0.001;%dt = 0.02;
%finalTime_Dynamic = 5; finalTime_Static = NR; dt = 0.001;
rhoFluid = 1e3; muFluid = 1e-3; nuFluid = muFluid/rhoFluid; 
massratio = rho/rhoFluid ;
ma = rhoFluid * pi * (d/2)^2; 
St = 0.2; %epsilon=0.3; A = 12; 

nameLiftFunc = 'liftCoefWOM'; % 0.3
%nameDragFunc = 'dragCoefWOM'; % 1.2
nameDragFunc = 'innerDragCoefCircular'; % = f(Re)
cL0 = feval(nameLiftFunc, 0, 1); cD0 = feval(nameDragFunc, 0, 1); 

materials.hyperElasParams = [ E nu ] ;
materials.density         = rho      ;
% Corot element for large displacements
materials.hyperElasModel  = '1DrotEngStrain' ; 

elements_static(1).elemType = 'node' ;
elements_static(2).elemType = 'frame' ;
% cross section params
elements_static(2).elemCrossSecParams{1,1} = 'circle' ;
elements_static(2).elemCrossSecParams{2,1} = d          ;

numGaussPoints  = 4 ;
elements_static(2).aeroCoefs   = {nameDragFunc; []; [] }   ;
%  chord vector and gauss points
elements_static(2).elemTypeAero = [0 0 -d numGaussPoints true] ; % [chordVec1 chordVec2 chordVec3 numGauss aeroTangBool ]
% mass element formulation
elements_static(2).massMatType = 'consistent' ; 

elements_dynamic = elements_static;
elements_dynamic(2).aeroCoefs   = {nameDragFunc; nameLiftFunc; [] }   ;
elements_dynamic(2).elemTypeAero = [0 0 -d numGaussPoints false] ;

boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;

initialConds_static = struct() ;
initialConds_dynamic = struct() ;
%
% mesh 
CoralMeshEx5
% mesh.nodesCoords = [  zeros(numElements+1,1) zeros(numElements+1,1) (0:(numElements))'*l/numElements ] ;
% mesh.conecCell = { } ;
% mesh.conecCell{ 1, 1 } = [ 0 1 1 1 ] ;
% for i=1:numElements,
%   mesh.conecCell{ i+1,1 } = [ 1 2 0 i i+1 ] ;
% end

analysisSettings_stat.fluidProps = {rhoFluid; nuFluid; 'windVelLinear'} ;
%If drag reconfiguration then analysisSettings.geometricNonLinearAero  = true!! also it if is false then the lift direction will be constant
analysisSettings_stat.geometricNonLinearAero = true;
% time parameters
analysisSettings_stat.finalTime   = NR ;
analysisSettings_stat.deltaT      = 1  ; 
% numerical parameters
analysisSettings_stat.methodName     =   'newtonRaphson' ;  %'alphaHHT'; %
analysisSettings_stat.stopTolIts     =   50        ; 
analysisSettings_stat.stopTolDeltau  =   1e-10      ; 
analysisSettings_stat.stopTolForces  =   1e-8      ; 
% load settings
analysisSettings_stat.booleanSelfWeight = false ;

%Dynamic
analysisSettings_dyn = analysisSettings_stat;
analysisSettings_dyn.fluidProps = {rhoFluid; nuFluid; 'windVelUniform5'} ; % 'windUniformFaded'
analysisSettings_dyn.deltaT      = dt        ; 
analysisSettings_dyn.finalTime   = finalTime_Dynamic ;
analysisSettings_dyn.methodName     =   'newmark' ;  %'alphaHHT'; %
analysisSettings_dyn.stopTolIts     =   50        ; 
analysisSettings_dyn.stopTolDeltau  =   1e-8      ; 
analysisSettings_dyn.stopTolForces  =   1e-5      ; 

otherParams_stat.problemName = 'staticReconfiguration';
otherParams_dyn.problemName = 'dynamicVIV';
% % geometry in meters (coral)
% %
% 7 cm in heigh
%l = 0.07 ; d = 0.002; I = pi * d^4 / 64 ;
% % geometry in meters (coral)
% %
% Printed geometry
%l = 0.15 ; d = 0.005; I = pi * d^4 / 64 ;
% geometry Vortex-induced vibrations of cylinders bent by the flow Tristan Leclercq *, Emmanuel de Langre
%
%Re = d*va/nuFluid;
%Cy = rhoFluid*(va^2)*d*(l^3)/(2*E*I);
%CyLeclerc = rhoFluid*cD0*(va^2)*d*(l^3)/(2*E*I);
Gamma = l/d ;
ms = rho * pi * (d/2)^2; % Structural mass
% From Leclerc DeLangre 2018
ma = rhoFluid * pi * (d/2)^2; % Fluid added mass
va = 0.282;
UrLeclercq = ((l^2)*va/d)*sqrt((ms+ma)/(E*I))% proportional to gamma
f1Exp = 4.05; %Hz fft measures legh
UrExp = va/(f1Exp*d);
% stiffness
ku = 3*E*I/(l^3); 
% nodal lumped mass %m  = rho*pi*(d^2)*l/(4*2) ;
%
% van der pol equation parameters
%
epsilon = 0.3;
B = 2 * pi * St * va / d;
cq = epsilon * B; kq = B^2; 
A = 12; C2 = A/d; 
%--------
%
% ODE option
options = odeset('RelTol',1e-6);
%
%--------
%
% Plot parameters
%
spanPlot = 20 ; lw = 2.0 ; plotfontsize = 20 ;
%%
% TPE young modulus computation Should be Etpe = 6MPa according to website
dtpe = 0.005; ltpe = 0.15; 
%4 branches coral 
m4BranchesCoraldry = 0.00629; % Dry weight
m4BranchesCoralwet = 0.00730; % Wet weight
V4BranchesCoral = 9213e-9; % volume according to Fusion360
rhotpedry = m4BranchesCoraldry/V4BranchesCoral;
rhotpewet = m4BranchesCoralwet/V4BranchesCoral;
% rod without branch
Itpe = pi * dtpe^4 / 64 ;
Ttpe = 0.25; %10 periods of vibrations in 2,5 sec
beta1 = 1.87510; beta2=4.69409; beta1test = 1.87510/ltpe; beta2test=4.69409/ltpe;
Etpe = (2*pi/beta2^2*Ttpe)^2*(rhotpewet*2*pi*(dtpe/2)^2/Itpe); % with vibrations % omegai = betai^2 sqrt(EI/rho)
EItpe = 419e-6; % with 3 Points bending test (E = 13.6 Mpa)
Etpe = EItpe/Itpe; %  1.3657e+07 Pa
l = 0.14; Etpe = 23e6; % fits fMode1 = 4.05 Hz
% For cantilever beams clamped-free end only
betavect = [1.87510 4.69409 7.857476 10.9955 14.1372]; % beta*l for cantilever clamped-free end beam!
%freqStruct = (betavect.^2)*sqrt(Etpe*Itpe/(ms+ma)/(2*pi*l^2)) 
freqStruct = (betavect.^2)/(2*pi*l^2)*sqrt(Etpe*Itpe/(ms+ma)) ;