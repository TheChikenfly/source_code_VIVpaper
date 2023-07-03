global massratio; global velocities;

E = 5e10; nu = 0.3; rho = 1000; cu = 0;
l = 1; d = 0.001; 
I = pi * d^4 / 64 ;
ms = rho * pi * (d/2)^2;
massratio=1;

numElements = 100; 

velocities = 0.0005:0.001:0.3;%1.4; % Urmax = 4.8, 30 points
vwindMax = velocities(end);
NR = length(velocities);
dt = 0.0005 ; finalTime_Dynamic = 2; finalTime_Static = NR;
dt = 0.0005 ; finalTime_Dynamic = 5; 

rhoFluid = 1000; nuFluid = 1e-6; 
ma = rhoFluid * pi * (d/2)^2; 
epsilon=0.3; A = 12; St = 0.2;

betavect = [1.87510 4.69409 7.857476 10.9955 14.1372]; % beta*l for cantilever clamped-free end beam!
freqStruct = (betavect.^2)*sqrt(E*I/(ms+ma))/(2*pi*l^2)
%fs1= sqrt(E*I/(ms+ma))/l^2 %1.25 Hz
fs1= freqStruct(1) %0.7 Hz

nameLiftFunc = 'liftCoefWOM'; % 0.3
nameDragFunc = 'dragCoefAmplified'; % 2.0
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
mesh.nodesCoords = [  zeros(numElements+1,1) zeros(numElements+1,1) (0:(numElements))'*l/numElements ] ;
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 1 ] ;
for i=1:numElements,
  mesh.conecCell{ i+1,1 } = [ 1 2 0 i i+1 ] ;
end

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
analysisSettings_dyn.fluidProps = {rhoFluid; nuFluid; 'windVelUniform'} ; % 'windUniformFaded'
analysisSettings_dyn.deltaT      = dt        ; 
analysisSettings_dyn.finalTime   = finalTime_Dynamic ;
analysisSettings_dyn.methodName     =   'newmark' ;  %'alphaHHT'; %
analysisSettings_dyn.stopTolIts     =   50        ; 
analysisSettings_dyn.stopTolDeltau  =   1e-8      ; 
analysisSettings_dyn.stopTolForces  =   1e-5      ; 

otherParams_stat.problemName = 'staticReconfiguration';
otherParams_dyn.problemName = 'dynamicVIV';
 