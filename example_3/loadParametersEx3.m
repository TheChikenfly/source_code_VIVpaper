% This file contains the parameters fo the hydro_frame_force example
% For Example 3: Experimental investigation of vortex-induced vibration of long marine risers A.D.Trim & al
%
% Time parameters
%
E = 36.2e9; massratio = 1.6; rho = massratio*1e3; nu = .3; 
d = 0.027; l = d*1400 ; % 37.8 m
dint = 0.021; % Internal radius
I = pi *((d/2)^4 - (dint/2)^4) / 4 ; %I = pi * d^4 / 64 ;
cu = 0   ; % Should be 3% in Trim but not in Gao&al
T0 = 5e3 ; % Tension (Newton) 
pretension_strain = T0/(E*pi*d^2/4) ; 
%pretension_strain = T0/(E*pi*(d^2 - dint^2)/4) ;
numElements = nelem; 
dt = 0.005 ; finalTime = 7.14;%5; To have 20 periods of motion

rhoFluid = 1e3; muFluid = 1e-3; nuFluid = muFluid/rhoFluid; 
vwindMax = 0.4; 

nameLiftFunc = 'liftCoefWOM'; % 0.3
nameDragFunc = 'dragCoefWOM'; % 1.2
%nameDragFunc = 'dragCoefAmplified'; % 2
cL0 = feval(nameLiftFunc, 0, 1); cD0 = feval(nameDragFunc, 0, 1); 
St = 0.2;

materials.hyperElasParams = [ E nu ] ;
materials.density         = rho      ;
% Corot element for large displacements
materials.hyperElasModel  = '1DrotEngStrain' ; 

elements(1).elemType = 'node' ;
elements(2).elemType = 'frame' ;
elements(2).elemCrossSecParams{1,1} = 'pipe' ;
elements(2).elemCrossSecParams{2,1} = [d dint];

numGaussPoints  = 4 ;
elements(2).aeroCoefs   = {nameDragFunc; nameLiftFunc; [] }   ;
aeroTangBool = false;
elements(2).elemTypeAero = [0 0 -d numGaussPoints aeroTangBool] ; % [chordVec1 chordVec2 chordVec3 numGauss aeroTangBool ]
elements(2).massMatType = 'consistent' ; 

boundaryConds(1).imposDispDofs = [ 1 3 5 6] ;  % pinned + rotation along z blocked
boundaryConds(1).imposDispVals = [ 0 0 0 0] ;

initialConds = struct() ;
% mesh 
mesh.nodesCoords = [  zeros(numElements+1,1) zeros(numElements+1,1) (0:(numElements))'*l/numElements ] ;
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 1 ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1 numElements+1 ] ; 
for i=1:numElements,
  mesh.conecCell{ i+2,1 } = [ 1 2 0 i i+1 ] ;
end

analysisSettings.geometricNonLinearAero = true;
analysisSettings.fluidProps = {rhoFluid; nuFluid; 'windVelUniform'} ; %along x
analysisSettings.deltaT      = dt        ; 
analysisSettings.finalTime   = finalTime ;
analysisSettings.methodName     =   'newmark' ;  %'alphaHHT'; %
analysisSettings.stopTolIts     =   50        ; 
analysisSettings.stopTolDeltau  =   1e-10      ; 
analysisSettings.stopTolForces  =   1e-5      ; 

otherParams.problemName = 'ONSAS_Ex3';

qvect =  zeros(numElements*2,round(finalTime/dt)+1);
qvect(1:2:end,1) = (2*rand(numElements, 1)-1); %(2*rand(numElements, 1)-1)*0.001 ;
pvect =  zeros(numElements*2,round(finalTime/dt)+1);
pvect(1:2:end,1) = (2*rand(numElements, 1)-1); %(2*rand(numElements, 1)-1)*0.001 ;

%Static params
elements_static(1).elemType = 'node' ;
elements_static(2).elemType = 'frame' ;
elements_static(2).elemCrossSecParams{1,1} = 'pipe' ;
elements_static(2).elemCrossSecParams{2,1} = [d dint];
elements_static(2).aeroCoefs   = {nameDragFunc; []; [] }   ;
elements_static(2).elemTypeAero = [0 0 -d numGaussPoints false] ; % [chordVec1 chordVec2 chordVec3 numGauss aeroTangBool ]
elements_static(2).massMatType = 'consistent' ; 
initialConds_static = struct() ;
analysisSettings_stat.fluidProps = {rhoFluid; nuFluid; 'windVelUniform'} ;
%If drag reconfiguration then analysisSettings.geometricNonLinearAero  = true!! also it if is false then the lift direction will be constant
analysisSettings_stat.geometricNonLinearAero = true;
% time parameters
analysisSettings_stat.finalTime   = 1 ;
analysisSettings_stat.deltaT      = 1  ; 
analysisSettings_stat.methodName     =   'newtonRaphson' ;  %'alphaHHT'; %
analysisSettings_stat.stopTolIts     =   50        ; 
analysisSettings_stat.stopTolDeltau  =   1e-10      ; 
analysisSettings_stat.stopTolForces  =   1e-8      ; 
analysisSettings_stat.booleanSelfWeight = false ;
otherParams_stat.problemName = 'staticReconfiguration';

numNodes = numElements+1;
midNode = floor(numNodes/2); 

fw0 = St* vwindMax /d ;%natural shedding frequency
Re = d* vwindMax *rhoFluid / muFluid;

% UrLeclercq = (St*(l^2)* vwindMax /d)*sqrt((ms+mf)/(E*I)); % Not relevent
% mu = (ms+mf)/(rhoFluid*d^2);
%
% van der pol equation parameters
%
epsilon = 0.3; A = 12;
%--------
% Plot parameters
%
spanPlot = 20 ; lw = 2.0 ; plotfontsize = 20 ; %ms = 11 ; 
%--------
