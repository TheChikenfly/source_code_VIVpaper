% Static ONSAS
close all, clear all ;
%global uniformUdot;
% add ONSAS path
accDir = pwd ;
onsas_path = getenv('ONSAS_PATH')
addpath( genpath( onsas_path ) );
% add shared files across examples 
addpath( genpath("./../shared/") );
%
loadParametersEx3
nameFuncVel = 'windVelUniform3';
numElements = 1;
materials.hyperElasParams = [ E nu ] ;
materials.density         = rho      ;
materials.hyperElasModel  = 'linearElastic' ;
%
% elements
%
% node
elements_static(1).elemType = 'node' ;
%
elements_static(2).elemType = 'frame' ;
% cross section params
elements_static(2).elemCrossSecParams{1,1} = 'circle' ;
elements_static(2).elemCrossSecParams{2,1} = d          ;
% hydro cross-section props
numGaussPoints  = 4 ;
%  only a distributed drag
elements_static(2).aeroCoefs   = {nameDragFunc; []; [] }   ;
elements_static(2).elemTypeAero = [0 0 -d numGaussPoints false] ; % [chordVec1 chordVec2 chordVec3 numGauss  ]
%
elements_static(2).massMatType = 'lumped' ; %'consistent' ;%
%
% boundaryConds
%md The first and unique BC corresponds to a welded condition for a cantilever beam
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%mdThe coordinates of the mesh nodes are given by the matrix:
mesh.nodesCoords = [  zeros(numElements+1,1) zeros(numElements+1,1) (0:(numElements))'*l/numElements] ;
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of nodes that compose the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first welded node is defined with material (M) zero since nodes don't have material, the first element (E) type (the first entry of the `elements` struct), and (B) is the first entry of the the `boundaryConds` struct. For (I) no non-homogeneous initial condition is considered (then zero is used) and finally the node is assigned:
mesh.conecCell{ 1, 1 } = [ 0 1 1  1 ] ;
%md Next the frame elements MEBI parameters are set. The frame material is the first material of `materials` struct, then $1$ is assigned. The second entry of the `elements` struct correspond to the frame element employed, so $2$ is set. Finally no BC and no IC is required for this element, then $0$ is used.  Consecutive nodes build the element so then the `mesh.conecCell` is:
for i=1:numElements,
  mesh.conecCell{ i+1,1 } = [ 1 2 0  i i+1 ] ;
end
initialConds_static = struct() ; % 0
% time parameters
analysisSettings_stat.fluidProps = {rhoFluid; nuFluid; nameFuncVel} ;
analysisSettings_stat.geometricNonLinearAero = false;
analysisSettings_stat.finalTime   = 1 ;
analysisSettings_stat.deltaT      = 1  ; 
analysisSettings_stat.methodName     =   'newtonRaphson' ;  %'alphaHHT'; %
analysisSettings_stat.stopTolIts     =   50        ; 
analysisSettings_stat.stopTolDeltau  =   1e-10      ; 
analysisSettings_stat.stopTolForces  =   1e-8      ; 
analysisSettings_stat.booleanSelfWeight = false ;
otherParams_stat.problemName = 'staticReconfiguration';
%
[matUsStat] = ONSAS( materials, elements_static, boundaryConds, initialConds_static, mesh, analysisSettings_stat, otherParams_stat ) ; 
%