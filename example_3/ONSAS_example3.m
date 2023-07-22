%
%md# WOM cantilever beam example
%
%close all, clear all ;
%tic ;
%
% declare global variables
%
global qvect; global pvect; %VIV boolean is called inside hydroFrameForces 
global VIVBool; %VIV boolean is called inside hydroFrameForces
global constantLiftDir; %constantLiftDir is called inside hydroFrameForces
global uniformUdot; %constantLiftDir is called inside hydroFrameForces
VIVBool = true;
if VIVBool
  constantLiftDir = true; uniformUdot = false; % constantLiftDir = true for validation !
end
% BE SURE THAT:
% Consistent not lumped, 
%
% materials
%
materials.hyperElasParams = [ E nu ] ;
materials.density         = rho      ;
% Corot element for large displacements
materials.hyperElasModel  = '1DrotEngStrain' ; 
%
% node
elements(1).elemType = 'node' ;
%
% hydro frame
%
elements(2).elemType = 'frame' ;
% cross section params
elements(2).elemCrossSecParams{1,1} = 'pipe' ;
elements(2).elemCrossSecParams{2} = [d dint] ;
% hydro cross-section props
numGaussPoints  = 4 ;
elements(2).aeroCoefs   = {nameDragFunc; nameLiftFunc; [] }   ;
%  chord vector and gauss points
elements(2).elemTypeAero = [0 0 -d numGaussPoints false] ; % [chordVec1 chordVec2 chordVec3 numGauss computeAeroBool]
% mass element formulation
elements(2).massMatType = 'consistent' ; 
%
% boundaryConds
%
%md The first and unique BC corresponds to a welded condition for a cantilever beam
boundaryConds(1).imposDispDofs = [ 1 3 4 5 ] ; % pinned + rotation along y blocked
boundaryConds(1).imposDispVals = [ 0 0 0 0 ] ;
boundaryConds(2).imposDispDofs = [ 1 3 5 ] ; % pinned
boundaryConds(2).imposDispVals = [ 0 0 0 ] ;
%
initialConds = struct() ;
%
%mdThe coordinates of the mesh nodes are given by the matrix:
mesh.nodesCoords = [  zeros(numElements+1,1) (0:(numElements))'*l/numElements zeros(numElements+1,1) ] ;
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of nodes that compose the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first and last welded nodes are defined :
mesh.conecCell{ 1, 1 } = [ 0 1 1 1 ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1 numElements+1 ] ; %pinned
%md Next the frame elements MEBI parameters are set
for i=1:numElements,
  mesh.conecCell{ i+2,1 } = [ 1 2 0  i i+1 ] ; % +2 to skip the 2 first BC
end
% %
% %md homogeneous initial conditions are considered, then an empty struct is set:
% initialConds = struct() ;
% %
% % mesh 
% meshgeometry
% myCell2Mat( mesh.conecCell ) % To visualise connection
% Initialize qvect pvect
qvect =  zeros(numElements*2,round(finalTime/dt)+1);
qvect(1:2:end,1) = 2*(2*rand(numElements, 1)-1); %(2*rand(numElements, 1)-1)*0.001 ;
pvect =  zeros(numElements*2,round(finalTime/dt)+1);
pvect(1:2:end,1) = (2*rand(numElements, 1)-1); %(2*rand(numElements, 1)-1)*0.001 ;
% analysisSettings
%
% fluid properties
analysisSettings.fluidProps = {rhoFluid; nuFluid; nameFuncVel} ;
%If drag reconfiguration then analysisSettings.geometricNonLinearAero  = true!! also it if is false then the lift direction will be constant
analysisSettings.geometricNonLinearAero = true;
% time parameters
analysisSettings.finalTime   = finalTime ;
analysisSettings.deltaT      = dt        ; 
% numerical parameters
analysisSettings.methodName     =   'newmark' ;
analysisSettings.stopTolIts     =   50        ;
analysisSettings.stopTolDeltau  =   1e-8      ; 
analysisSettings.stopTolForces  =   1e-5      ; 
% load settings
analysisSettings.booleanSelfWeight = false ;
% % Folder path
%
% analysisSettings
%
% otherParams
otherParams.nodalDispDamping = cu ;
% Do NOT write vtk !
otherParams.problemName      = strcat('ONSAS_example3') ;
%otherParams.plotsFormat      = 'vtk' ;
%
% Run ONSAS
%
[ matUsDynLD2D, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 
%matUsDynLD2D
% Extract numerical solution
%linear disp first node
%zdefNum1 = mesh.nodesCoords(:,3) + matUsDynLD2D(5:6:end, :) ;
%xdefNum1 = mesh.nodesCoords(:,1) + matUsDynLD2D(1:6:end, :) ;
%linear disp last node
zdefNumall = matUsDynLD2D(5:6:end, :) ;% Z for all nodes
xdefNumall = matUsDynLD2D(1:6:end, :) ;% X for all nodes
%times vector
times = linspace(0,analysisSettings.finalTime, size(matUsDynLD2D, 2));
% % xlswrite('valsNewmark.xlsx',ydefNum(2,:));
% % xlswrite('valsNewmark.xlsx',qvect(:, 1));
% DeltaZ = (max(zdefNumlast(1:end)) - min(zdefNumlast(1:end)))/d;
%DeltaX = (max(xdefNumlast(1:end)) - min(xdefNumlast(1:end)))/d;
%elapTime = toc ;
%fprintf( strcat('\n The elapsed execution time was ', num2str(elapTime), ' seconds', '\n' ) );
%-----------------------------------------------------------------------------------------------------
