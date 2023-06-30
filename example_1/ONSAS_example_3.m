%
%md# WOM cantilever beam example
%
%close all, clear all ;
%
% declare global variables
%
global pvect; global qvect; %VIV boolean is called inside hydroFrameForces 
global ILVIVBool; global VIVBool; %VIV boolean is called inside hydroFrameForces
global constantLiftDir; %constantLiftDir is called inside hydroFrameForces
global uniformUdot;
VIVBool = true;
if VIVBool
  constantLiftDir = false; % True ??
  uniformUdot = true ;
end
%
% add ONSAS path
%
% accDir = pwd ;
% onsas_path = getenv('ONSAS_PATH')
% addpath( genpath( onsas_path ) );
% % add shared files across examples 
% addpath( genpath("./../shared/") );
% %
% load parameters script
%
loadParametersEx3
% nameFuncVel = 'windVelUniform3';
% ONSAS MEBI parameters
% % materials
% %
% materials.hyperElasParams = [ E nu ] ;
% materials.density         = rho      ;
% % Corot element for large displacements
% %materials.hyperElasModel  = '1DrotEngStrain' ; % Both work
% materials.hyperElasModel  = 'linearElastic' ;
% elements
% node
elements(1).elemType = 'node' ;
elements(2).elemType = 'frame' ;
% cross section params
elements(2).elemCrossSecParams{1,1} = 'circle' ;
elements(2).elemCrossSecParams{2,1} = d          ;
% hydro cross-section props
numGaussPoints  = 4 ;
%  only a distributed lift
%elements(2).aeroCoefs   = {[]; nameLiftFunc; [] }   ;
%  only a distributed lift and drag
elements(2).aeroCoefs   = {nameDragFunc; nameLiftFunc; [] }   ;
%  chord vector and gauss points
elements(2).elemTypeAero = [0 0 -d numGaussPoints false] ; % [chordVec1 chordVec2 chordVec3 numGauss  ]
% mass element formulation
elements(2).massMatType = 'lumped' ; %'consistent' ;%
%
% boundaryConds
%
%md The first and unique BC corresponds to a welded condition for a cantilever beam
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%
% initialConds
%
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
% From static reconfiguration
if staticBool
    initialConds.U = matUsStat(:,2);
end
%
% mesh 
%
%mdThe coordinates of the mesh nodes are given by the matrix:
% mesh.nodesCoords = [  zeros(numElements+1,1) zeros(numElements+1,1) (0:(numElements))'*l/numElements] ;
% %mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of nodes that compose the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
% mesh.conecCell = { } ;
% %md then the first welded node is defined with material (M) zero since nodes don't have material, the first element (E) type (the first entry of the `elements` struct), and (B) is the first entry of the the `boundaryConds` struct. For (I) no non-homogeneous initial condition is considered (then zero is used) and finally the node is assigned:
% mesh.conecCell{ 1, 1 } = [ 0 1 1  1 ] ;
% %md Next the frame elements MEBI parameters are set. The frame material is the first material of `materials` struct, then $1$ is assigned. The second entry of the `elements` struct correspond to the frame element employed, so $2$ is set. Finally no BC and no IC is required for this element, then $0$ is used.  Consecutive nodes build the element so then the `mesh.conecCell` is:
% for i=1:numElements,
%   mesh.conecCell{ i+1,1 } = [ 1 2 0  i i+1 ] ;
% end
% Initialize qvect
qvect =  zeros(numElements*2,round(finalTime/dt)+1);
qvect(1:2:end,1) = q0;
pvect =  zeros(numElements*2,round(finalTime/dt)+1);
pvect(1:2:end,1) = p0;
% analysisSettings
% fluid properties
analysisSettings.fluidProps = {rhoFluid; nuFluid; nameFuncVel} ;
%If drag reconfiguration then analysisSettings.geometricNonLinearAero  = true!! also it if is false then the lift direction will be constant
% analysisSettings.geometricNonLinearAero = true;
analysisSettings.geometricNonLinearAero = false;
% time parameters
analysisSettings.finalTime   = finalTime ;
analysisSettings.deltaT      = dt        ;
% numerical parameters
analysisSettings.methodName     =   'newmark' ;
analysisSettings.stopTolIts     =   10        ; 
analysisSettings.stopTolDeltau  =   1e-16     ; 
analysisSettings.stopTolForces  =   1e-15     ; %1e-5
% load settings
analysisSettings.booleanSelfWeight = false ;
% otherParams
otherParams.nodalDispDamping = cu ;
otherParams.problemName      = strcat('HydroForcetest') ;
%otherParams.plotsFormat      = 'vtk' ;
%
% Run ONSAS
%
[ matUsDynLD2D, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 
% %
% % Extract numerical solution
% %linear disp
xdefNum = mesh.nodesCoords(:,1) + matUsDynLD2D(1:6:end, :) ;
ydefNum = mesh.nodesCoords(:,2) + matUsDynLD2D(3:6:end, :) ;
%times vector
times = linspace(0,analysisSettings.finalTime, size(matUsDynLD2D, 2));
% % xlswrite('valsNewmark.xlsx',zdefNum(2,:));
% % xlswrite('valsNewmark.xlsx',qvect(:, 1));
%%
figure(1)
hold on, grid on
yyaxis left
plot(times(1:spanPlot:end), ydefNum(2,1:spanPlot:end)' ,'b-*', 'linewidth', lw,'markersize', ms )
hold on 
%plot(tsC(1:spanPlot:end), UsC(1:spanPlot:end,1) ,'b-o', 'linewidth', lw,'markersize', ms );
ylabel('Uy')
hold on
yyaxis right
ylabel('q')
plot(times(1:spanPlot:end-1), qvect(1, 1:spanPlot:end-1) ,'r-*', 'linewidth', lw,'markersize', ms )
hold on 
%plot( tsC(1:spanPlot:end), UsC(1:spanPlot:end,5) ,'r-o', 'linewidth', lw,'markersize', ms );
xlabel('times')
legend('Uy ONSAS','Uy Monolithic',...
    'q ONSAS','q Monolithic '); %, 'Uz Iterative','q Iterative'
title(sprintf('Lift validation'));

figure(2)
yyaxis left
plot(times(1:spanPlot:end), xdefNum(2,1:spanPlot:end) ,'b-*', 'linewidth', lw,'markersize', ms )
hold on;
%plot(tsC(1:spanPlot:end), UsC(1:spanPlot:end,3) ,'b-o', 'linewidth', lw,'markersize', ms );
ylabel('Ux')
hold on;
yyaxis right
plot(times(1:spanPlot:end-1), pvect(1, 1:spanPlot:end-1) ,'r-*', 'linewidth', lw,'markersize', ms )
hold on;
%plot( tsC(1:spanPlot:end), UsC(1:spanPlot:end,7) ,'r-o', 'linewidth', lw,'markersize', ms );
ylabel('p')
xlabel('times')
legend('Ux ONSAS','Ux Monolithic',...
    'p ONSAS','p Monolithic ');
title(sprintf('Drag validation'));

% %labelTitle = ['formulation = ' elements(2).massMatType 'cq = ' cq 'kq = ' kq] ;
% %title(labelTitle)
% %title(sprintf('formulation = %s , dt = %.3d, cq = %d, kq = %d', elements(2).massMatType,deltat, cq, kq))
% box off;

% elapTime = toc ;
% fprintf( strcat('\n The elapsed execution time was ', num2str(elapTime), ' seconds', '\n' ) );
