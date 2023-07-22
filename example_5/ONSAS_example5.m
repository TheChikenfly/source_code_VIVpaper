%
%md# WOM coral example
%
%close all; clear all ;
tic ;
%
% declare global variables
%
global qvect; %VIV boolean is called inside hydroFrameForces 
global VIVBool; %VIV boolean is called inside hydroFrameForces
global constantLiftDir; %constantLiftDir is called inside hydroFrameForces
global uniformUdot; %constantLiftDir is called inside hydroFrameForces
global AMBool;
VIVBool = true;
if VIVBool
  constantLiftDir = false; uniformUdot = false; % constantLiftDir = true for validation !
end
global exportFirstMatrices;
exportFirstMatrices = false; % saves matrix in output
global modalAnalysisBoolean
modalAnalysisBoolean = false; % writes modes vtk
% add ONSAS path
%
accDir = pwd ;
addpath( genpath( [ accDir '/../../../ONSAS.m/src'] ) );
addpath( genpath( [ accDir '/../source'] ) );
%
% BE SURE THAT:
% Consistent not lumped, 
%
loadParametersCoral
%betavect = [1.8751 4.6941 7.8547 10.9955 14.13717];
%freqanalytical = (betavect.^2)*sqrt(E*I/(ms))/(2*pi*l^2); % ma contained in ms
%
% materials
%
materials.hyperElasParams = [ E nu ] ;
materials.density         = rho      ;
% Corot element for large displacements
materials.hyperElasModel  = '1DrotEngStrain' ; 
%
% elements
%
% node
elements(1).elemType = 'node' ;
%
% hydro frame
%
elements(2).elemType = 'frame' ;
% cross section params
elements(2).elemCrossSecParams{1,1} = 'circle' ;
elements(2).elemCrossSecParams{2,1} = d          ;
% hydro cross-section props
numGaussPoints  = 4 ;
elements(2).aeroCoefs   = {nameDragFunc; nameLiftFunc; [] }   ;
%  chord vector and gauss points
elements(2).elemTypeAero = [0 0 -d numGaussPoints false] ; % [chordVec1 chordVec2 chordVec3 numGauss  ]
% mass element formulation
elements(2).massMatType = 'consistent' ; %'lumped'; %
%
% boundaryConds
%
%md The first and unique BC corresponds to a welded condition for a cantilever beam
%{}
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%
% initialConds
%
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
%
% mesh 
CoralMesh
%meshgeometry
%meshgeom
myCell2Mat( mesh.conecCell );% To visualise connection
% Initialize qvect
qvect =  zeros(numElements*2,round(finalTime/dt)+1);
qvect(1:2:end,1) = (2*rand(numElements, 1)-1)*0.001 ;
% analysisSettings
%
% fluid properties
analysisSettings.fluidProps = {rhoFluid; nuFluid; nameFuncVel} ;
%If drag reconfiguration then analysisSettings.geometricNonLinearAero  = true!! also it if is false then the lift direction will be constant
analysisSettings.geometricNonLinearAero = true;
% time parameters
analysisSettings.finalTime   = finalTime ;
analysisSettings.deltaT      = dt        ; %/400 ;
% numerical parameters
analysisSettings.methodName     =   'newmark' ;
analysisSettings.stopTolIts     =   50        ; %;10
analysisSettings.stopTolDeltau  =   1e-10      ; 
analysisSettings.stopTolForces  =   1e-5      ; 
% load settings
analysisSettings.booleanSelfWeight = false ;
% % Folder path
%
% analysisSettings
%
% otherParams
otherParams.nodalDispDamping = cu ;
otherParams.problemName      = strcat('ONSAS_Ex4_Coral') ;
otherParams.plots_format      = 'vtk' ;
%
% Run ONSAS
%
%sprintf('matUsCoral_FT=%d_vmax=%.1f_freq=%.1f_epsilon=%f_A=%.d.mat',finalTime, vwindMax, freq, epsilon, A)
[ matUs, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 
% Extract numerical solution
%linear disp first node
zdefNum = mesh.nodesCoords(:,3) + matUs(5:6:end, :) ;
xdefNum = mesh.nodesCoords(:,1) + matUs(1:6:end, :) ;
%linear disp last node
zdefNum = mesh.nodesCoords(:,3) + matUs(end-1:6:end, :) ;
xdefNum = mesh.nodesCoords(:,1) + matUs(end-5:6:end, :) ;
%for all nodes
zdefNumall = matUs(5:6:end, :) ;% Z 
xdefNumall = matUs(1:6:end, :) ;% X 
zdefNumallbis = mesh.nodesCoords(:,3) + matUs(5:6:end, :) ;
xdefNumallbis = mesh.nodesCoords(:,1) + matUs(1:6:end, :) ;
%
save(sprintf('matUsconstantBeam_FT=%d_vmax=%.5f_epsilon=%f_A=%.d.mat',finalTime, vwindMax, epsilon, A), 'matUs', 'mesh')
%times vector
times = linspace(0,analysisSettings.finalTime, size(matUs, 2));
length(times);
% 
% spanPlot = 1 ; lw = 2.0 ; ms = 11 ; plotfontsize = 20 ;
% figure(2)
% hold on, grid on
% yyaxis left
% plot(times(1:spanPlot:end), zdefNum(2,1:spanPlot:end) ,'b-*', 'linewidth', lw,'markersize', ms )
% ylabel('Uz')
% hold on
% yyaxis right
% % q last element
% plot(times(1:spanPlot:end-1), qvect(end-1, 1:spanPlot:end) ,'r-*', 'linewidth', lw,'markersize', ms )
% % plot(times(1:spanPlot:end-1), qvect(1, 1:spanPlot:end) ,'r-*', 'linewidth', lw,'markersize', ms )
% xlabel('times')
% ylabel('q')
% 
% figure(5)
% plot(xdefNum(2,end-100:spanPlot:end), zdefNum(2,end-100:spanPlot:end) ,'-', 'linewidth', lw,'markersize', ms );
% %plot(xdefNum(2,1:spanPlot:end), zdefNum(2,1:spanPlot:end) ,'-', 'linewidth', lw,'markersize', ms );
% hold on
% xlabel('Ux')
% ylabel('Uz')
% %labelTitle = ['formulation = ' elements(2).massMatType 'cq = ' cq 'kq = ' kq] ;
% %title(labelTitle)
% %title(sprintf('formulation = %s , dt = %.3d, cq = %d, kq = %d', elements(2).massMatType,deltat, cq, kq))

elapTime = toc ;
fprintf( strcat('\n The elapsed execution time was ', num2str(elapTime), ' seconds', '\n' ) );
%-----------------------------------------------------------------------------------------------------
