% Fork N=1
close all, clear all ;
global vwindMax;global qvect; global pvect; global VIVBool; global constantLiftDir; global uniformUdot; 
global AMBool; global fluidFlowBool; global finalTime;
global ILVIVBool;
% add onsas block
onsas_path = getenv('ONSAS_PATH');
addpath( genpath( [ pwd '/../shared'] ) );
if isempty( onsas_path )
  % default alexandre path
  onsas_path = ['C:\Users\alvilh\ONSAS.m\sharedFiles\ONSAS.m'] ;
end
addpath( genpath( onsas_path ) );
%
loadParametersCoralEx5
% Static
[matUsStat] = ONSAS( materials, elements_static, boundaryConds, initialConds_static, mesh, analysisSettings_stat, otherParams_stat ) ; 
%% Dynamic
global Vrelvect;
VIVBool = true; ILVIVBool = true;
constantLiftDir = false; uniformUdot = false; 
AMBool = false; fluidFlowBool = false;
finalTime = finalTime_Dynamic;
Yrms = zeros(1, NR); Urvect = zeros(1, NR);
for indvel = 17:NR 
    tic;
    vwindMax = velocities(indvel)
    initialConds_dynamic.U= matUsStat(:,indvel+1);
    Urvect(1,indvel) = vwindMax/(d*4.05);
    Vrelvect = zeros(6*(numElements+1), finalTime/dt);
    qvect =  zeros(numElements*2,round(finalTime/dt)+1); pvect =  zeros(numElements*2,round(finalTime/dt)+1);
    qvect(1:2:end,1) = 0.001;%2*(2*rand(numElements, 1)-1); 
    pvect(1:2:end,1) = 0.001;%2*(2*rand(numElements, 1)-1); 
    otherParams_dyn.nodalDispDamping = 0;
    [ matUsDynLD2D, ~ ] = ONSAS( materials, elements_dynamic, boundaryConds, initialConds_dynamic, mesh, analysisSettings_dyn, otherParams_dyn ) ; 
    ydefNumall = matUsDynLD2D(3:6:end, :); xdefNumall = matUsDynLD2D(1:6:end, :) ; zdefNumall = matUsDynLD2D(5:6:end, :) ;
    Yrms(1, indvel) = rms(ydefNumall(end, floor(end/2):end)./d);
    save(sprintf('Ysolutions\\YCoralN=1VR_va=%.3f_Nelem=%d_FT%d_dt=%.4f.mat', vwindMax, numElements, finalTime,dt), 'ydefNumall', 'xdefNumall', 'zdefNumall', 'Vrelvect')
    elapTime2 = toc ;
    fprintf( strcat('\n The elapsed execution time was ', num2str(elapTime2/60), ' minutes or ', num2str(elapTime2/3600), ' hours',  '\n' ) );
end
%% 
nelem = [50 100 200 500]% 5 10 20 1000 [10 20 50 100 200 500]
%Yrmsmidnode = zeros(1,length(nelem)); 
timeHistoryInt = zeros(1,length(nelem));
load(sprintf('Ysolutions\\YCoralN=0MeshStudy_va=%.3f_Nelem=%d_FT%d_dt=%.4f.mat', 0.280, 1000, 0.5,dt))
Nref = 1000;
%ref = ydefNumall(Nref/2+1, end/2:end);
refx = xdefNumall(:, 1:end);
refy = ydefNumall(:, 1:end);
for j = 1:length(nelem)
    ratio = Nref/nelem(j);
    load(sprintf('Ysolutions\\YCoralN=0MeshStudy_va=%.3f_Nelem=%d_FT%d_dt=%.4f.mat', 0.280, numElements, 0.5,dt))
    signalx = xdefNumall(:, 1:end); % all signal 
    signaly = ydefNumall(:, 1:end); % all signal 
    s = linspace(0,1,nelem(j)+1);
    Yrms= rms(signaly(:, end/2:end)./d, 2);
    %Yrmsmidnode(1,j) = (rms(signaly(midNode,:)./d, 2) - rms(ref./d, 2))/rms(ref./d, 2);
    %midnodeHistory(1,j) = vecnorm((signaly(:,:) - ref(1:ratio:end,:))/ref(1:ratio:end,:));
    %nodal_differencex = vecnorm(abs(signalx(:,:)' - refx(1:ratio:end,:)'),2);
    %nodal_differencey = vecnorm(abs(signaly(:,:)' - refy(1:ratio:end,:)'),2);
    nodal_differencex = sum(abs(signalx(:,:)' - refx(1:ratio:end,:)'),1);
    nodal_differencey = sum(abs(signaly(:,:)' - refy(1:ratio:end,:)'),1);
    nodal_difference = sum( nodal_differencex + nodal_differencey ) ;
    %intUref = sum( vecnorm(refx(1:ratio:end,:)',2) + vecnorm(refy(1:ratio:end,:)',2) ) ;
    intUref = sum( sum(abs(refx(1:ratio:end,:)'),1) + sum(abs(refy(1:ratio:end,:)'),1) ) ;
    timeHistoryInt(1,j) = nodal_difference/intUref ;
    figure(3)
    plot(s, Yrms)
    xlabel('Yrms'); ylabel('z')
    hold on     
end
%%
legend('N = 10','N = 20','N = 50','N = 100','N = 200', 'N = 500', 'N = 100')
figure(5)
loglog(nelem, timeHistoryInt, 'k-o')
%plot(nelem, timeHistoryInt, 'k-o')
xlabel('Number of elements'); ylabel('Relative error')
title('Relative error on the displacement over 1000 time steps')
% 
% figure(12)
% plot(Urvect, Yrms, 'k-o')
% xlabel('St Ur'); ylabel('Yrms/d');  