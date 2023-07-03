% Validation IL VIV with Trim (Exp) and Gao (WOM)
close all, clear all; %clearvars -except nelem;
nelem = 1000;
%
global pretension_strain;
global vwindMax;global qvect; global pvect; global VIVBool; global constantLiftDir; global uniformUdot; 
global AMBool; global fluidFlowBool; global finalTime; global pretension_strain; global massratio;
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
loadParametersEx3
% Static
[matUsStat] = ONSAS( materials, elements_static, boundaryConds, initialConds_static, mesh, analysisSettings_stat, otherParams_stat ) ; 
% Dynamic
VIVBool = true; ILVIVBool = true;
constantLiftDir = false; uniformUdot = false; % Lift along Vrelperp
AMBool = false; fluidFlowBool = false;
for ps = [2.4e-4]
    tic;
    pretension_strain = ps
    initialConds.U = matUsStat(:,2);
    qvect =  zeros(numElements*2,round(finalTime/dt)+1);
    qvect(1:2:end,1) = 0.001;%0.001*(2*rand(numElements, 1)-1); 
    pvect =  zeros(numElements*2,round(finalTime/dt)+1);
    pvect(1:2:end,1) = 0.001;%0.001*(2*rand(numElements, 1)-1); 
    otherParams.nodalDispDamping = 0;
    [ matUsDynLD2D, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 
    ydefNumall = matUsDynLD2D(3:6:end, :); xdefNumall = matUsDynLD2D(1:6:end, :) ;
    save(sprintf('YEx3NelemV16.4=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,pretension_strain), 'ydefNumall', 'xdefNumall')
    elapTime2 = toc ;
    fprintf( strcat('\n The elapsed execution time was ', num2str(elapTime2/60), ' minutes = ', num2str(elapTime2/3600), ' hours',  '\n' ) );
%end
end
%%
%Plot beam reconfiguration
for i = 1:1
    figure(2)
    plot(matUsStat(1:6:end,i),matUsStat(5:6:end,i)+mesh.nodesCoords(:,3))
    %axis equal
    xlabel('x'); ylabel('z')
    hold on 
end
%%
nelem = [5 10 20 50 100 200 500]%1000 [10 20 50 100 200 500]
%Yrmsmidnode = zeros(1,length(nelem)); 
timeHistoryInt = zeros(1,length(nelem));
load(sprintf('YEx3NelemV16.3bis=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', 1000, finalTime, vwindMax,ps))
Nref = 1000;
%ref = ydefNumall(Nref/2+1, end/2:end);
refx = xdefNumall(:, 1:end);
refy = ydefNumall(:, 1:end);
for j = 1:length(nelem)
    ratio = Nref/nelem(j);
    load(sprintf('YEx3NelemV16.3bis=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', nelem(j), finalTime, vwindMax,ps))
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
%-----------------------------------------------------------------------------------------------------
