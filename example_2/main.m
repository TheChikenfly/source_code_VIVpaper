% Validation Leclercq
close all, clear all ;
global vwindMax;global qvect; global VIVBool; global constantLiftDir; global uniformUdot; 
global AMBool; global fluidFlowBool; global finalTime; global velocities;
% add onsas block
onsas_path = getenv('ONSAS_PATH');
addpath( genpath( [ pwd '/../shared'] ) );
if isempty( onsas_path )
  % default alexandre path
  onsas_path = ['C:\Users\alvilh\ONSAS.m\sharedFiles\ONSAS.m'] ;
end
addpath( genpath( onsas_path ) );
%%
ONSAS_static
[matUsStat] = ONSAS( materials, elements_static, boundaryConds, initialConds_static, mesh, analysisSettings_stat, otherParams_stat ) ;
%% Dynamic
VIVBool = true;
constantLiftDir = true; uniformUdot = false;
AMBool = false; fluidFlowBool = true;
finalTime = finalTime_Dynamic;
Yrms = zeros(1, NR); Urvect = zeros(1, NR);
%Av = 96:NR; Bv = 96:10:NR;
Av = 174:NR; Bv = 96:10:NR;
idx=ismember(Av,Bv);
Av(idx)=[]; % remove idx from Av
for indvel = [141 142] %Av % 55
    tic;
    vwindMax = velocities(indvel)
    Urvect(1,indvel) = (St*(l^2)*vwindMax/d)*sqrt((ms+ma)/(E*I));
    initialConds_dynamic.U = matUsStat(:,indvel+1);
    qvect =  zeros(numElements*2,round(finalTime/dt)+1);
    qvect(1:2:end,1) = 2*(2*rand(numElements, 1)-1); %(2*rand(numElements, 1)-1)*0.001 ;
    [ matUsDynLD2D, ~ ] = ONSAS( materials, elements_dynamic, boundaryConds, initialConds_dynamic, mesh, analysisSettings_dyn, otherParams_dyn ) ; 
    ydefNumall = matUsDynLD2D(3:6:end, :); xdefNumall = matUsDynLD2D(1:6:end, :) ;
    Yrms(1, indvel) = rms(ydefNumall(end, floor(end/2):end)./d);
    %save(sprintf('Ysolutions\\YflowBoolTrueV3_va=%.3f_Nelem=%d_FT%d_dt=%.3f.mat', vwindMax, numElements, finalTime,dt), 'ydefNumall')
%     save(sprintf('Ysolutions\\YflowBoolTrueV2_va=%.3f_Nelem=%d_FT%d.mat', vwindMax, numElements, finalTime), 'ydefNumall')
    save(sprintf('Ysolutions\\YflowBoolTrueV6Examine_va=%.4f_Nelem=%d_FT%d_dt=%.4f.mat', vwindMax, numElements, finalTime,dt), 'ydefNumall')
    elapTime2 = toc ; fprintf( strcat('\n The elapsed execution time was ', num2str(elapTime2/60), ' minutes or ', num2str(elapTime2/3600), ' hours',  '\n' ) );
end
figure(12)
plot(Urvect, Yrms, 'k-o')
xlabel('St Ur'); ylabel('Yrms/d'); 
%% 
% times = linspace(0,finalTime, size(ydefNumall, 2));
% zdefNumlast = matUsDynLD2D(end-1, :) ;
% figure()
% plot(times, ydefNumlast)
% Yrms = rms(ydefNumlast(floor(end/2):end)./d)
% for i = 1:10:300
% figure()
% plot(times, zdefNumlast)
% hold on 
% end

% Plot beam reconfiguration
for i = 1:30
    figure(2)
    plot(matUsStat(1:6:end,i),matUsStat(5:6:end,i)+mesh.nodesCoords(:,3))
    %axis equal
    xlabel('x'); ylabel('z')
    hold on 
end


%for va=0.0045:0.001:0.0045; % 15 simulations Ur form 0.08 to 2.4,pick at va = 0.005
% for ne = [15 30 50 100 200]
%     tic;
%     numElements = ne
% %     vwindMax = va
% %     loadParametersValidation
% %     sprintf('ZUniform_udotFlowNode_locktdtl_Ur=%.2f.mat', UrLeclercq)
%     ONSAS_example5
%     %sprintf('Zsolutions\\ZLinear_va=%.2f_Nelem=%d_FT%d.mat', vwindMax, numElements,finalTime)
%     save(sprintf('Zsolutions\\Zwinstep_va=%.2f_Nelem=%d_FT%d_Nstep=20_coefacc=0.2.mat', vwindMax, numElements,finalTime), 'zdefNumall')
%     elapTime2 = toc ;
%     fprintf( strcat('\n The elapsed execution time was ', num2str(elapTime2/60), ' minutes or ', num2str(elapTime2/3600), ' hours',  '\n' ) );
% end

%{
for n = [50]% 25 simulations Ur form 0.08 to 4
for ft = [10 30 50 100]% 25 simulations Ur form 0.08 to 4
    tic;
    finalTime = ft
    numElements = n
    ONSAS_example5
    save(sprintf('ZlowUrNelem=%d_FT=%d.mat', numElements, finalTime), 'zdefNumall')
    elapTime2 = toc ;
    fprintf( strcat('\n The elapsed execution time was ', num2str(elapTime2/60), ' minutes or ', num2str(elapTime2/3600), ' hours',  '\n' ) );
end
end
%save(sprintf('ZsolV4Nelem=%d_FT=%d_E=%.0s_dt=%.4f_vmax=%.1f_cu=%d.mat', numElements, finalTime,E,dt,vwindMax,cu), 'zdefNumall')
%-----------------------------------------------------------------------------------------------------
%}