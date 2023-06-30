% Main scrip comparing MOnlitic, Iterative and Onsas for a beam with lift
% and no drag
close all; clear all;
% add onsas block
onsas_path = getenv('ONSAS_PATH');
addpath( genpath( [ pwd '/../shared'] ) );
if isempty( onsas_path )
  % default alexandre path
  onsas_path = ['C:\Users\alvilh\ONSAS.m\sharedFiles\ONSAS.m'] ;
end
addpath( genpath( onsas_path ) );
% load parameters script
%
loadParametersEx3
global ILVIVBool;
% MUST have CONSTANT velocity va = 0.004 in windVelUniform, WOMV2 
% constantLiftDir = false; uniformUdot = true ;

ONSAS_example_3_Static
ONSAS_example_3

MonolithicODE3

save(sprintf('resultsEx1V2FT=%.1f,dt=%.6f.mat', finalTime, dt), 'times', 'xdefNum', 'qvect', 'pvect', 'ydefNum', 'UsC', 'spanPlot', 'finalTime', 'dt', 'staticBool');
%IterativeODE3
%% Plot 
figure(1)
hold on, grid on
yyaxis left
plot(times(1:spanPlot:end), ydefNum(2,1:spanPlot:end)' ,'b-*', 'linewidth', lw,'markersize', ms )
hold on 
plot(tsC(1:spanPlot:end), UsC(1:spanPlot:end,1) ,'b-o', 'linewidth', lw,'markersize', ms );
ylabel('Uy')
hold on
yyaxis right
ylabel('q')
plot(times(1:spanPlot:end-1), qvect(1, 1:spanPlot:end-1) ,'r-*', 'linewidth', lw,'markersize', ms )
hold on 
plot( tsC(1:spanPlot:end), UsC(1:spanPlot:end,5) ,'r-o', 'linewidth', lw,'markersize', ms );
xlabel('times')
legend('Uy ONSAS','Uy Monolithic',...
    'q ONSAS','q Monolithic '); %, 'Uz Iterative','q Iterative'
title(sprintf('Lift validation'));

figure(2)
yyaxis left
plot(times(1:spanPlot:end), xdefNum(2,1:spanPlot:end) ,'b-*', 'linewidth', lw,'markersize', ms )
hold on;
plot(tsC(1:spanPlot:end-1), UsC(1:spanPlot:end-1,3) ,'b-o', 'linewidth', lw,'markersize', ms );
ylabel('Ux')
hold on;
yyaxis right
plot(times(1:spanPlot:end-1), pvect(1, 1:spanPlot:end-1) ,'r-*', 'linewidth', lw,'markersize', ms )
hold on;
plot( tsC(1:spanPlot:end-1), UsC(1:spanPlot:end-1,7) ,'r-o', 'linewidth', lw,'markersize', ms );
ylabel('p')
xlabel('times')
legend('Ux ONSAS','Ux Monolithic',...
    'p ONSAS','p Monolithic ');
title(sprintf('Drag validation'));
% tdfixed = drag direction is constant
%V2: force tol was 1e-5, we decrease it to 1e-7 odetol from 1e-6 to 1e-8
%% Error
coef = 0.99
errorX = norm(abs(xdefNum(2,1:coef*end)-UsC(1:coef*end,3)'))/norm(UsC(1:coef*end,3))
errorY = norm(abs(ydefNum(2,1:coef*end)-UsC(1:coef*end,1)'))/norm(UsC(1:coef*end,1))
errorQ = norm(abs(qvect(1,1:coef*end)-UsC(1:coef*end,5)'))/norm(UsC(1:coef*end,5))
errorP = norm(abs(pvect(1,1:coef*end)-UsC(1:coef*end,7)'))/norm(UsC(1:coef*end,7))