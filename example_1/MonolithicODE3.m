accDir = pwd ;
addpath( genpath( [ accDir '/../source'] ) );
%Iteratively solves oscillator equation and non linear Vander Pol equation
loadParametersEx3
initialConds.U = matUsStat(:,2);
if staticBool
    X0 = matUsStat(7,2);
    % Need to include initial moment on Y%%
end

%global VIVBool;
tsB = (0:dt:finalTime)' ;
if ILVIVBool
[tsC, UsC] = ode45(@(t,y) springBeam(t,y,m,ku,cu,kq,cq,l,d,rhoFluid,cL0,cD0,cDi,Cx,Cy,va), ...   
     tsB, [Z0; dZ0; X0; dX0; q0; dq0; p0; dp0], options) ;
else
[tsC, UsC] = ode45(@(t,y) springBeam(t,y,m,ku,cu,kq,cq,l,d,rhoFluid,cL0,cD0,0,0,C2,va), ...   
     tsB, [Z0; dZ0; X0; dX0; q0; dq0], options) ;
end
%[tsC, UsC] = ode45(@(t,y) springBeam(t,y,m,ku,cu,kq,cq,l,d,rhoFluid,cL0,cD0,C2,va), ...   
%     tsB, [Z0; dZ0; X0; dX0; q0], options) ;


% if isThisOctave
%      x2 = tsC(1:spanPlot:end);
%      y3 = UsC(1:spanPlot:end,1);
%      y4 =  UsC(1:spanPlot:end,5);
%      [ax, h3, h4] = plotyy(x2,y3,x2,y4);
%      set ([h3], "color", "b", "marker", "o", "linewidth", lw, "markersize", ms+4);
%      set ([h4], "color", "r", "marker", "o", "linewidth", lw, "markersize", ms+4);
%      ylabel(ax(1), 'Y displacement')
%      ylabel(ax(2), 'q')
%  else 
%      yyaxis left
%      plot(tsC(1:spanPlot:end), UsC(1:spanPlot:end,1) ,'b-o', 'linewidth', lw,'markersize', ms );
%      xlabel('times')
%      ylabel('Ux')
%      hold on
%      yyaxis right
%      %plot( tsA, ysA(:,1));
%      %hold on 
%      %plot( tsB, ysB(:,1));
%      plot( tsC(1:spanPlot:end), UsC(1:spanPlot:end,5) ,'r-o', 'linewidth', lw,'markersize', ms );
%      ylabel('q')
%      title('Lift')
% end
%%
figure(2)
yyaxis left
plot(tsC(1:spanPlot:end), UsC(1:spanPlot:end,3) ,'b-o', 'linewidth', lw,'markersize', ms );
ylabel('Ux')
hold on;
yyaxis right
plot( tsC(1:spanPlot:end), UsC(1:spanPlot:end,7) ,'r-o', 'linewidth', lw,'markersize', ms );
ylabel('p')
xlabel('Time')
title('Drag')

figure(1), grid on, hold on
yyaxis left
plot(tsC(1:spanPlot:end), UsC(1:spanPlot:end,1) ,'b-o', 'linewidth', lw,'markersize', ms );
hold on
ylabel('Ux')
yyaxis right
plot( tsC(1:spanPlot:end), UsC(1:spanPlot:end,5) ,'r-o', 'linewidth', lw,'markersize', ms );
ylabel('q')
hold on

%% FFT
sig = UsC(:,3);
Ns = length(sig);
xhat = fft(sig,Ns);
PSD = xhat.*conj(xhat)/Ns ;
freq = 1/(dt*Ns)*(0:Ns);
L = 1:floor(Ns/1000); % Only plot the first half of frequencies
% figure(50) % plot FFT
% plot(freq(L), PSD(L))

% figure(5)
% plot(UsC(1:spanPlot:end,3), UsC(1:spanPlot:end,1) ,'-', 'linewidth', lw,'markersize', ms );
% hold on
% plot([-delta -delta], [min(UsC(1:spanPlot:end,1)) max(UsC(1:spanPlot:end,1))])
% hold on;
% xlabel('Ux')
% ylabel('Uz')

%IterativeODELift(d, l, m, ku, cu, cq, kq, finalTime, dt,Z0, dZ0, ddZ0, q0, dq0, ddq0);

