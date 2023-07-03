% Post processing signal Ydefnum
loadParametersCoral
Y_rmsvectExp = zeros(numElements+1, NR); UrExp = zeros(1, NR); fExp = zeros(1, NR); gamma = zeros(1, NR);
%velocities = 0.06:0.02:1.4;
velocities = [0.01:0.01:0.01+8*0.01 0.12:0.02:1.4]; %NR=74
NR = length(velocities);
%velocities = 0.06:0.01:1.4;
l = 0.14;
stop
FY_tip=NaN(1002,NR);PY_tip=NaN(1002,NR);
for j = 1:NR%[1 18 23 30 NR]%1:NR% [8:3:17] % idx(1, 1:NR)%[1 21 39 56 NR]
    j
   %ind = (sp-60)/5;
   U = velocities(j); UrExp(1, j) = ((l^2)*U/d)*sqrt((ms+ma)/(E*I)); % U/(d*f1Exp);%
%    if (j <21); dt = 0.01; else; dt = 0.005;end
%    load(sprintf('Ysolutions\\YCoralN=0_va=%.3f_Nelem=%d_FT%d_dt=%.4f_cu=%.4f.mat', U, numElements, finalTime, dt, cu));
   %Correction on windvel function that was wrong.
   %if (j <18); dt = 0.005; finalTime=5; else; dt = 0.0005;finalTime=2;end
   %load(sprintf('Ysolutions\\YCoralN=0V2_va=%.3f_Nelem=%d_FT%d_dt=%.4f_cu=%.5f.mat', U, numElements, finalTime,dt, cu));
   % epsilony = 0.04 not 0.3. finalTime=2; NR = 8
   %load(sprintf('Ysolutions\\YCoralN=0V3_va=%.3f_Nelem=%d_FT%d_dt=%.4f_cu=%.5f.mat', U, numElements, finalTime,dt, cu));
   % 'innerDragCoefCircular'; cu = 0.00024, Ay = 12; epsilony = 0.04; 
%    load(sprintf('Ysolutions\\YCoralN=0V4_va=%.3f_Nelem=%d_FT2_dt=%.4f_cu=%.5f.mat', U, numElements, dt, 0.00024));
   % 'innerDragCoefCircular'; cu = 0.000006, Ay = 12; epsilony = 0.3; velocities = 0.06:0.02:1.4; % Best
   if  UrExp(1, j)> 4.5
       load(sprintf('Ysolutions\\YCoralN=0V5_va=%.3f_Nelem=%d_FT2_dt=%.4f_cu=%.5f.mat',U, numElements, dt, cu));        
       %l=0.14;%
       %UrExp(1, j) = ((l^2)*U/d)*sqrt((ms+ma)/(E*I))
   %same, %48 dots
   elseif  UrExp(1, j)<= 4.5 % Finer dUr
       load(sprintf('Ysolutions\\YCoralN=0V5bis_va=%.3f_Nelem=%d_FT2_dt=%.4f_cu=%.5f.mat', U, numElements, dt, cu));
       %UrExp(1, j) = U/(d*f1Exp)
   end
   % Changing epsilony
   %load(sprintf('Ysolutions\\YCoralN=0V6_va=%.3f_Nelem=%d_FT%d_dt=%.4f_cu=%.5f_Ay=%d_epsy=%.4f.mat', U, numElements, finalTime,dt, cu, 12, epsilony));
   % epsilony = 0.6
   %load(sprintf('Ysolutions\\YCoralN=0V7_va=%.3f_Nelem=%d_FT%d_dt=%.4f_cu=%.5f_epsy=%.4f.mat', U, numElements, finalTime,dt, cu, 0.6));
   % Printing drag reconfiguration
   load(sprintf('Ysolutions\\YCoralN=0V8_va=%.3f_Nelem=%d_FT%d_dt=%.4f.mat', U, numElements, 2,dt))
   sig = ydefNumall(end, floor(end/2):end)./d;
   Y_rmsvectExp(:,j) = rms( ydefNumall(:, floor(end/2):end)./d,2 ); % Computes at all nodes
   %fExp(1,j) = Yfft(sig,length(sig), dt, true, 0.1, Etpe, I, ms, ma, l); % 
   %l = 0.15;
   Nt = length(sig); % 2001
%    fi = Yfft(sig,Nt, dt, true, 0.1, UrExp(1,j), Etpe, I, ms, ma, l);
%     plotspectro(UrExp(1,j), j, fi, '')
    %fw0i = (St*U/d);
    fw0i = 1;
    [fY,pY] = psdlec(sig,1,1/dt);
    size(fY.'); %= size(pY)  %size(FY_tip(:,j));
    FY_tip(:,j)=fY.'./fw0i; PY_tip(:,j)=pY.';
    %[~, gamma(1,j)] = lemniscate(fExp(1,j), dt, xdefNumall, ydefNumall, false, UrExp(1, j),d);
%     figure(3)
%     xdefNumall(end,end)
%     s=mesh.nodesCoords(:,3) + zdefNumall(:,end);
%     ur=((l^2)*U/d)*sqrt((ms+ma)/(E*I))
%     plot(ur+xdefNumall(:,end)./d,s./d,'k-','linewidth',2)
%     hold on 
%     ylabel('s'); xlabel('x/d');title('In-line displacement')
end
%% plots
%NR = 68
figure(12)
%Y_rmsvectExp(Y_rmsvectExp==0) = []; UrExp(UrExp==0) = [];
%fig = openfig('ExpRms.fig'); hold on
plot(UrExp(1:NR), Y_rmsvectExp(end,1:NR), 'k-')
xlabel('Ur =  U/(d f1)'); ylabel('Yrms/d'); title('Yrms amplitude');
legend('Experimental','ONSAS')
relative_error = [(0.496-0.399)/0.496 ; (0.721-0.714)/0.714]; % Peak 1 and 2
%
betavect = [1.87510 4.69409 7.857476]; freqStruct = (betavect.^2)*sqrt(E*I/(ms+ma))/(2*pi*l^2)  % In water  (beta^2 / L^2) sqrt(EI / rho A)
PY_tipmaxtimeslog = max(10*log(PY_tip),[],1); % Max for each time: vector length times
PY_tipmaxtimes = max(PY_tip,[],1); % Max for each time: vector length times
FY_tip(:,1) = zeros(1002, 1);
PY_tip(:,1)= zeros(1002, 1);
figure()
%h=pcolor(UrExp,FY_tip./freqStruct(1),10*log(PY_tip) ); 
h=pcolor(UrExp,FY_tip,10*log(PY_tip) ); 
%h=pcolor(UrExp,FY_tip,PY_tip./repmat(PY_tipmaxtimes,length(FY_tip),1));  
caxis([-200 1]);
%ylim([0 60]); xlim([0 50]);
%ylim([0 19]); xlim([0 50]);
ylim([0 50]); xlim([0 34]);
hold on
c = colorbar; %c.Label.String = 'PSD';
set(h,'edgecolor','none','facecolor','interp');
%yline(freqStruct(1)/freqStruct(1), '-.');yline(freqStruct(2)/freqStruct(1), '-.');yline(freqStruct(3)/freqStruct(1), '-.');
%plot(UrExp, St*(UrExp*sqrt((E*I)/(ms+ma))*d/l^2/d)/freqStruct(1), 'k--','MarkerSize',2) % Plot Strouhal frequency
plot(UrExp, St*(UrExp*sqrt((E*I)/(ms+ma))*d/l^2/d), 'k--','MarkerSize',2) % Plot Strouhal frequency
%xlabel('Ur');ylabel('f/fw0');title('ONSAS PSD of Y at the tip')
% figure() % fExp
% plot(UrExp(1:NR), fExp(1:NR), 'r*')
% xlabel('Ur =  U/(d f1)'); ylabel('f');
%figure()
%plot(ydefNumall(end, :)./d)
%% Static reconfiguration
% first run only the static with velocities = 0.06:0.02:1.4;
velocities = 0.06:0.02:1.4; NR = length(velocities )
csvname = 'Ex4Reconfiguration.csv';
t = readtable(csvname); 
for i = 1%:2:NR
    U = velocities(i); Ur = ((l^2)*U/d)*sqrt((ms+ma)/(E*I))
    x = matUsStat(1:6:end,i);
    z = mesh.nodesCoords(:,3) + matUsStat(5:6:end,i);
    figure(10)
    plot(Ur+x./d,z./d,'k-','linewidth',2)
    hold on
    NColumn = size(t,2)
    t.(NColumn+1) = Ur+x./d; t.Properties.VariableNames{NColumn+1} = sprintf('ur%.0fx',Ur);
    t.(NColumn+2) = z./d;t.Properties.VariableNames{NColumn+1} = sprintf('ur%.0fz',Ur);
end
writetable(t, 'Ex4Reconfigurationbis.csv')
%% SPACE-TIME EVOLUTION
s = linspace(0,1,numElements+1);
Yrmsallmaxtimes = max(Y_rmsvectExp,[],1); % Max for each time: vector length times
figure()
h=pcolor(UrExp,s,Y_rmsvectExp(:,:)); % Yrms not normalized
%h=pcolor(UrExp,s,Y_rmsvectExp./repmat(Yrmsallmaxtimes,length(s),1)); % Yrms normalized
c = colorbar; %c.Label.String = 'Yrms/d';
%ylabel('s') ; xlabel('Ur');title('Y space-time evolution, ONSAS')
set(h,'edgecolor','none','facecolor','interp');
%% Shape of the beam at the end
figure()
Ns = length(xdefNumall(:,end));
%plot(-xdefNumall(:,end)./d,s,'linewidth',2)
%loc = 3.92/30*Ns;
loc = 0.99*Ns;
if ~isempty(xdefNumall)
%    subplot(1,2,1);
    figure(3)
    plot(-xdefNumall(:,loc:loc)./d,s.*0.15/0.005,'k-','linewidth',2)
    ylabel('s'); xlabel('x/d');title('In-line displacement')
%     subplot(1,3,2);
%     plot(-xdefNumall(:,loc-10:loc)./d+mean(xdefNumall(:,end/2:end),2)./d,s,'linewidth',2)
%     ylabel('s'); xlabel('x/d');
%     subplot(1,2,2);
%     plot(ydefNumall(:,loc-200:loc-200)./d,s,'linewidth',2)
    ylabel('s'); xlabel('z/d');title('Cross-flow displacement')
else 
    plot(zdefNumall(:,loc-100:loc)./d,s,'linewidth',2)
end
%%  FFT and signal     
%{
% FT = fft(zmid(N:end));
%from steve brunton (weird it offsets on freq axis)
xhat = fft(ztiptrunk(1:end),Ns); %same as xhat = fft(zmid(1:end));
PSD = xhat.*conj(xhat)/Ns ;
freq = 1/(dt*Ns)*(0:Ns);
L = 1:floor(Ns/2); % Only plot the first half of frequencies
figure(10)
plot(freq(L), PSD(L))
title('FFT'); xlabel('f(Hz)')
hold on
%}
%% Spectrogram 
Nfft = 2^nextpow2(Ns);
Nspec = 512;
Noverlap = Nspec/2;
% dur = 1/(UR(2)-UR(1)); % For fft
[~,~,~,pxx,fc,uc] = spectrogram(ymid, Nspec, Noverlap, Nspec, 1/dt,'MinThreshold',-135); 
Fvect = fc(pxx>0)./pi; %careful spectogram function in Matlab normalizes frquencies by pi 
Uri = uc(pxx>0);
for i = 1:length(fc(pxx>0))
    fw0i = Uri(i)/((l^2)*sqrt(mtot/(E*I))); % without added mass
    Fvect(i) = Fvect(i)/fw0i; 
end
figure()
imagesc(uc(pxx>0),fc(pxx>0)./pi, 10*log10(pxx(1:end,1:end)+eps)); 
h = colorbar;
h.Label.String = 'Power/frequency (dB/Hz)';
ylabel('PSD') ;
% Leclrecq: 
%FY_tip=NaN(Nt/4+1,NR);
%[FY(end,:),PY(end,:)] = psdlec((Y(Ntip,floor(end/2)+1:end)),1,1/dt);
%FY_tip(:,ku)=FY(end,:).';
%h=pcolor(fs_range,FY_tip,10*log(PY_tip)); % To plot 
%% Space time evolution
% windowlenght = 50; %floor(Ns/10); 
% nlast = Ns - mod(Ns, windowlenght);
% Zrmsw = sqrt(mean(reshape(zdefNumall(midNode,1:nlast), windowlenght, []) .^ 2))/d; % Counting every data only once
% movrms = sqrt(movmean(zdefNumall(midNode,:).^ 2, windowlenght))/d;  % Sliding window on the first peak
% size(movrms)
% figure()
% plot(times, movrms)
% 
time1 = floor(20/dt); % which second ?
time2 = time1 + floor(1/dt); % + 1 second
figure()
c = colorbar; c.Label.String = 'Zrms/d';
Nspace = 100;
%h=pcolor(times,s,Zrms./repmat(Zrmsallmaxtimes,length(s),1)); % Zrms normalized
h=pcolor(times(time1:time2),s.*(l/d),zdefNumall(:,time1:time2));
c = colorbar; c.Label.String = 'Zrms/d';
ylabel('s') ; xlabel('Ur');
set(h,'edgecolor','none','facecolor','interp');
hold on
title('Z space-time evolution, ONSAS')
% figure()
% contour(times(time1:time2),s.*(l/d),zdefNumall(:,time1:time2)) % Contour plot 
%}
function connectNodes(fig, connection, mesh)
    for i = 2:length(connection)
        node1 = connection(i, 1);
        node2 = connection(i, 2);
        x1 = mesh.nodesCoords(node1,3)
        y1 = mesh.nodesCoords(node1,2);
        x2 = mesh.nodesCoords(node2,3);
        y2 = mesh.nodesCoords(node2,2);
        fig
        plot([x1,x2], [y1,y2], 'm-o', 'LineWidth', 12);
        hold on
    end
end
%
function lk = Yfft(ydefNumlast,Ns, dt, plotBool, MinPeakProminence, Ur, E, I, ms, ma, l)
    xhat = fft(ydefNumlast,Ns);
    PSD = xhat.*conj(xhat)/Ns ;
    freq = 1/(dt*Ns)*(0:Ns);
    L = 1:floor(Ns/7); % Only plot the first half of frequencies
    [pkf,lkf] = findpeaks(PSD(L),freq(L),'MinPeakProminence',MinPeakProminence, 'MinPeakDistance', 2);%,'MinPeakProminence',3e-12
    if plotBool
        betavect = [1.87510 4.69409 7.857476]; freqStruct = (betavect.^2)*sqrt(E*I/(ms+ma))/(2*pi*l^2)  % In water  (beta^2 / L^2) sqrt(EI / rho A)
        figure(50) % plot FFT
        plot(freq(L), PSD(L))
        hold on 
        plot(lkf,pkf,'ro','MarkerSize',12)
        hold on         
        xline(freqStruct(1), '-.');xline(freqStruct(2), '-.');xline(freqStruct(3), '-.');
        title('FFT'); xlabel('f(Hz)')
        hold on
        xlim([0 40])
    end
    if Ur<=10
        [val, idx] = max(pkf);
        lk = lkf(idx);
    else
        lk = [];
        while max(pkf) > 5
            pkf
            [val, idx] = max(pkf);
            lk(end+1) = lkf(idx);
            pkf(idx) = -Inf;
        end
     end
end
% 8 shape
function [XYlist, gamma] = lemniscate(f, dt, xdefNumall, ydefNumall, plotBool, Ur,d)
if (f==0); Nsteps = 666; else; Nsteps = floor(1/f/dt)+100; end;
    xlemniscate = xdefNumall(end, end-Nsteps:end)./d;
    ylemniscate = ydefNumall(end, end-Nsteps:end)./d;
    Xmax = max(xlemniscate)-min(xlemniscate);
    Ymax = max(ylemniscate)-min(ylemniscate);
    gamma = Ymax/Xmax;
    if plotBool
        figure(30)
        plot(xlemniscate+Ur, ylemniscate)
        hold on
    end
    XYlist=[xlemniscate;ylemniscate]';
end
%
function plotspectro(Ur,indvel, fi, fig)
    if indvel == 1
        %openfig(fig);
        figure(51)
        ylabel('f'); xlabel('Ur'); title('Simulations Tip frequencies');
        plot(Ur, fi, 'k*','MarkerSize',3)
        hold on
    else
        figure(51)
        plot(Ur, fi, 'k*','MarkerSize',3)
        hold on
    end
end