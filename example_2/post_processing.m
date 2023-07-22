% Post processing signal Zdefnum
%zdefNumlast = SolZTF50s ; % = Z for last node
%clear all; close all;
global finalTime; global qvect
mattolatexfig = ['C:\Users\alvilh\matlab2tikz-master\matlab2tikz-master\src'] ;
addpath( genpath( mattolatexfig ) );
Y_rmsvect = zeros(numElements+1, NR);  Ym_rmsvect = zeros(1, NR); Urvect = zeros(1, NR); fvect = {};%zeros(1, NR);
Ns = floor(finalTime_Dynamic/dt);
% path =  'D:\Alexandre\session Poly\Experimental\Data\LEGH\Rod\0902mainppump4K120fps';
% load( path, filename));
%%
stop
A = 1:182; B = 186:10:NR; %A = 1:135; B = 96:10:NR; 
C = cat(2, A, B);
for indvel = [1:145 147:NR] %[88 89 125 130]
    vwindMax = velocities(indvel);
    %sprintf('Ysolutions\\Ysdrd_udotGfull_va=%.3f_Nelem=%d_FT%d.mat', vwindMax, numElements, finalTime_Dynamic)
    % Best for 1:43
    if ismember(indvel,[88 89 125 130])
        dt = 0.0005; finalTime_Dynamic = 5;
        load(sprintf('Ysolutions\\YflowBoolTrueV6_va=%.4f_Nelem=%d_FT%d_dt=%.4f.mat', vwindMax, numElements, finalTime_Dynamic ,dt));
    elseif ismember(indvel, [5 7 11 15 19]) % add 34,37,38, 41?
        dt = 0.02; finalTime_Dynamic = 5;
        load(sprintf('Ysolutions\\YflowBoolTrue_va=%.4f_Nelem=%d_FT20.mat',vwindMax,numElements));
    elseif indvel>=1 && indvel<=20
        dt = 0.02; finalTime_Dynamic = 5;
        load(sprintf('Ysolutions\\YflowBoolTrue_va=%.3f_Nelem=%d_FT20.mat',vwindMax,numElements)); finalTime_Dynamic=20; % (2*rand(numElements,1)-1)*0.001; % Best for 1:42
    elseif indvel>=45 && indvel<=95
        finalTime_Dynamic=0.8;dt = 0.0005;
        load(sprintf('Ysolutions\\YflowBoolTrueV4_va=%.3f_Nelem=%d_FT%d_dt=%.4f.mat',vwindMax, numElements, finalTime_Dynamic,dt)); % 45:95,
    elseif indvel>=95 && indvel<=135
        finalTime_Dynamic=0.5; dt = 0.0005;
        load(sprintf('Ysolutions\\YflowBoolTrueV4_va=%.3f_Nelem=%d_FT%d_dt=%.4f.mat',vwindMax, numElements, finalTime_Dynamic,dt)); % 95:135,
    % Best for 21:44
    elseif indvel>=21 && indvel<=44
        finalTime_Dynamic=2; dt = 0.0005;
        %sprintf('Ysolutions\\YflowBoolTrueV4_va=%.3f_Nelem=%d_FT%d_dt=%.4f.mat', vwindMax, numElements, finalTime_Dynamic,dt)
        %stop
        load(sprintf('Ysolutions\\YflowBoolTrueV4_va=%.3f_Nelem=%d_FT%d_dt=%.4f.mat', vwindMax, numElements, finalTime_Dynamic,dt));  % 2*(2*rand(numElements, 1)-1);
    elseif indvel>=136 && indvel<=142 || ~isempty(find(96:10:NR == indvel))
        finalTime_Dynamic = 5; dt = 0.0005;
        load(sprintf('Ysolutions\\YflowBoolTrueV4_va=%.3f_Nelem=%d_FT%d_dt=%.4f.mat', vwindMax, numElements, finalTime_Dynamic,dt));  % 2*(2*rand(numElements, 1)-1);
    elseif indvel>=143 && indvel<=NR %&& isempty(find(96:10:NR == indvel))
        finalTime_Dynamic = 2; dt = 0.0005;
        load(sprintf('Ysolutions\\YflowBoolTrueV5_va=%.3f_Nelem=%d_FT%d_dt=%.4f.mat', vwindMax, numElements, finalTime_Dynamic,dt));  % 2*(2*rand(numElements, 1)-1);
    end
    %Urvect(1,indvel) = (St*(l^2)*vwindMax/d)*sqrt((ms+ma)/(E*I));
    Urvect(1,indvel) = ((l^2)*vwindMax/d)*sqrt((ms+ma)/(E*I));
    Y_rmsvect(:,indvel) = rms( (ydefNumall(:, floor(end/2):end) - mean(ydefNumall(:, floor(end/2):end),2))./d,2 );
    ydefNumlast = ydefNumall(end, :)'; %xdefNumall = matUsDynLD2D(1:6:end, :) ;
    %RMSall = rms(ydefNumall(:, end/2:end)./d, 2); % size(RMSall) = 101 1
    %Ym_rmsvect(1,indvel) = max(RMSall);
    sig = (ydefNumlast(end/2:4*end/4) - mean(ydefNumlast(end/2:4*end/4)))./d;
%     figure(3)
%     plot(sig)
%     hold on
%     rms(sig)
    Ym_rmsvect(1,indvel) = rms(sig);
    %fw0i = (St*vwindMax/d); 
    lk = Yfft(ydefNumlast(end/2:4*end/4),Ns, dt, false, 1e-18, Urvect(1,indvel));
    fs1 = 1;
    if isempty(lk)
        indvel
        %fvect(1,indvel) = NaN; % Should not achieve
        fi = lk./fs1;
    else
        %fvect(1,indvel) = Yfft(ydefNumlast(end/2:end),Ns, dt, false, 1e-18, Urvect(1,indvel))./fw0i; 
        fi = lk./fs1;
    end
    plotspectro(Urvect(1,indvel), indvel, fi, 'PSD3.fig')
    % Spectrogram
    %[uvect Fvect]= spectroDots(ydefNumlast, 512, 512/2, (St*(l^2)*0.001/d)*sqrt((ms+ma)/(E*I)), -90, l, ms+ma, E, I)
end
times = linspace(0,finalTime_Dynamic, size(ydefNumall, 2));
%% Plot RMS
% figure(1)
% hold on 
%fig = openfig('Yplot2.fig');
figure()
plot(Urvect(1:end),Ym_rmsvect, 'r-o') % Y_rmsvect(end,1:end)
%xlim([0 Urvect(1, 120)])
%xlabel('St Ur'); ylabel('Yrms/d'); 
%legend('Leclercq & al', 'ONSAS')
%title('Tip RMS transverse displacement')
%% SPACE-TIME EVOLUTION
% amplitude lin plot
Yrmsallmaxtimes = max(Y_rmsvect(:,1:2:end),[],1); % Max for each time: vector length times
s = linspace(0,1,numElements+1);
openfig('Amplitude2.fig');
title(''); xlabel('');ylabel('');
subplot(2,1,2)
%figure()
%h=pcolor(Urvect,s,Y_rmsvect(:,:)); % Yrms not normalized
h=pcolor(Urvect(1:2:end),s,Y_rmsvect(:,1:2:end)./repmat(Yrmsallmaxtimes,length(s),1)); % Yrms normalized
c = colorbar; %c.Label.String = 'Yrms/d';
%ylabel('s') ; xlabel('Ur');title('Y space-time evolution, ONSAS')
set(h,'edgecolor','none','facecolor','interp');
%hold on
%title('Z space-time evolution, ONSAS')
%% Plot signal
figure()
plot(times, ydefNumlast(:)./d, 'r*')
rms(ydefNumall(end,end/2:end)./d)
rms(ydefNumall(end,end/4:2*end/4)./d)
fm = Yfft(ydefNumlast(end/2:end),Ns, dt, true, 1e-18, Urvect(1,indvel))
%% Plot frequencies
% Overplotting on Leclercq
%openfig('PSD2.fig');
%plot(Urvect, fvect, 'k*')
matlab2tikz('figure(1).tex');
%% Plot Ampmlitude and q along s for Ur constant
NdotsZ = 1600;
AmpAndQplotalongs(numElements, NdotsZ, zdefNumall, times, qvect, UrLeclercq, d)
%% Plot tip Amplitude and q of tip along t for Ur constant
figure(10)
NdotsZ = 200; UrLeclercq = 2;
AmpAndQplotalongt(numElements, NdotsZ, ydefNumall, times, qvect, UrLeclercq, d)
%% Structural Modes
% From Fundamentals of vibrations by Leonard Meirovotch
% For cantilever beams clamped-free end only
betavect = [1.87510 4.69409 7.857476 10.9955 14.1372]; % beta*l for cantilever clamped-free end beam!
%if l = 1 !: 
%freqStruct = (betavect.^2)*sqrt(E*I/(rho*pi*(d/2)^2))/(2*pi) ; % In air
%freqStruct = (betavect.^2)*sqrt(E*I/((rhoFluid+rho)*pi*(d/2)^2))/(2*pi) ; % In a fluid
freqStruct = (betavect.^2)*sqrt(E*I/(ms+ma))/(2*pi*l^2)
%% Ur list for abscisse
UR=zeros(Ns,1);
for i = 1:Ns
    va_vector = feval(nameFuncVel, [0, 0, 0], times(i));
    u         = norm(va_vector);
    UR(i)     = (St*(l^2)*u/d)*sqrt(mtot/(E*I));    
end
dur = 1/(UR(2)-UR(1)); % For fft
%title(sprintf('signal, Ur from 0 to %.2d', UrLeclercq)); xlabel('Ur');
figure(11)   
plot(UR, ydefNumlast)
%{
%% Compute Zrms at each node
windowlenghtmean = 400;
windowlenghtmax = 50;
nlast = Ns - mod(Ns, windowlenghtmax); % So length(zdefNumlast(1:nlast))%windowlenght = 1
Zrmsall = sqrt(movmean(zdefNumall.^ 2, windowlenghtmean, 2))/d;  % Sliding window on  columns of zdefNumall (on time) %Zrmsall(end, :) == movrms % OUI
%Zrmsallmax = max(Zrmsall,[],2); % Max for each node: vector length numNodes
Zrmsallmaxtimes = max(Zrmsall,[],1); % Max for each time: vector length times
Zrmsmax = max(reshape(Zrmsallmaxtimes(1:nlast), windowlenghtmax, [])); % Counting every data only once
openfig('Yplot2.fig');
%plot(UR, Zrmsallmaxtimes)
hold on
plot(UR(1:windowlenghtmax:end-windowlenghtmax), Zrmsmax)
hold on
plot(UR, Zrmsall(end,:))
legend('Leclercq and de Langre','ONSAS')
%}
%% Test Zrms low Ur
windowlenght = 300; %floor(Ns/10); %60
nlast = Ns - mod(Ns, windowlenght);
Zrms1 = sqrt(mean(reshape(ydefNumlast(1:nlast), windowlenght, []) .^ 2))/d; % Counting every data only once
%plot(UR(1:windowlenght:end-windowlenght), Zrms1)
figure(11)
plot(UR(1:windowlenght:end-windowlenght), Zrms1)
hold on
stop
%% Compute Zrms at each period
% Zrms = rms(Zvect(1, :)); % RMS = @(x) sqrt(mean(x.^2)); %RMS of amplitude
% Moving rms of Z displacement
windowlenghtmean1 = 300; %60
windowlenghtmean2 = 100;%60
windowlenghtmean3 = 100; %200
windowlenghtmax = 100; %300
nmaxlast = Ns - mod(Ns, windowlenghtmax);
movrms = sqrt(movmean(ydefNumlast.^ 2, windowlenghtmean1))/d;  % Sliding window on the first peak
% movrms(end/10:9*end/10) = sqrt(movmean(zdefNumlast(end/10:9*end/10).^ 2, windowlenghtmean2))/d;  % Sliding window
% movrms(9*end/10:end) = sqrt(movmean(zdefNumlast(9*end/10:end).^ 2, windowlenghtmean1))/d;  % Sliding window
% movrms(end/20:19*end/20) = sqrt(movmean(zdefNumlast(end/20:19*end/20).^ 2, windowlenghtmean2))/d;  % Sliding window
% movrms(4*end/7:end) = sqrt(movmean(zdefNumlast(4*end/7:end).^ 2, windowlenghtmean3))/d;  % Sliding window on the end
% Taking max 
Zrmsmax = max(reshape(movrms(1:nmaxlast), windowlenghtmax, [])); % Counting every data only once
%nlast = Ns - mod(Ns, windowlenght); % So length(zdefNumlast(1:nlast))%windowlenght = 1
% Zrms1 = sqrt(mean(reshape(zdefNumlast(1:nlast), windowlenght, []) .^ 2))/d; % Counting every data only once
%Envelop
%[up,lo] = envelope(Zrms,15,'peak');
%y = envelope(zdefNumlast,100,'rms')/d;
%{no
%openfig('Yplot2.fig');
%hold on
plot(UR(1:windowlenghtmax:end-windowlenghtmax), Zrmsmax)
hold on
%plot(UR, movrms)
legend('Leclercq and de Langre','ONSAS')
%{
figure(12)
title('Dimensionless RMS Amplitude of vibration'); xlabel('Ur');ylabel('Zrms/d') ;
%plot(UR(100:end), movrms(100:end)) % Sliding window 
plot(UR(1:windowlenghtmax:end-windowlenghtmax), Zrmsmax)
% Plot dots from the 50 simulations
% hold on 
%plot(URvect, Zrmsvect, 'ro')
hold on;
%plot(xldl, yldl)
plot(xLDL, yLDL)
%}
%% Mesh analysis
for n = [15 30 50 100]
    load(sprintf('Zwinstep_va=0.04_Nelem=%d_FT50_Nstep=20_coefacc=0.2.mat', n));
    ydefNumlast = zdefNumall(end,:);
    Ns = length(ydefNumlast);
    windowlenghtmean1 = 300;
    movrms = sqrt(movmean(ydefNumlast.^ 2, windowlenghtmean1))/d;
    Zrmsmax = max(reshape(movrms(1:nmaxlast), windowlenghtmax, []));
    figure(21)
    plot(UR(1:windowlenghtmax:end-windowlenghtmax), Zrmsmax)
    hold on
end
legend('Nelem = 15','Nelem = 30', 'Nelem = 50', 'Nelem = 100')
%%  FFT and signal     
% N = 1; % Better period if transient not cut
% FT = fft(zdefNumlast(N:end));
%from steve brunton (weird it offsets on freq axis)
xhat = fft(ydefNumlast,Ns);
PSD = xhat.*conj(xhat)/Ns ;
freq = 1/(dt*Ns)*(0:Ns);
L = 1:floor(Ns/7); % Only plot the first half of frequencies
%{ no
figure(10) % plot FFT
plot(freq(L), PSD(L))
title('FFT'); xlabel('f(Hz)')
for i=1:4
    xline(freqStruct(i),'-.',sprintf('Structural Mode %d', i));
end
hold on
%}
%% Plot flow
Nstep = 20; coefacceleration = 0.2; 

lflat = floor((1-coefacceleration)*finalTime/Nstep/dt); %length(times(floor(coefacceleration*finalTime/Nstep/dt)));
flatInt = zeros(Nstep, lflat);
for j = 1:Nstep % indices on the jth flat part of the step flow function
  flatInt(j,:) = floor(j*finalTime/Nstep/dt)-lflat+1: floor(j*finalTime/Nstep/dt);
end
figure()
plot(flatInt, 'k-o')
%% Spectrogram
dur = 1/(Urvect(2)-Urvect(1));
Nfft = 2^nextpow2(Ns);
Nspec = 512;
Noverlap = Nspec/2;
%fintv = 0.1:0.001:0;
[~,F,T,P] = spectrogram(ydefNumlast, Nspec, Noverlap, Nspec, dur,'Yaxis');
%[~,~,~,pxx2,fc2,uc2] = spectrogram(zdefNumlast, Nspec, Noverlap, Nspec, dur,'Yaxis');
%{no
figure() % Plot spectrogram
%imagesc(T(1:138), F(1:138)./pi, 10*log10(P(1:138,1:138)+eps)); % add eps like pspectrogram does and nondimensionalize by pi
imagesc(T(1:end), F(1:size(T))./pi, 10*log10(P(1:length(T),1:end)+eps)); 
fw0vect = zeros(length(F), 1);
for i = 1:length(length(F))
    va_vector = feval(nameFuncVel, [0, 0, 0], times(i));
    u         = norm(va_vector);
    fw0vect(i) = St*u/d; 
end
F2 = F./fw0vect;
%F2 = F./fw0vect.' ;
% C = 1:length(F); % only changes the axis not the values
%imagesc(T(1:end), F(1:end)./pi, 10*log10(P(1:end,1:end)+eps)); % add eps like pspectrogram does and nondimensionalize by pi
%imagesc(uc2(pxx2>0),fc2(pxx2>0)./pi, 10*log10(pxx2(1:138,1:138)+eps)); % add eps like pspectrogram does 
%axis xy
h = colorbar;
h.Label.String = 'Power/frequency (dB/Hz)';
ylabel('PSD') ;
%axis([T(1) T(end) 0 32])
%axis([0 UR(end) 0 35])
%caxis([-150,-80]);
title('Spectrogram'); ylabel('f (Hz)'); xlabel('Ur');
for i=1:4
    yline(freqStruct(i),'-.',sprintf('Structural Mode %d', i));
    hold on
end
%}
%% Spectrogram dimensionless dots
%UR(i) = (St*(l^2)*u/d)*sqrt((ms+ma)/(E*I));
%{
figure(13)
plot(uc(pxx>0),Fvect,'.')
hold on;
yline(1,'-.');
axis([0 35 0 1.5])
xlabel('Ur');ylabel('f/fw0'); 
%}
%% Shape of the beam
%{no
lockInstimes = [4.7 12.6 33.3 68 120];
lockInsN = floor(lockInstimes./dt); % 1678        4500       11892       24285
print_mode = 5;
figure()
%subplot(2,1,1)
plot(zdefNumall(:,51000)./d,s,'linewidth',2)
%axis equal
grid on
%axis([-0.1 1.1 -0.1 1.1])
axis([-1 1 -0.1 1.1])
xlabel('z/d')
ylabel('s')
title('shape of the beam')
%}
%% Plot Peaks on spectrogram
%Need to run FFT in spectrogramlooptest
% figure(5)
% hold on
% for col = 1:length(peaksvect(1,1:end)) %50
%   for line = 1:length(peaksvect(1:end,col)) % 4 if 4 peaks detected once
%     if peaksvect(line,col)>0
%       plot(URvect(col), peaksvect(line,col), 'ro') 
%       hold on
%     end
%   end
% end
%% Spectro function
function [uvect Fvect]= spectroDots(ydefNumlast, Nspec, Noverlap, dur, MinThreshold, l, mtot, E, I)
[~,~,~,pxx,fc,uc] = spectrogram(ydefNumlast, Nspec, Noverlap, Nspec, dur,'MinThreshold',MinThreshold); 
%{
figure(11) % Plot spectrogram dots and movrms
plot(uc(pxx>0),fc(pxx>0)./pi,'.') 
yyaxis left; ylabel('f (Hz)') ; xlabel('Ur');
for i=1:4
    yl = yline(freqStruct(i),'-.',sprintf('Mode %d', i));
    yl.LabelHorizontalAlignment = 'left';
    hold on
end
hold on
yyaxis right; ylabel('Zrms/d') ;
plot(UR(100:end), movrms(100:end), 'r')
axis([0 UR(end) 0 0.5])
%}
% Should non dimensionalize the frequencies by natural shedding frequency
% fw0 = St U / d
Fvect = fc(pxx>0)./pi; %careful spectogram function in Matlab normalizes by pi frquencies
Uri = uc(pxx>0);
for i = 1:length(fc(pxx>0))
    fw0i = Uri(i)/((l^2)*sqrt(mtot/(E*I))); % without added mass
    Fvect(i) = Fvect(i)/fw0i; 
end
uvect = uc(pxx>0);
% Overplotting on Leclercq
openfig('PSD3.fig');
plot(uc(pxx>0),Fvect,'.')
hold on
end
% Watch wind
% wind=feval(nameFuncVel, [0, 0, 0], times);
% figure(4)
% plot(times, wind(1:end-2))
function velvect = computestepVel(Ns, Nstep, finalTime, coefacceleration)
    velvect = zeros(Ns);
    for i = 1:100:Ns
      t= times(i);
      step = floor(Nstep*t/finalTime);
      velvect(i) =  (vwindMax/ (coefacceleration*finalTime) *(t-finalTime*step/Nstep) +coefacceleration*finalTime * vwindMax * (step)/(Nstep *coefacceleration* finalTime))*(t <= finalTime/Nstep *(step+ coefacceleration)) ...
                  +coefacceleration*finalTime * vwindMax * (step+1)/(Nstep *coefacceleration* finalTime)*(t > finalTime/Nstep *((step+ coefacceleration)));
    end
    % figure(1)
    % plot(times,velvect, 'r*')
end    
%
function AmpAndQplotalongs(numElements, NdotsY, zdefNumall, times, qvect, Ur, d)
    snodes = linspace(0,1,numElements+1);
    c = colorbar; c.Label.String = 'Zrms/d';
    subplot(2,1,1)
    Nspace = 100;
    %Z = zdefNumall(:,end-200:end)./repmat(Zrmsallmaxtimes,length(snodes),1)
    Z = zdefNumall(:,end-NdotsY:end);
    h = pcolor(times(end-NdotsY:end),snodes,Z); % Zrms normalized
    c = colorbar; c.Label.String = 'Zrms/d';
    ylabel('s') ; xlabel('Times (s)');
    set(h,'edgecolor','none','facecolor','interp');
    hold on
    title(sprintf('Z space-time evolution, ONSAS, Ur = %.2f', Ur))
    %h=pcolor(times(end-200:end),s,Y(:,end-200:end));
    colorbar
    set(h,'edgecolor','none','facecolor','interp');
    subplot(2,1,2)
    selem = linspace(0,1,numElements);
    size(times)
    size(selem)
    size(qvect(1:2:end,:))
    h=pcolor(times(1:end-1),selem ,qvect(1:2:end,:));
    colorbar
    set(h,'edgecolor','none','facecolor','interp');
    xlabel('t');ylabel('s');
    title('q')
end
%
function AmpAndQplotalongt(numElements, NdotsZ, zdefNumall, times, qvect, Ur, d)
    qtip = qvect(end-1,:);
    Ztip = zdefNumall(end, :);
    subplot(2,1,1)
    plot(times(end-NdotsZ:end), Ztip(end-NdotsZ:end)./d)
    xlabel('times');ylabel('Z/d');
    title(sprintf('Tip Amplitude, ONSAS, Ur = %.2f', Ur))
    subplot(2,1,2)
    size(times)
    size(qtip)
    plot(times(end-NdotsZ:end), qtip(end-NdotsZ:end))
    xlabel('times');ylabel('q');
end
%
function lk = Yfft(ydefNumlast,Ns, dt, plotBool, MinPeakProminence, Ur)
    xhat = fft(ydefNumlast,Ns);
    PSD = xhat.*conj(xhat)/Ns ;
    freq = 1/(dt*Ns)*(0:Ns);
    L = 1:floor(Ns/7); % Only plot the first half of frequencies
    [pkf,lkf] = findpeaks(PSD(L),freq(L),'MinPeakProminence',MinPeakProminence, 'MinPeakDistance', 2);%,'MinPeakProminence',3e-12
    if plotBool
        figure(50) % plot FFT
        plot(freq(L), PSD(L))
        hold on 
        plot(lkf,pkf,'ro','MarkerSize',12)
        title('FFT'); xlabel('f(Hz)')
        hold on
        xlim([0 40])
    end
    if Ur<=149
        [val, idx] = max(pkf);
        lk = lkf(idx);
    else
        lk = [];
        while max(pkf) > 1e-7
            pkf
            [val, idx] = max(pkf);
            lk(end+1) = lkf(idx);
            pkf(idx) = -Inf;
        end
     end
end
%
function plotspectro(Ur,indvel, fi, fig)
    if indvel == 1
        openfig(fig);
        %title('Tip frequency Power Spectral Density')
        %ylabel('f')
        h.Label.String = '';
        plot(Ur, fi, 'k*','MarkerSize',3)
        hold on
    else
        %figure(9)
        plot(Ur, fi, 'k*','MarkerSize',3)
        hold on
    end
end