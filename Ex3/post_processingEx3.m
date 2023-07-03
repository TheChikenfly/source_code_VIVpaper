% Post processing signal ydefnum
%ymid = SolZTF50s ; % = Z for last node
global finalTime;
%mattolatexfig = ['C:\Users\alvilh\matlab2tikz-master\matlab2tikz-master\src'] ;
%addpath( genpath( mattolatexfig ) );
nelem = 1000;
loadParametersEx3
for ps = [T0/(E*pi*d^2/4)] % 5e-4 6.1e-4 T0/(E*pi*(d^2-dint^2)/4) T0/(E*pi*d^2/4)
% drag in norm( VpiRelGflow) * VpiRelGflow Too large !
%finalTime = 5;
%load(sprintf('YEx3Nelem=%d_FT=%.0d_vmax=%.1f_pretension=%.0d.mat', numElements, finalTime, vwindMax,ps));
% drag in norm( VpiRelG) * VpiRelG
finalTime = 5; % 
%load(sprintf('YEx3Nelemv2=%d_FT=10_vmax=%.1f_pretension=%.1d.mat', numElements, vwindMax,ps)); % 3*end/10:4*end/10
% flowBool false constdirlift false   GOOD RESULTS here modes 3 and 5 or 2 and 5
%load(sprintf('YEx3NelemV3=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps));
% flowBool false constdirlift false drag in norm( VpiRelG) VpiRelG
%load(sprintf('YEx3NelemV4=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% With amplification factor drag norm( VpiRelG) * VpiRelG
%load(sprintf('YEx3NelemV5=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% Try lift in norm(VpiRelG)^2, drag in norm( VpiRelG) * VpiRelG and amp factor GOOD RESULTS here, modes 3 and 5
%load(sprintf('YEx3NelemV6=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% Begin with static deformed matrix, modes 3 and 7 & 4 and 2.5
%load(sprintf('YEx3NelemV7=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% constliftdir true, fll =  ...* norm( VpiRelG )^2 * tlift_defCoords c_d_il = 0.2; BEST RESULTS ! ps = T0/(E*pi*(d^2)/4)
%load(sprintf('YEx3NelemV8=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% Try constliftdir false lift in norm( VpiRelG ) * VpiRelGperp , Ax = 96; epsilonx = 0.02, no famp, c_d_il = 0.1! St0.2
%load(sprintf('YEx3NelemV9cdlFalseAx96=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% cu = 0.03
%load(sprintf('YEx3NelemV9bis=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% Try same (qvect = 2*) but St = 0.17; Cd0= 2 MODE 2
% load(sprintf('YEx3NelemV10=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% qvect(1:2:end,1) = 4*..., St = 0.17; Cd0= 1.2
%load(sprintf('YEx3NelemV11=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% qvect(1:2:end,1) = 0.1*..., St = 0.17; Cd0= 1.2 MODE 2??
% finalTime = 20;
% load(sprintf('YEx3NelemV12=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% cu =0, qvect(1:2:end,1) = 2*..., St = 0.17; MODE 2!
%load(sprintf('YEx3NelemV13=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% cu =0.03, qvect(1:2:end,1) = 1*..., St = 0.2; ?
%load(sprintf('YEx3NelemV14=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% cu = 0 
%load(sprintf('YEx3NelemV14bis=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% cu =0, qvect(1:2:end,1) = 4*..., St = 0.2; BEST one maybe
%load(sprintf('YEx3NelemV15=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% 1000 elem
%load(sprintf('YEx3NelemStudy=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', nelem, finalTime, vwindMax,ps))
%  Ay = 12; epsilony = 0.04!!,St = 0.2 cu=0 results too high BEST Results with FT=5
%finalTime=20;
%load(sprintf('YEx3NelemV16=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% n=100, qvect(1:2:end,1) = 0.001*... weird
%load(sprintf('YEx3NelemV16=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% n=1000, qvect(1:2:end,1) = 0.001*... BEST
%load(sprintf('YEx3NelemV16.2=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% n=1000, qvect(1:2:end,1) = 0.001*..., cu =  0.0002
%load(sprintf('YEx3NelemV16.3=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% n=1000, qvect(1:2:end,1) = 0.001, cu =  0
%load(sprintf('YEx3NelemV16.3bis=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% fixed tlift in womV4
load(sprintf('YEx3NelemV16.4=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
% changed final time to have 10 periods
finalTime = 7.14;
load(sprintf('YEx3NelemV16.5=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
%  Ay = 12; epsilony = 0.04, St=0.17 cu=0 results too high
%load(sprintf('YEx3NelemV17=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
%  epsilony = 0.04, St=0.2 cu=0.003 MODE 4??
%load(sprintf('YEx3NelemV18=%d_FT=%.0d_vmax=%.1f_pretension=%.1d.mat', numElements, finalTime, vwindMax,ps))
ymid = ydefNumall(numElements/2, :); xmid = xdefNumall(numElements/2, :);
s = linspace(0,1,numNodes); 
Ns = length(ymid);
times = linspace(0,finalTime, Ns);
%% Plot signal
figure()
plot(times, ydefNumall(50, :)./d, '*')
%plot(times, xdefNumall(50, :)./d, '*')
% sig = (ydefNumall(10, :) - mean(ydefNumall(10, :), 2))./d;
% hold on 
% plot(times, sig)
%% Trim & al. & Holmes & al. solution
xtrimCF2 = [0	0.094308943	0.214634146	0.32195122	0.393495935	0.429268293	0.435772358	0.416260163	0.380487805	0.33495935	0.318699187	0.318699187	0.315447154	0.32195122	0.331707317	0.354471545	0.383739837	0.419512195	0.448780488	0.461788618	0.448780488	0.416260163	0.390243902	0.357723577	0.331707317	0.331707317	0.357723577	0.41300813	0.461788618	0.520325203	0.575609756	0.624390244	0.630894309	0.617886179	0.598373984	0.565853659	0.507317073	0.455284553	0.393495935	0.318699187	0.234146341	0.162601626	0.104065041	0];
ytrimCF2 = [0	0.018779343	0.042253521	0.072769953	0.103286385	0.129107981	0.159624413	0.187793427	0.215962441	0.262910798	0.281690141	0.295774648	0.305164319	0.323943662	0.356807512	0.392018779	0.422535211	0.448356808	0.478873239	0.511737089	0.539906103	0.582159624	0.617370892	0.661971831	0.692488263	0.708920188	0.730046948	0.758215962	0.774647887	0.79342723	0.816901408	0.842723005	0.85915493	0.88028169	0.896713615	0.910798122	0.927230047	0.936619718	0.948356808	0.962441315	0.971830986	0.983568075	0.995305164	1.004694836];
xtrimIL2 = [0	0.042276423	0.084552846	0.123577236	0.146341463	0.159349593	0.149593496	0.130081301	0.104065041	0.087804878	0.07804878	0.07804878	0.094308943	0.110569106	0.130081301	0.136585366	0.136585366	0.123577236	0.100813008	0.094308943	0.091056911	0.087804878	0.094308943	0.104065041	0.123577236	0.130081301	0.123577236	0.110569106	0.097560976	0.110569106	0.120325203	0.136585366	0.143089431	0.143089431	0.136585366	0.120325203	0.104065041	0.100813008	0.104065041	0.123577236	0.139837398	0.169105691	0.191869919	0.208130081	0.208130081	0.195121951	0.175609756	0.152845528	0.123577236	0.091056911	0.058536585	0.032520325	0];
ytrimIL2 = [0	0.016431925	0.03286385	0.056338028	0.079812207	0.093896714	0.11971831	0.143192488	0.17370892	0.194835681	0.208920188	0.225352113	0.241784038	0.262910798	0.281690141	0.298122066	0.316901408	0.335680751	0.356807512	0.368544601	0.389671362	0.401408451	0.417840376	0.436619718	0.464788732	0.495305164	0.521126761	0.544600939	0.563380282	0.58685446	0.605633803	0.629107981	0.647887324	0.67370892	0.694835681	0.713615023	0.727699531	0.751173709	0.767605634	0.781690141	0.795774648	0.816901408	0.833333333	0.852112676	0.877934272	0.896713615	0.910798122	0.929577465	0.948356808	0.962441315	0.976525822	0.988262911	1.004694836];
xCFDCF = [0	0.045528455	0.097560976	0.182113821	0.260162602	0.308943089	0.38699187	0.416260163	0.422764228	0.403252033	0.357723577	0.263414634	0.201626016	0.143089431	0.107317073	0.123577236	0.156097561	0.185365854	0.234146341	0.302439024	0.367479675	0.41300813	0.445528455	0.455284553	0.448780488	0.419512195	0.380487805	0.338211382	0.273170732	0.221138211	0.156097561	0.123577236	0.107317073	0.126829268	0.162601626	0.221138211	0.250406504	0.318699187	0.377235772	0.422764228	0.445528455	0.455284553	0.455284553	0.432520325	0.34796748	0.240650407	0.091056911	0 ];
yCFDCF = [0	0.011737089	0.023474178	0.044600939	0.0657277	0.079812207	0.112676056	0.138497653	0.154929577	0.192488263	0.220657277	0.26056338	0.281690141	0.302816901	0.321596244	0.342723005	0.356807512	0.366197183	0.384976526	0.408450704	0.4342723	0.453051643	0.481220657	0.5	0.521126761	0.544600939	0.563380282	0.582159624	0.600938967	0.61971831	0.643192488	0.659624413	0.67370892	0.690140845	0.708920188	0.725352113	0.737089202	0.758215962	0.786384977	0.812206573	0.830985915	0.842723005	0.863849765	0.887323944	0.922535211	0.948356808	0.985915493	1.002347418];
xCFDIL = [0.003252033	0.045528455	0.091056911	0.117073171	0.143089431	0.136585366	0.117073171	0.097560976	0.094308943	0.094308943	0.094308943	0.120325203	0.136585366	0.149593496	0.149593496	0.123577236	0.104065041	0.091056911	0.113821138	0.130081301	0.130081301	0.110569106	0.110569106	0.087804878	0.087804878	0.110569106	0.136585366	0.139837398	0.123577236	0.104065041	0.091056911	0.091056911	0.123577236	0.136585366	0.130081301	0.100813008	0.084552846	0.110569106	0.136585366	0.149593496	0.149593496	0.130081301	0.117073171	0.097560976	0.091056911	0.104065041	0.123577236	0.152845528	0.156097561	0.149593496	0.130081301	0.097560976	0.058536585	0.019512195	-0.003252033];
yCFDIL = [0	0.018779343	0.03286385	0.046948357	0.070422535	0.093896714	0.115023474	0.131455399	0.14084507	0.157276995	0.157276995	0.178403756	0.192488263	0.211267606	0.232394366	0.258215962	0.276995305	0.295774648	0.323943662	0.34741784	0.370892019	0.394366197	0.394366197	0.422535211	0.436619718	0.460093897	0.481220657	0.509389671	0.53286385	0.549295775	0.568075117	0.582159624	0.615023474	0.64084507	0.661971831	0.690140845	0.711267606	0.737089202	0.76056338	0.779342723	0.791079812	0.812206573	0.823943662	0.842723005	0.854460094	0.873239437	0.887323944	0.910798122	0.920187793	0.941314554	0.955399061	0.971830986	0.985915493	0.997652582	1.004694836];
%%  Yrms along s
signaly = ydefNumall(:, 5*end/10:10*end/10);
Yrms = rms(signaly./d, 2);
subplot(1,2,1)
plot(Yrms, s, 'k-'); 
hold on 
plot(xtrimCF2, ytrimCF2, 'b--')
hold on 
plot(xCFDCF, yCFDCF, 'm--')
title('Cross-flow RMS displacements'); ylabel('z/l');xlabel('Yrms/d') ;
legend('Present wake-oscillator', 'Experimental data (Trim & al)', 'CFD data (Holmes & al)'); 
relativeErrorY = (abs(max(Yrms)-max(xtrimCF2)))/max(xtrimCF2)
%% In line rms 
signalx = xdefNumall(:, 9*end/10:10*end/10)./d;
sig = signalx - mean(signalx, 2); 
%sig = (signalx - mean(signalx, 2));
Xrms = rms(sig, 2);
subplot(1,2,2)
plot(Xrms, s, 'k-'); 
hold on 
plot(xtrimIL2, ytrimIL2, 'b--')
hold on 
plot(xCFDIL, yCFDIL, 'm--')
title('In-line RMS displacements'); ylabel('z/l');xlabel('Xrms/d') ;
hold on 
legend('Present wake-oscillator', 'Experimental data (Trim & al)', 'CFD data (Holmes & al)'); 
%relativeErrorX = (abs(max(xtrimIL2)-max(Yrms)))/max(Yrms)
relativeErrorX = (abs(max(Xrms)-max(xtrimIL2)))/max(xtrimIL2)
end
%matlab2tikz('myfile.tex');
%% Structural Modes 
% From Fundamentals of vibrations by Leonard Meirovotch
% For cantilever beams pinned-pinned: 
betavect = [pi 2*pi 3*pi 4*pi]; % beta*l for cantilever pinned pinned beam!
%freqStruct = (betavect.^2)*sqrt(E*I/(ms/l))/(2*pi);  % In air 
freqStruct = (betavect.^2)*sqrt(E*I/(ms+ma)/l^4)/(2*pi);  % In water  (beta^2 / L^2) sqrt(EI / rho A)
%f1 = sqrt((E*I)/(ms+mf))/(St*(l^2)); 
%%  FFT and signal     
%{no
% FT = fft(ymid(N:end));
%from steve brunton (weird it offsets on freq axis)
xhat = fft(ymid(1:end),Ns); %same as xhat = fft(ymid(1:end));
PSD = xhat.*conj(xhat)/Ns ;
freq = 1/(dt*Ns)*(0:Ns);
L = 1:floor(Ns/20); % Only plot the first half of frequencies
figure(10)
plot(freq(L), PSD(L))
title('FFT'); xlabel('f(Hy)')
for i=1:4
    xline(freqStruct(i),'-.',sprintf('Structural Mode %d', i));
end
hold on
xline(f1,'-.','f1 in fluid');
% %% Space time evolution
% windowlenght = 50; %floor(Ns/10); 
% nlast = Ns - mod(Ns, windowlenght);
% yrmsw = sqrt(mean(reshape(ydefNumall(midNode,1:nlast), windowlenght, []) .^ 2))/d; % Counting every data only once
% movrms = sqrt(movmean(ydefNumall(midNode,:).^ 2, windowlenght))/d;  % Sliding window on the first peak
% siye(movrms)
% figure()
% plot(times, movrms)
% 
% figure()
% c = colorbar; c.Label.String = 'yrms/d';
% Nspace = 100;
% %h=pcolor(times,s,yrms./repmat(yrmsallmaxtimes,length(s),1)); % yrms normaliyed
% h=pcolor(times,s,yrms);
% c = colorbar; c.Label.String = 'yrms/d';
% ylabel('s') ; xlabel('Ur');
% set(h,'edgecolor','none','facecolor','interp');
% hold on
% title('y space-time evolution, ONSAS')
%}
% %% Spectrogram 
% Nfft = 2^nextpow2(Ns);
% Nspec = 512;
% Noverlap = Nspec/2;
% % dur = 1/(UR(2)-UR(1)); % For fft
% [~,~,~,pxx,fc,uc] = spectrogram(ymid, Nspec, Noverlap, Nspec, 1/dt,'MinThreshold',-135); 
% Fvect = fc(pxx>0)./pi; %careful spectogram function in Matlab normalizes frquencies by pi 
% Uri = uc(pxx>0);
% for i = 1:length(fc(pxx>0))
%     fw0i = Uri(i)/((l^2)*sqrt(mtot/(E*I))); % without added mass
%     Fvect(i) = Fvect(i)/fw0i; 
% end
% figure()
% imagesc(uc(pxx>0),fc(pxx>0)./pi, 10*log10(pxx(1:end,1:end)+eps)); 
% h = colorbar;
% h.Label.String = 'Power/frequency (dB/Hz)';
% ylabel('PSD') ;