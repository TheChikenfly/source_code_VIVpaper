% Post processing signal Ydefnum
addpath(genpath('../../../ONSAS.m/src')); % to have 'myCell2Mat'.
connec = myCell2Mat( mesh.conecCell ); 
connection = connec(:, 4:5);
loadParametersCoralEx5
tb1 = 76; tb2 = 101; fork = 51; %Tips at 76 and 101 and 150 fork = 51
j0 = 19;
%%
Y_rmsvectExp = zeros(numElements+1, NR); UrExp = zeros(1, NR); ftips = zeros(3, NR); %gamma = zeros(1, NR);
Y_rmsBranches = zeros(6, NR);
%velocities = 0.06:0.02:1.4;
%velocities = [0.01:0.01:0.01+8*0.01 0.12:0.02:1.4]; %NR=74
NR = length(velocities);
%velocities = 0.06:0.01:1.4;
l = 0.15;%l = 0.14;
FY_tip=NaN(1877,NR); PY_tip=NaN(1877,NR);
FY_b1=NaN(1877,NR); PY_b1=NaN(1877,NR); FY_b2=NaN(1877,NR); PY_b2=NaN(1877,NR); % CF
FY_b1il=NaN(1877,NR); PY_b1il=NaN(1877,NR); % IL
%
for j = 1:NR %[5 16 25 33 45 60], [5 16 21 25 33 45 60 68] (Ur = [1.8 6.5 8.5 13 19.5 29.3 41.5 48]) to do: 
    j
   U = velocities(j); UrExp(1, j) = ((l^2)*U/d)*sqrt((ms+ma)/(E*I)); % U = Ur*sqrt((E*I)/(ms+ma))*d/l^2;
   % lb in Ur ??
   if U < 0.215 %j<=j0 %
       dt = 0.001; finalTime = 5; Npsd = 1252;%finalTime/dt/2;, dt = 0.02; Npsd = 64;
       load(sprintf('Ysolutions\\YCoralN=1_va=%.3f_Nelem=%d_FT%d_dt=%.4f.mat',U, numElements, finalTime, dt));  
   elseif U < 1.03 %(j<=60)
       dt = 0.001; finalTime = 2; Npsd = 502;
       load(sprintf('Ysolutions\\YCoralN=1_va=%.3f_Nelem=%d_FT%d_dt=%.4f.mat',U, numElements, finalTime, dt));
   else %(j>60)
       dt = 0.0002; finalTime = 1.5; Npsd = 1877;
       load(sprintf('Ysolutions\\YCoralN=1_va=%.3f_Nelem=%d_FT%d_dt=%.4f.mat',U, numElements, finalTime, dt));
   end
   if ismember(j, [5 16 25 33 45 60]) %[5 16 21 25 33 45 49 68]
       rose = [1 0.5 0.8]; violet = [0.494 0.184 0.556];
       params{5} = {200 1 17 50}; %{200 3 11}; %tiini, step, nsteps, Nval 
       params{16} = {210 1 8 50};params{21} = {1810 5 13 300};params{25} = {1810 5 13 300};params{33} = {1925 4 10 300};params{45} = {1929 4 12 300};params{49} = {1930 5 10 300}; params{60} = {1930 5 10 150};
       %plotEnvelope(params{j}{1}, params{j}{2}, params{j}{3}, figure(j),connection, mesh,xdefNumall,ydefNumall,zdefNumall,false, false)
       %xz
       %plotEnvelope(params{j}{1}, params{j}{2}, params{j}{3}, figure(j+1),connection, mesh,ydefNumall,xdefNumall,zdefNumall,false, true)
       %3D
       %plotEnvelope(params{j}{1}, params{j}{2}, params{j}{3}, figure(j+2),connection, mesh,xdefNumall,ydefNumall,zdefNumall,true, false)
       %plotsigend(finalTime, dt, params{j}{4}, sig, sigb1, sigb2,xdefNumall,j , tb1, tb2)
   end
   %plotsigend(finalTime, dt, 200, sig, sigb1, sigb2,xdefNumall,j , tb1, tb2)
   % Remove static displacement to compute oscillations around steady deform
   xdefNumall = xdefNumall(:,:)-xdefNumall(:,1);ydefNumall = ydefNumall(:,:)-ydefNumall(:,1);zdefNumall = zdefNumall(:,:)-zdefNumall(:,1);
   sig = (ydefNumall(end, floor(end/2):end))./d; % end = tip of the trunk
   Y_rmsvectExp(:,j) = rms( (ydefNumall(:, floor(end/2):end))./d,2); % Computes at all nodes
   sigb1 = cos(angle)*ydefNumall(tb1, floor(end/2):end)./d - sin(angle)*zdefNumall(tb1, floor(end/2):end)./d;
   sigb2 = cos(angle)*ydefNumall(tb2, floor(end/2):end)./d + sin(angle)*zdefNumall(tb2, floor(end/2):end)./d;
   sigfork = ydefNumall(fork, floor(end/2):end)./d;
   Y_rmsBranches(1,j) = rms(sigb1); Y_rmsBranches(3,j) = rms(sigb2); Y_rmsBranches(5,j) = rms(sigfork);    
   Y_rmsBranches(2,j) = rms((xdefNumall(tb1, floor(3*end/4):end)- mean(xdefNumall(tb1, floor(3*end/4):end)))./d); 
   Y_rmsBranches(4,j) = rms((xdefNumall(tb2, floor(3*end/4):end)- mean(xdefNumall(tb2, floor(3*end/4):end)))./d); 
   Y_rmsBranches(6,j) = rms((xdefNumall(end, floor(3*end/4):end)- mean(xdefNumall(end, floor(3*end/4):end)))./d);    
   %ftips(1,j) = Yfft(sig,length(sig), dt, false, 0.1, UrExp(1, j), Etpe, I, ms, ma, l); %
   %ftips(2,j) = Yfft(sigb1,length(sigb1), dt, false, 0.1, UrExp(1, j), Etpe, I, ms, ma, l); %
   %ftips(3,j) = Yfft(sigb2,length(sigb2), dt, false, 0.1, UrExp(1, j), Etpe, I, ms, ma, l); %
   %l = 0.15;
%    fi = Yfft(sig,length(sig), dt, true, 0.1, UrExp(1,j), Etpe, I, ms, ma, l);
%    plotspectro(UrExp(1,j), j, fi, '')
    fw0i = 1;%(St*U/d);
    [fY,pY] = psdlec(sig,1,1/dt);
    FY_tip(1:Npsd,j)=fY.'./fw0i; PY_tip(1:Npsd,j)=pY.';
    [fYb1,pYb1] = psdlec(sigb1,1,1/dt);[fYb2,pYb2] = psdlec(sigb2,1,1/dt);[fYb1il,pYb1il] = psdlec(xdefNumall(tb1, floor(end/2):end)./d,1,1/dt);
    FY_b1(1:Npsd,j)=fYb1.'./fw0i; PY_b1(1:Npsd,j)=pYb1.';FY_b2(1:Npsd,j)=fYb2.'./fw0i; PY_b2(1:Npsd,j)=pYb2.';
    FY_b1il(1:Npsd,j)=fYb1il.'./fw0i; PY_b1il(1:Npsd,j)=pYb1il.';
end
%% plots
%NR = 68
figure(12)
%Y_rmsvectExp(Y_rmsvectExp==0) = []; UrExp(UrExp==0) = [];
%fig = openfig('ExpRms.fig'); hold on
plot(UrExp(1:NR), Y_rmsvectExp(end,1:NR), 'k-')
hold on 
plot(UrExp(1:NR), Y_rmsBranches(1,1:NR), 'r-')
hold on 
plot(UrExp(1:NR), Y_rmsBranches(3,1:NR), 'b-')
%hold on 
%plot(UrExp(1:NR), Y_rmsBranches(5,1:NR), 'g-')
hold on 
plot(UrExp(1:NR), Y_rmsBranches(2,1:NR), 'r*') 
hold on 
plot(UrExp(1:NR), Y_rmsBranches(4,1:NR), 'b*')
hold on 
plot(UrExp(1:NR), Y_rmsBranches(6,1:NR), 'k*') % IL trunk
xlabel('Ur =  U/(d f1)'); ylabel('Yrms/d'); title('Yrms amplitude');
legend('Experimental','ONSAS')
%
%IL
betavect = [1.87510 4.69409 7.857476]; freqStruct = (betavect.^2)*sqrt(E*I/(ms+ma))/(2*pi*l^2); % In water  (beta^2 / L^2) sqrt(EI / rho A)
figure(10)
%h=pcolor(UrExp,FY_tip,10*log(PY_tip)); 
h=pcolor(UrExp,FY_b1il./freqStruct(1),10*log(PY_b1il)); 
xlim([0 50]); ylim([0 50/freqStruct(1)]); %ylim([0 50]);%
hold on
c = colorbar; %c.Label.String = 'PSD';
caxis([-200 1]);
set(h,'edgecolor','none','facecolor','interp');
hold on
plot(UrExp, St*(UrExp*sqrt((E*I)/(ms+ma))*d/l^2/d./freqStruct(1)), 'k--','MarkerSize',2) % Plot Strouhal frequency
hold on
plot(UrExp, 2*St*(UrExp*sqrt((E*I)/(ms+ma))*d/l^2/d./freqStruct(1)), 'k--','MarkerSize',2) % Plot Strouhal frequencyc = colorbar; %c.Label.String = 'PSD';
Urtoplot = [1.8 6.5 13 19.5 29.3 41.5]; 
for u =Urtoplot; xline(u, '-.'); end;
%CF
figure(11)
h=pcolor(UrExp,FY_b1,10*log(PY_b1)); 
xlim([0 50]); ylim([0 50]);
fmodal =[9.6038 2.2688 11.8080 10.1257 9.2411 2.8134 14.0784 41.7157 44.4973 32.9001];
Urtoplot = [1.8 6.5 13 19.5 29.3 41.5]; 
%Urtoplot = [1.8 6.5 19.5]; 
fmodal =[9.6038 2.2688 11.8080 10.1257 9.2411 2.8134 ];
%xline(1.8, ':');xline(6.5, ':'); xline(17, ':');
%xline(25.5, ':');xline(28.3, ':'); xline(48, ':');
for u =Urtoplot; xline(u, '-.'); end;
%for f =fmodal; yline(f, ':'); end;
hold on
plot(UrExp, St*(UrExp*sqrt((E*I)/(ms+ma))*d/l^  2/d), 'k--','MarkerSize',2) % Plot Strouhal frequency
hold on
plot(UrExp, 2*St*(UrExp*sqrt((E*I)/(ms+ma))*d/l^  2/d), 'k--','MarkerSize',2) % Plot Strouhal frequencyc = colorbar; %c.Label.String = 'PSD';
c = colorbar; %c.Label.String = 'PSD';
caxis([-200 1]);
set(h,'edgecolor','none','facecolor','interp');
%set(gca, 'LooseInset', [0,0,0,0]);
stop
%xlabel('Ur');ylabel('f/fw0');title('ONSAS PSD of Y at the tip')
% ylim([0 2])
%ylim([0 50])
figure(15)
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
nexttile
%title('CF Branch 1')
h1=pcolor(UrExp,FY_b1,10*log(PY_b1));
hold on
ylim([0 50])
%subplot(1,2,2)
nexttile
h3=pcolor(UrExp,FY_b1il,10*log(PY_b1il)); 
set(gca,'YTickLabel',{' '})
%title('IL Branch 1')
c = colorbar;
%set(h1,'edgecolor','none','facecolor','interp');
%set(h2,'edgecolor','none','facecolor','interp');
set(h3,'edgecolor','none','facecolor','interp');
%xlabel('Ur');ylabel('f/fw0');title('ONSAS PSD of Y at the tip')
% ylim([0 2])
ylim([0 50])
set(gca, 'LooseInset', [0,0,0,0]);
% figure()
% plot(ydefNumall(end, :)./d)
%
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
connectNodes(figure(4), connection, mesh)
%plot(-xdefNumall(:,end)./d,s,'linewidth',2)
%loc = 3.92/30*Ns;
loc = 0.99*Ns;
if ~isempty(xdefNumall)
    subplot(1,2,1);
    plot(-xdefNumall(:,loc:loc)./d,s,'k-','linewidth',2)
    ylabel('s'); xlabel('x/d');title('In-line displacement')
%     subplot(1,3,2);
%     plot(-xdefNumall(:,loc-10:loc)./d+mean(xdefNumall(:,end/2:end),2)./d,s,'linewidth',2)
%     ylabel('s'); xlabel('x/d');
    subplot(1,2,2);
    plot(ydefNumall(:,loc-200:loc-200)./d,s,'linewidth',2)
    ylabel('s'); xlabel('z/d');title('Cross-flow displacement')
else 
    plot(zdefNumall(:,loc-100:loc)./d,s,'linewidth',2)
end
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
%%
function plotsigend(finalTime, dt, Nval, sig, sigb1, sigb2,xdefNumall, j, tb1, tb2)
       I13 = finalTime-Nval*dt:dt:finalTime;
        figure(j-1) % Transverse
       plot(I13, sig(end-Nval:end)-mean(sig), 'k-')
       hold on 
       plot(I13, sigb1(end-Nval:end)-mean(sigb1), 'b-')
       hold on 
       plot(I13, sigb2(end-Nval:end)-mean(sigb2), 'r-')
       title('Transverse displacements')
       figure(j-2) % Il
       plot(I13, xdefNumall(end, end-Nval:end)-mean(xdefNumall(end, end-Nval:end)), 'k-')
       hold on 
       plot(I13, xdefNumall(tb1, end-Nval:end)- mean(xdefNumall(tb1, end-Nval:end)), 'b-')
       hold on         
       plot(I13, xdefNumall(tb2, end-Nval:end)- mean(xdefNumall(tb2, end-Nval:end)), 'r-')
       title('in-line displacements')
end
%
function plotEnvelope(tiini, step, nsteps, fig,connection, mesh,xdefNumall,ydefNumall,zdefNumall,xyzBool, xzBool)
%tiini = 1000; step = 10; nsteps = 20; 
rose = [1 0.5 0.8]; violet = [0.494 0.184 0.556]; ampfactor=1;
for ti = linspace(tiini,tiini+(nsteps-1)*step,nsteps)
    colorMarker = rose + (violet-rose)*(ti-tiini)/(tiini+(nsteps-1)*step-tiini); 
    %colorMarker = 'k'; % mesh.nodesCoords(:,3) +
    connectNodes(fig, connection, mesh, ampfactor*xdefNumall(:, ti), ampfactor*ydefNumall(:, ti), ampfactor*zdefNumall(:, ti), xyzBool, xzBool, colorMarker)
end
% angle of view
view([-1 -1 1])
set(gca, 'LooseInset', [0,0,0,0]);
end
%
function connectNodes(fig, connection, mesh, dispx, dispy, dispz, xyzBool, xzBool,  colorMarker)
    fig;
    axis equal;
    grid off; box off; 
    set(findobj(gcf, 'type','axes'), 'Visible','off') % hide axis
    set (gca,'Position',[0 0 1 1]);  
    axis off;    
    for i = 2:length(connection)
        node1 = connection(i, 1);
        node2 = connection(i, 2);
        if xzBool; mn1=0; mn2=0; else; mn1=mesh.nodesCoords(node1,2); mn2=mesh.nodesCoords(node2,2); end;
        y1 = mn1 +  dispy(node1);%mesh.nodesCoords(node1,2) +
        y2 = mn2 +  dispy(node2);%mesh.nodesCoords(node2,2) +
        z1 = mesh.nodesCoords(node1,3) + dispz(node1);
        z2 = mesh.nodesCoords(node2,3) + dispz(node2);
        if xyzBool
            x1 = mesh.nodesCoords(node1,1) + dispx(node1);
            x2 = mesh.nodesCoords(node2,1) + dispx(node2);
            plot3([x1,x2], [y1,y2], [z1,z2],'-','Color',colorMarker,'MarkerSize',80,'MarkerFaceColor', colorMarker, 'LineWidth', 0.8);
            %xlabel('x');ylabel('y');zlabel('z');
            if i ==length(connection)
                %hold on
                %plot3(mesh.nodesCoords(101,1) + dispx(101),mesh.nodesCoords(101,2) + dispy(101),mesh.nodesCoords(101,3) + dispz(101),'-o','Color',colorMarker,'MarkerSize',5,'MarkerFaceColor',colorMarker)
            end
        else
            plot([y1,y2], [z1,z2], '-','Color',colorMarker,'MarkerSize',80,'MarkerFaceColor', colorMarker, 'LineWidth', 0.8);
        end
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
    if isempty(lkf)
        lk = NaN;
    elseif Ur<=50
        [val, idx] = max(pkf);
        lk = lkf(idx);
    else
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