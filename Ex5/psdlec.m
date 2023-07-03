%% psd - Power spectral density
%
% Compute power spectral density of a time signal
%
% [freq,pow] = psd(X,n_partitions,Fs)
%
% INPUT:   time signal X
%          number of partitions of the whole signal n_partitions
%          sampling frequency Fs
% OUTPUT:  frequency vector freq
%          amplitude vector pow
% Fs = 10000;   
% t = 0:1/Fs:2.96;
% x = cos(2*pi*t*10)+ randn(size(t));
% nfft = 2^nextpow2(length(x));
% Pxx = abs(fft(x,nfft)).^2/length(x)/Fs;
% 
% [freq,pow] = psd2(Pxx(1:length(Pxx)/2),1,Fs); 
% plot(freq,pow)

function [freq,pow] = psdlec(X,n_partitions,Fs)

    N0=length(X);
    % loop on the partitions of the signal
    for p=1:n_partitions
        x=X((p-1)*N0/n_partitions+1:p*N0/n_partitions);% isolate portion of the signal to consider
        N = length(x);
        xdft = fft(x);                                  % fft of the signal portion
        xdft = xdft(1:floor(N/2)+1);                           % consider only left side of the spectrum
        xpsd = (1/(Fs*N)) * abs(xdft).^2;               % power magnitude
        xpsd(2:end-1) = 2*xpsd(2:end-1);
        % add power spectrums of the different portions
        if p==1;
            pow=xpsd;
        else
          pow=pow+xpsd;  
        end
    end
    
    freq = 0:Fs/N:Fs/2;                                 % frequency vector

end

