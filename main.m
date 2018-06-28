function main()

clc;
clear all;
% close all;

%% Single Microphone System --- Generating Noisy Speech Signals

SNR = 10;
[cs1, fs] = audioread('clean_speech.wav');
cs2 = audioread('clean_speech_2.wav');
aNoise = audioread('aritificial_nonstat_noise.wav');
bNoise = audioread('babble_noise.wav');
cNoise = audioread('Speech_shaped_noise.wav');

noise = aNoise;

% noisySS = cs1 + cNoise(1:length(cs1));
% noisySS = noisySS(1:800000); %Using shorter version of noisy signal
% soundsc(noisySS, fs)

%Power of signal
sig_Pow = sum((abs(cs1).^2)./length(cs1));
noi_Pow = sum((abs(noise).^2)./length(noise));

snrOrig = 10*log10(sig_Pow./noi_Pow) %original SNR
noi_Des = sig_Pow./(10^(SNR/10));
noiSNR = sqrt(noi_Des/noi_Pow).*noise;
noisySS = cs1 + noiSNR;

audiowrite('noisySS.wav',noisySS,fs)
%% Framing using 50% Overlapping Window

wN = 1;                                                                     %Window Number, starting from 1
wT = 0.025;                                                                 %Window time in seconds
wS = wT*fs;                                                                 %Window size = Number of Samples in Window
win = hann(wS+1);                                                           %Generating window, that can be hanning or hamming, of size wS samples   
numOfWins = ceil(length(noisySS)/wS)*2;                                     %Total number of windows 

%Appending zeros at the end of Audio Signal to avoid dimesion mismatch for last window
if(numel(noisySS) < ((numOfWins+1)/2)*wS)
    noisySS(numel(noisySS)+1:((numOfWins+1)/2)*wS) = 0;   
end

%Storing the frames of noisy signal, y(t), into columns of outWO
outWO(:,1) = noisySS(1:wS+1).*win;                                          %Storing the first frame                                          

for i = 1:numOfWins-1                                                       %Looping to store other frames
        outWO(:,i+1) = noisySS(i*wS/2:i*wS/2+wS).*win;
%         plot((i*wS/2:i*wS/2+wS),outWO(:,i+1));
%         hold on;    
end

%Plotting the spectrogram of the noisySpeechSignal
% figure()
% spectrogram(noisySS,win,'yaxis')
% title('Spectrogram of unprocessed Noisy Speech Signal')

%% FFT
size(outWO)
outWFFT = fft(outWO);                                                      %FFT works on columns so take transpose of outWO
% sizeoutWO = size(outWO)
% sizeoutWFFT = size(outWFFT)

%% Noise PSD Estimator

yPSD = (abs(outWFFT).^2);                                                   %Calculating PSD for each frame

%--------------------------------------------------------------------------
%Reducing Variance in signal Yk(l) using Exponential Smoother
%--------------------------------------------------------------------------
tSm = 0.2;                                                                  %Smoothing window time = 0.2 seconds as given in paper
alpha = ((tSm*fs / (wS/2)) - 1) / ((tSm*fs / (wS/2)) + 1)                   %Alpha determines the smoothing factor calculated using formula from the paper in ref [38] of the book
yPSD_expS(:,1) = (1-alpha).*yPSD(:,1);                                      %Reducing Variance

for i = 2:size(yPSD,2)
    yPSD_expS(:,i) = alpha.*yPSD_expS(:,i-1) + (1-alpha).*yPSD(:,i);
end

%-------------------------------------------------------------------------- 
%Noise PSD Estimation using Minimum Statistics
%--------------------------------------------------------------------------
sWTime = 1.00;                                                              %Sliding window time should be in between 1-2 seconds (page 31 of Book)
M = (sWTime / wT);                                                          %Sliding Window Size (Span of the sliding window calculated by dividing the sliding time by the time of each window)
%In the paper the authors took M = 96 and tracing back the equations they
%have wT = 0.032 seconds giving sWTime approximately equal to 3 seconds but
%we chose as mentioned in slides between 1-2 seconds

nPSD(:,1:M-1) = yPSD_expS(:,1:M-1);

for k = 1:size(yPSD_expS,1)                                                 %For all the frequency bins
    for l = M:size(yPSD_expS,2)                                             %For time frames in window M
        nPSD(k,l) = min(yPSD_expS(k,l-M+1:l));                              %Noise PSD estimation using min of PSD of Yk(l)                
    end    
end

% y_PSD_Size = size(yPSD_expS)
% noise_PSD_Size = size(nPSD)

figure();
plot(10.*log(abs(yPSD(25,:))),':k');
hold on;
plot(10.*log(abs(yPSD_expS(25,:))),'k');
hold on;
plot(10.*log(abs(nPSD(25,:))),'k','LineWidth',2);
hold on;

%-------------------------------------------------------------------------- 
%Calculating the bias compensation
%--------------------------------------------------------------------------

%Optimizing alpha for each frame and frequency bin

yPSD_expS(:,1) = (1-alpha).*yPSD(:,1);
yPSD_expS2(:,1) = (1-0.8^2)*yPSD_expS(:,1);
yPSD_expSq(:,1) = yPSD_expS(:,1).^2;

alpha_max = 0.96;
alpha_min = 0.30;
alpha_Opt(1:size(yPSD,1),1) = 0.50;
alpha_Corr(1) = 0.5;
alpha_CorrRes(1) = 0.5;

for l = 2:size(yPSD,2)                                                      %For all the frames
    average_STPSD_PrevFr = mean(yPSD_expS(:,l-1));
    average_Periodogram = mean(yPSD(:,l));
    alpha_Corr(l) = 1 / (1 + (average_STPSD_PrevFr / average_Periodogram - 1).^2);
    alpha_CorrRes(l) = 0.7*alpha_CorrRes(l-1)+ 0.3*(max(alpha_Corr(l),0.7));
    for k = 1:size(yPSD,1)                                                  %For all the frequency bins
        alpha_Opt(k,l) = min(alpha_max,(alpha_CorrRes(l))/(1+(yPSD_expS(k,l-1)/nPSD(k,l-1)-1).^2));                  %Recalculating smoothing factor using the formula (7) given in paper
        alpha_Opt(k,l) = alpha_max*(alpha_CorrRes(l))/(1+(yPSD_expS(k,l-1)/nPSD(k,l-1)-1).^2);                  %Recalculating smoothing factor using the formula (7) given in paper
        alpha_Opt(k,l) = max(alpha_Opt(k,l),alpha_min);
        beta = min((alpha_Opt(k,l).^2),0.8);
        yPSD_expS(k,l) = alpha_Opt(k,l).*yPSD_expS(k,l-1) + (1-alpha_Opt(k,l)).*yPSD(k,l);
        yPSD_expS2(k,l) = beta.*yPSD_expS2(k,l-1) + (1-beta).*yPSD_expS(k,l);
        yPSD_expSq(k,l) = beta.*yPSD_expSq(k,l-1) + (1-beta).*(yPSD_expS(k,l).^2);
    end    
end

%Bias Correction Bmin

nPSD(:,1:M-1) = yPSD_expS(:,1:M-1);

for k = 1:size(yPSD_expS,1)                                                 %For all the frequency bins
    for l = M:size(yPSD_expS,2)                                             %For time frames in window M
        nPSD(k,l) = min(yPSD_expS(k,l-M+1:l));                              %Noise PSD estimation using min of PSD of Yk(l)                
    end    
end

D = M; %For us D = 75 that is close to 80 so using the table given in paper we find MD and HD
MD = 0.865;
HD = 3.25;

var_NS = yPSD_expSq - yPSD_expS2.^2;
Qeq(:,1) = 2;
QeqE(:,1) = 2;
nNewPSD(:,1) = D.*nPSD(:,1);    %In paper it says Bmin = D for Qeq = 2
 
for l = 2:size(yPSD,2)                                                      %For all the frames
    for k = 1:size(yPSD,1)                                                  %For all the frequency bins
        Qeq(k,l) = (2*(nNewPSD(k,l-1)).^2) / (var_NS(k,l));
%         Qeq(k,l) = (2*(nNewPSD(k,l-1)).^2) / abs(var_NS(k,l));        
        QeqE(k,l) = (Qeq(k,l) - 2*MD) / (1 - MD);
        Bmin(k,l) = 1 + ((D-1)*2/(QeqE(k,l)));%*((gamma(1 + 2/Qeq(k,l))).^HD);
        nNewPSD(k,l) = nPSD(k,l);
        nNewPSD(k,l) = Bmin(k,l)*nPSD(k,l);
    end
end
% 
% %--------------------------------------------------------------------------
plot(10.*log(abs(yPSD_expS(25,:))),'g');
hold on;
plot(10.*log(abs(nPSD(25,:))),'r');
hold on;
plot(10.*log(abs(nNewPSD(25,:))),'r');
hold on;    
% title('Power Spectral Density, frequency bin = 25')
% xlabel('frames')
% ylabel('dB')
% legend('periodogram (k = 25)','smoothed periodogram (k = 25)','noise estimate without bias compensation (k = 25)','noise estimate with bias compensation (k = 25)')

% plot(alpha_Opt(25,:),'b');
% hold on;
%% Gain Function
%--------------------------------------------------------------------------
%Speech PSD Estimation using Wiener Filter 
%--------------------------------------------------------------------------

for r = 1:size(yPSD,1)
    for c = 1:size(yPSD,2)
        filt_Response(r,c) = max((1 - abs(nNewPSD(r,c)) / yPSD(r,c)),0);
        PSD_Speech(r,c) = (filt_Response(r,c).*abs(yPSD(r,c)));
        speech_est_Wiener(r,c) = ((filt_Response(r,c))*abs(outWFFT(r,c)))*exp(complex(0,angle(outWFFT(r,c))));
    end
end



%Spectral Subtraction

sEstSS = (max((outWFFT-abs((nNewPSD))),0)).*exp(complex(0,angle(outWFFT)));


%% IFFT

% outWI = real(ifft(sEstSS));                                                    %Using transpose to get frames in rows
outWI = real(ifft(speech_est_Wiener));                                                    %Using transpose to get frames in rows
% size(outWI)
%% Overlap and Add

procSS = zeros(size(noisySS));                                              %Initializing the size of processed speech signal
procSS(1:wS+1) = outWI(:,1);                                                %Getting back the first frame

for i = 1:numOfWins-1
    procSS(i*wS/2:i*wS/2+wS) = outWI(:,i+1) + procSS(i*wS/2:i*wS/2+wS);     %Overlapping and Adding the received frames to recover the processed speech in procSS
end

% figure()
% spectrogram(procSS,win,'yaxis')
% title('Spectrogram of Enhanced Speech Signal')

figure()
subplot(2,1,1)
plot(noisySS)
ylim([-0.5 0.5])
subplot(2,1,2)
plot(procSS)
ylim([-0.5 0.5])

audiowrite('procSS.wav',procSS,fs)

end