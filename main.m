function main()

clc;
clear all;
close all;

%% Single Microphone System --- Generating Noisy Speech Signals

[cs1, fs] = audioread('clean_speech.wav');
cs2 = audioread('clean_speech_2.wav');
aNoise = audioread('aritificial_nonstat_noise.wav');
bNoise = audioread('babble_noise.wav');

% noisySS = cs1 + bNoise(1:length(cs1));
noisySS = cs1 + aNoise;
noisySS = noisySS(1:500000); %Using shorter version of noisy signal

%% Framing using 50% Overlapping Window

wN = 1;                                                                     %Window Number, starting from 1
wT = 0.020;                                                                 %Window time in seconds
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
sWTime = 1.50;                                                              %Sliding window time should be in between 1-2 seconds (page 31 of Book)
M = (sWTime / wT);                                                          %Sliding Window Size (Span of the sliding window calculated by dividing the sliding time by the time of each window)
%In the paper the authors took M = 96 and tracing back the equations they
%have wT = 0.032 seconds giving sWTime approximately equal to 3 seconds but
%we chose as mentioned in slides between 1-2 seconds

nPSD(:,1:M-1) = yPSD_expS(:,1:M-1);
% for k = 1:size(yPSD_expS,1)                                                 %For all the frequency bins
%     for l = 1:M                                                             %For time frames in window M
%         nPSD(k,l) = (yPSD_expS(k,l));                                       %Noise PSD estimation using min of PSD of Yk(l)                
%     end    
% end

for k = 1:size(yPSD_expS,1)                                                 %For all the frequency bins
    for l = M:size(yPSD_expS,2)                                             %For time frames in window M
        nPSD(k,l) = min(yPSD_expS(k,l-M+1:l));                              %Noise PSD estimation using min of PSD of Yk(l)                
    end    
end

% y_PSD_Size = size(yPSD_expS)
% noise_PSD_Size = size(nPSD)



%-------------------------------------------------------------------------- 
%Calculating the bias compensation
%--------------------------------------------------------------------------

%Optimizing alpha for each frame and frequency bin
yPSDSq = yPSD_expS.^2;

yPSD_expS2(:,1) = (0.96).*yPSD(:,1);
yPSD_expSq(:,1) = (0.96).*yPSDSq(:,1);
alpha_max = 0.96;
alpha_min = 0.3;
alpha_Opt = zeros();
alpha_Corr(1) = 1;
alpha_CorrRes(1) = 1;

for l = 2:size(yPSD,2)                                                      %For all the frames
    average_STPSD_PrevFr = sum(yPSD_expS2(:,l-1));
    average_Periodogram = sum(yPSD(:,l));
    alpha_Corr(l) = 1 / (1 + (average_STPSD_PrevFr / average_Periodogram).^2);
    alpha_CorrRes(l) = 0.7*alpha_CorrRes(l-1)+ 0.3*(max(alpha_Corr(l),0.7));
    for k = 1:size(yPSD,1)                                                  %For all the frequency bins
        alpha_Opt(k,l) = (alpha_max*alpha_CorrRes(l))/(1+(yPSD(k,l-1)/nPSD(k,l-1)).^2);                  %Recalculating smoothing factor using the formula (7) given in paper
        alpha_Opt(k,l) = min(alpha_Opt(k,l),alpha_min);
        beta = min((alpha_Opt(k,l)).^2,0.8);
        yPSD_expS2(k,l) = beta.*yPSD_expS2(k,l-1) + (1-beta).*yPSD_expS(k,l);
        yPSD_expSq(k,l) = beta.*yPSD_expSq(k,l-1) + (1-beta).*yPSDSq(k,l);
    end    
end

%Bias Correction Bmin

var_NS = yPSD_expSq - yPSD_expS2.^2;
Qeq(:,1) = 0.5;
nNewPSD(:,1) = nPSD(:,1);
% (2.*(nPSD.^2))./var_NS;

D = M; %For us D = 75 that is close to 80 so using the table given in paper we find MD and HD
MD = 0.865;
HD = 3.25;

for l = 2:size(yPSD,2)                                                      %For all the frames
    for k = 1:size(yPSD,1)                                                  %For all the frequency bins
        Qeq(k,l) = (2*(nNewPSD(k,l-1)).^2) / (var_NS(k,l));
        QeqE(k,l) = (Qeq(k,l) - 2*MD) / (1 - MD);
        Bmin(k,l) = 1 + (D-1)*2/(QeqE(k,l)); %*(gamma((1 + 2./Qeq(k,l)).^HD));
        nNewPSD(k,l) = Bmin(k,l)*nPSD(k,l);
    end
end

%--------------------------------------------------------------------------

figure();
plot(10.*log(abs(yPSD(25,:))),':k');
hold on;
plot(10.*log(abs(yPSD_expS(25,:))),'k');
hold on;
plot(10.*log(abs(nPSD(25,:))),'k','LineWidth',2);
hold on;
plot(10.*log(abs(nNewPSD(25,:))),'r');
hold off;
title('Power Spectral Density, frequency bin = 25')
xlabel('frames')
ylabel('dB')
legend('periodogram (k = 25)','smoothed periodogram (k = 25)','noise estimate without bias compensation (k = 25)','noise estimate with bias compensation (k = 25)')

% figure()
% plot(nNewPSD(25,:))

%% Gain Function



%% IFFT

outWI = (ifft(outWFFT));                                                    %Using transpose to get frames in rows

%% Overlap and Add

procSS = zeros(size(noisySS));                                              %Initializing the size of processed speech signal
procSS(1:wS+1) = outWI(:,1);                                                %Getting back the first frame

for i = 1:numOfWins-1
    procSS(i*wS/2:i*wS/2+wS) = outWI(:,i+1) + procSS(i*wS/2:i*wS/2+wS);     %Overlapping and Adding the received frames to recover the processed speech in procSS
end

% figure()
% spectrogram(procSS,win,'yaxis')
% title('Spectrogram of Enhanced Speech Signal')

% figure()
% subplot(2,1,1)
% plot(noisySS)
% subplot(2,1,2)
% plot(procSS)

end