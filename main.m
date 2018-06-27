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
noisySS = noisySS(1:50000); %Using shorter version of noisy signal

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
figure()
spectrogram(noisySS,win,'yaxis')
title('Spectrogram of unprocessed Noisy Speech Signal')

%% FFT

outWFFT = fft(outWO);                                                      %FFT works on columns so take transpose of outWO
% sizeoutWO = size(outWO)
% sizeoutWFFT = size(outWFFT)

%% Noise PSD Estimator

yPSD = (abs(outWFFT).^2);                                                   %Calculating PSD for each frame

%--------------------------------------------------------------------------
%Reducing Variance in signal Yk(l) using Exponential Smoother
%--------------------------------------------------------------------------
alpha = 0.90;                                                               %Alpha determines smoothness factor
yPSD_expS(:,1) = (1-alpha).*yPSD(:,1);                                      %Reducing Variance

for i = 2:size(yPSD,2)
    yPSD_expS(:,i) = alpha.*yPSD_expS(:,i-1) + (1-alpha).*yPSD(:,i);
end

%-------------------------------------------------------------------------- 
%Noise PSD Estimation using Minimum Statistics
%--------------------------------------------------------------------------
sWTime = 1.5;                                                               %Sliding window time should be in between 1-2 seconds (page 31 of Book)
M = (sWTime / wT);                                                          %Sliding Window Size

for k = 1:size(yPSD_expS,1)                                                 %For all the frequency bins
    for l = 1:size(yPSD_expS,2)-M+1                                         %For time frames in window M
        nPSD(k,l) = min(yPSD_expS(k,l:l+M-1));                              %Noise PSD estimation using min of PSD of Yk(l)                
    end    
end

%For the last frames that go outside of the sliding window size
for k = 1:size(yPSD_expS,1)
    nPSD(k,size(yPSD_expS,2)-M+2:size(yPSD_expS,2)) = min(yPSD_expS(k,size(yPSD_expS,2)-M+2:size(yPSD_expS,2)));
end

% y_PSD_Size = size(yPSD_expS)
% noise_PSD_Size = size(nPSD)

%% Gain Function



%% IFFT

outWI = (ifft(outWFFT));                                                    %Using transpose to get frames in rows

%% Overlap and Add

procSS = zeros(size(noisySS));                                              %Initializing the size of processed speech signal
procSS(1:wS+1) = outWI(:,1);                                                %Getting back the first frame

for i = 1:numOfWins-1
    procSS(i*wS/2:i*wS/2+wS) = outWI(:,i+1) + procSS(i*wS/2:i*wS/2+wS);     %Overlapping and Adding the received frames to recover the processed speech in procSS
end

figure()
spectrogram(procSS,win,'yaxis')
title('Spectrogram of Enhanced Speech Signal')

figure()
subplot(2,1,1)
plot(noisySS)
subplot(2,1,2)
plot(procSS)

end