function main()

clear all;
close all;
% figure()
%% Reading Audio Signals and Adding Noise

[cs1, fs] = audioread('clean_speech.wav');
cs2 = audioread('clean_speech_2.wav');
aNoise = audioread('aritificial_nonstat_noise.wav');
bNoise = audioread('babble_noise.wav');

noisySS = cs1 + aNoise;
noisySS = noisySS(1:20000); %Using shorter version of noisy signal
subplot(2,1,1)
plot(noisySS);

%% Framing using 50% Overlapping Window

wN = 1; %Window Number, starting from 1
wT = 0.025; %Window time in seconds
wS = wT*fs %Window size = Number of Samples in Window

%Selecting odd window size
% if (mod(wS,2) == 0)
%     wS = wS + 1
% end

dem = 1*ones(20000,1);
% noisySS = dem;
sd = size(dem);
% stem (dem)
% hold on;
% L = 17;
win = hann(wS+1);
% win = ones(wS+1,1);
% sw = size(win)
% stem(win)

numOfWins = ceil(length(noisySS)/wS)*2-1;

if(numel(noisySS) < ((numOfWins+1)/2)*wS)
    noisySS(numel(noisySS)+1:((numOfWins+1)/2)*wS) = 0;   
end
plot(noisySS)
hold on;

% subplot(2,1,2)
outWO(:,1) = noisySS(1:wS+1).*win;
% plot((1:wS+1),outWO(:,1));
% hold on;

for i = 1:numOfWins-1
        outWO(:,i+1) = noisySS(i*wS/2:i*wS/2+wS).*win;
        x_axi = (i*wS/2:i*wS/2+wS);
        plot((i*wS/2:i*wS/2+wS),outWO(:,i+1));
        hold on;    
end

%% FFT

outWFFT = fft(outWO);

%% Gain Function

%% IFFT
outWI = ifft(outWFFT);

%% Overlap and Add
procSS = zeros(size(noisySS));
procSS(1:wS+1) = outWI(:,1);

for i = 1:numOfWins-1
%     numel(i*wS/2:i*wS/2+wS)
%     size(outWI(:,i+1))
    procSS(i*wS/2:i*wS/2+wS) = outWI(:,i+1) + procSS(i*wS/2:i*wS/2+wS);
%         x_axi = (i*wS/2:i*wS/2+wS);
%         plot((i*wS/2:i*wS/2+wS),outWO(:,i+1));
%         hold on;    
end

figure()
subplot(2,1,1)
plot(noisySS)
subplot(2,1,2)
plot(procSS)


end