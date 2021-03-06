function main()

clc;
clear all;
close all;

%% Single Microphone System --- Generating Noisy Speech Signals

SNR = 30;
[cs1, fs] = audioread('clean_speech.wav');
cs2 = audioread('clean_speech_2.wav');
aNoise = audioread('aritificial_nonstat_noise.wav');
bNoise = audioread('babble_noise.wav');
cNoise = audioread('Speech_shaped_noise.wav');

% noise = cNoise(1:length(cs1));
speech = cs1;
noise = aNoise;

%Power of signal
sig_Pow = sum((abs(speech).^2)./length(speech));
noi_Pow = sum((abs(noise).^2)./length(noise));

snrOrig = 10*log10(sig_Pow./noi_Pow); %original SNR
noi_Des = sig_Pow./(10^(SNR/10));
noiSNR = sqrt(noi_Des/noi_Pow).*noise;
noisySS = speech + noiSNR;

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
end

%Plotting the spectrogram of the noisySpeechSignal
% figure()
% spectrogram(noisySS,win,'yaxis')
% title('Spectrogram of unprocessed Noisy Speech Signal')

%% FFT

outWFFT = fft(outWO);                                                      %FFT works on columns so take transpose of outWO

%% Noise PSD Estimator

yPSD = (abs(outWFFT).^2);                                                   %Calculating PSD for each frame

%--------------------------------------------------------------------------
%Reducing Variance in signal Yk(l) using constant Exponential Smoother
%--------------------------------------------------------------------------

tSm = 0.2;                                                                  %Smoothing window time = 0.2 seconds as given in paper
alpha = ((tSm*fs / (wS/2)) - 1) / ((tSm*fs / (wS/2)) + 1);                  %Alpha determines the smoothing factor calculated using formula from the paper in ref [38] of the book
yPSD_expS(:,1) = (1-alpha).*yPSD(:,1);                                      %Reducing Variance

for i = 2:size(yPSD,2)
    yPSD_expS(:,i) = alpha.*yPSD_expS(:,i-1) + (1-alpha).*yPSD(:,i);
end

%-------------------------------------------------------------------------- 
%Noise PSD Estimation using Minimum Statistics
%--------------------------------------------------------------------------

sWTime = 2.00;                                                              %Sliding window time should be in between 1-2 seconds (page 31 of Book)
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

nPSD = max(nPSD,0);

figure(1);
plot(10.*log(abs(yPSD(25,:))),':k');
hold on;
plot(10.*log(abs(yPSD_expS(25,:))),'k');
hold on;
plot(10.*log(abs(nPSD(25,:))),'k','LineWidth',2);
hold off;
title('Power Spectral Density, frequency bin = 25, before Bias Compensation')
xlabel('frames')
ylabel('dB')
legend('periodogram', (sprintf('smoothed periodogram using constant alpha = %.3g', alpha)), 'noise estimate without bias compensation')

%-------------------------------------------------------------------------- 
%Calculating the bias compensation
%--------------------------------------------------------------------------

%Optimizing alpha for each frame and frequency bin

yPSD_expSnew(:,1) = (1-alpha).*yPSD(:,1);
yPSD_expS2(:,1) = (1-0.8^2)*yPSD_expSnew(:,1);
yPSD_expSq(:,1) = (1-0.8^2).*yPSD_expSnew(:,1).^2;

alpha_max = 0.96;
alpha_min = 0.30;
alpha_Opt(1:size(yPSD,1),1) = 0.50;
alpha_Corr(1) = 0.5;
alpha_CorrRes(1) = 0.5;

for l = 2:size(yPSD,2)                                                      %For all the frames
    average_STPSD_PrevFr = mean(yPSD_expSnew(:,l-1));
    average_Periodogram = mean(yPSD(:,l));
    alpha_Corr(l) = 1 / (1 + (average_STPSD_PrevFr / average_Periodogram - 1).^2);
    alpha_CorrRes(l) = 0.7*alpha_CorrRes(l-1)+ 0.3*(max(alpha_Corr(l),0.7));
    for k = 1:size(yPSD,1)                                                  %For all the frequency bins
        alpha_Opt(k,l) = (alpha_CorrRes(l)).*min(alpha_max,1/(1+(yPSD_expS(k,l-1)/nPSD(k,l-1)-1).^2));                  %Recalculating smoothing factor using the formula (7) given in paper
        alpha_Opt(k,l) = max(alpha_Opt(k,l),alpha_min);
        beta = min((alpha_Opt(k,l).^2),0.8);
        yPSD_expSnew(k,l) = alpha_Opt(k,l).*yPSD_expSnew(k,l-1) + (1-alpha_Opt(k,l)).*yPSD(k,l);
        yPSD_expS2(k,l) = beta.*yPSD_expS2(k,l-1) + (1-beta).*yPSD_expS(k,l);
        yPSD_expSq(k,l) = beta.*yPSD_expSq(k,l-1) + (1-beta).*(yPSD_expS(k,l).^2);
    end    
end

D = M; %For us D = 80 that is close to 80 so using the table given in paper we find MD and HD
MD = 0.865;
HD = 3.25;

var_NS = yPSD_expSq - yPSD_expS2.^2;
Qeq(:,1) = 2;
QeqE(:,1) = 2;
nNewPSD(:,1) = D.*nPSD(:,1);    %In paper it says Bmin = D for Qeq = 2
 
for l = 2:size(yPSD,2)                                                      %For all the frames
    for k = 1:size(yPSD,1)                                                  %For all the frequency bins
        Qeq(k,l) = min((2*(nPSD(k,l-1)).^2) / (var_NS(k,l)), 2);            %Because mentioned in paper
%         Qeq(k,l) = min((2*(nNewPSD(k,l-1)).^2) / abs(var_NS(k,l)),50);        
        QeqE(k,l) = (Qeq(k,l) - 2*MD) / (1 - MD);
        Bmin(k,l) = 1 + ((D-1)*2/(QeqE(k,l)));%*(gamma((1 + 2/Qeq(k,l)).^HD));
        nNewPSD(k,l) = max(Bmin(k,l)*nPSD(k,l),0);
    end
end

%--------------------------------------------------------------------------
figure(2)
plot(10.*log(abs(yPSD(25,:))),':k');
hold on;
plot(10.*log(abs(yPSD_expSnew(25,:))),'k');
hold on;
plot(10.*log(abs(nNewPSD(25,:))),'r');
hold on
plot(10.*log(abs(nPSD(25,:))),'b');
hold off;
title('Power Spectral Density, frequency bin = 25, after Optimal Smoothing and Bias Compensation')
xlabel('frames')
ylabel('dB')
legend('periodogram', 'smoothed periodogram using optimal smoothing', 'noise estimate with bias compensation','noise estimate without bias compensation')

%% Gain Function

%Either choose Noise as estimation without bias or with bias compensation
%by uncommenting or commenting nNewPSD = nPSD below, respectively.

% nNewPSD = nPSD;


%--------------------------------------------------------------------------
%Speech PSD Estimation using Wiener Filter 
%--------------------------------------------------------------------------

PSD_Speech = max(yPSD_expSnew - nNewPSD,0);
        
for r = 1:size(yPSD,1)
    for c = 1:size(yPSD,2)
        filt_Res(r,c) = PSD_Speech(r,c) / (yPSD_expSnew(r,c));
        speech_est_Wiener(r,c) = filt_Res(r,c)*outWFFT(r,c);
    end
end

%--------------------------------------------------------------------------
%Simple Spectral Subtraction
%--------------------------------------------------------------------------

a = 2;
b = 3;

sEstSS = ((max((1 - (b.*abs(nNewPSD).^a)./((abs(outWFFT)).^a)),0)).^(1/a)).*outWFFT;

%--------------------------------------------------------------------------
%MMSE 
%--------------------------------------------------------------------------

kappa = yPSD./nPSD;
sigmaSNR(:,1) = kappa(:,1);
alphaSNR = 0.98;    %Typical Value
sigmaSNR(:,1) = (1-alphaSNR).*max(kappa(:,1)-1,0);
for i = 2:size(yPSD,2)
    sigmaSNR(:,i) = alphaSNR*(4/pi).*(PSD_Speech(:,i-1)./nPSD(:,i)) + (1-alphaSNR).*max(kappa(:,i)-1,0);
end

% for i = 2:size(yPSD,2)
%     sigmaSNR(:,i) = (yPSD(:,i)./nPSD(:,i-1)) - 1;                        %Slide 18 Lecture 4, Using Maximum Likelihood Estimator we get noise only at the output
% end

u = (sigmaSNR./(sigmaSNR + 1)).*kappa;
HMMSE = ((sqrt(pi*u))./(2*kappa)).*(exp(-u/2)).*((1+u).*besseli(0,u./2) + u.*besseli(1,u./2));
sEMMSE = max((HMMSE.*abs(outWFFT)).*exp(complex(0,angle(outWFFT))),0);
% sEMMSE = (outWFFT - max(HMMSE.*abs(outWFFT),0)).*exp(complex(0,angle(outWFFT)));

%% IFFT

outWI(:,:,1) = real(ifft(sEstSS));                                                    
outWI(:,:,2) = real(ifft(speech_est_Wiener));
outWI(:,:,3) = real(ifft(sEMMSE));                                                    

%% Overlap and Add

orig = stoi(speech, noisySS(1:length(speech)), fs);
% fprintf('STOI score unprocessed: %.4f\n', orig);
figure(3)
subplot(4,1,1)
plot(noisySS)
title(sprintf('Noisy Speech Signal with SNR = %.4g dB and STOI Score = %.4g', SNR, orig))
ylim([-0.5 0.5])
set(gca,'XTick',[]);

for filt = 1:3
    procSS = zeros(size(noisySS));                                              %Initializing the size of processed speech signal
    procSS(1:wS+1) = outWI(:,1,filt);                                                %Getting back the first frame

    for i = 1:numOfWins-1
        procSS(i*wS/2:i*wS/2+wS) = outWI(:,i+1,filt) + procSS(i*wS/2:i*wS/2+wS);     %Overlapping and Adding the received frames to recover the processed speech in procSS
    end

    e = (speech - procSS(1:length(speech))).^2;
    eNan = isnan(e);
    eI = find(eNan == 0);
    er = (sum(e(eI)));
    ee = (e(eI));
    
    error = sqrt(sum(e(eI)));
%     fprintf('Sum of squared errors: %.4f\n', error);
    intelligibility = stoi(speech, procSS(1:length(speech)), fs);
%     fprintf('STOI score   processed: %.4f\n', intelligibility);
    figure(3)
    subplot(4,1,filt+1)
    plot(procSS)
    ylim([-0.5 0.5])
    if filt < 4
      set(gca,'XTick',[]);
    end   
    switch filt
        case 1
            title(sprintf('Enhanced Speech Signal using Spectral Subtraction with STOI Score = %.4g and SSE = %.4g', intelligibility, error));
        case 2
            title(sprintf('Enhanced Speech Signal using Wiener Filtering with STOI Score = %.4g and SSE = %.4g', intelligibility, error));
        case 3
            title(sprintf('Enhanced Speech Signal using MMSE Filtering with STOI Score = %.4g and SSE = %.4g', intelligibility, error));
    end    
    audiowrite(strcat('procSS', num2str(filt),'.wav'),procSS,fs)
end
% figure()
% spectrogram(procSS,win,'yaxis')
% title('Spectrogram of Enhanced Speech Signal')


% %% Multi-Mic System 

% speech2 = cs2;
% noise2 = cNoise(1:length(cs2));

% sig_Pow2 = sum((abs(speech2).^2)./length(speech2));
% noi_Pow2 = sum((abs(noise2).^2)./length(noise2));

% snrOrig2 = 10*log10(sig_Pow2./noi_Pow2); %original SNR
% noi_Des2 = sig_Pow2./(10^(SNR/10));
% noiSNR2 = sqrt(noi_Des2/noi_Pow2).*noise2;
% noisySS2 = speech2 + noiSNR2;

% sound_speed = 330;                                                          %in meter / second
% alpha = 40;                                                                 %direct signal beam
% beta = 180;                                                                 %undesired signal
% theta = [alpha beta];
% dist = 0.002;                                                                  %distance between 2 mics, in meters 
% numOfMics = 2;

% %Appending zeros at the end of Audio Signal to avoid dimesion mismatch for last window
% if(numel(noisySS2) < ((numOfWins+1)/2)*wS)
    % noisySS2(numel(noisySS2)+1:((numOfWins+1)/2)*wS) = 0;   
% end

% %Storing the frames of noisy signal, y(t), into columns of outWO
% outWO2(:,1) = noisySS2(1:wS+1).*win;                                          %Storing the first frame                                          

% for i = 1:numOfWins-1                                                       %Looping to store other frames
        % outWO2(:,i+1) = noisySS2(i*wS/2:i*wS/2+wS).*win;   
% end

% outWFFT2 = fft(outWO2);                                                      %FFT works on columns so take transpose of outWO

% numOfSamples = min(size(outWO,2),size(outWO2,2));
% S = [outWFFT; outWFFT2];
% numFBins = size((outWFFT),1);

% X = gen_data(S,2,numOfSamples,dist,theta,sound_speed,numFBins,fs);

% Delta = (dist.*((0:numFBins-1)*fs/numFBins))./(sound_speed);
% angleS = theta;

% for t = 1:numel(angleS)
    % W_res(:,((t-1)*length(Delta)+1):t*length(Delta)) = gen_a(2,Delta,angleS(t));
% end

% w = W_res;
% % size(X)
% size(w)
% wH = ctranspose(w);

% stk = (wH(1:size((W_res),2)/2,:)*X);
% outWO = stk;

% % %% FFT
% % 
% outWFFT = fft(outWO);                                                      %FFT works on columns so take transpose of outWO

% %% Noise PSD Estimator

% yPSD = (abs(outWFFT).^2);                                                   %Calculating PSD for each frame

% %--------------------------------------------------------------------------
% %Reducing Variance in signal Yk(l) using constant Exponential Smoother
% %--------------------------------------------------------------------------

% tSm = 0.2;                                                                  %Smoothing window time = 0.2 seconds as given in paper
% alpha = ((tSm*fs / (wS/2)) - 1) / ((tSm*fs / (wS/2)) + 1);                  %Alpha determines the smoothing factor calculated using formula from the paper in ref [38] of the book
% yPSD_expS(:,1) = (1-alpha).*yPSD(:,1);                                      %Reducing Variance

% for i = 2:size(yPSD,2)
    % yPSD_expS(:,i) = alpha.*yPSD_expS(:,i-1) + (1-alpha).*yPSD(:,i);
% end

% %-------------------------------------------------------------------------- 
% %Noise PSD Estimation using Minimum Statistics
% %--------------------------------------------------------------------------

% % sWTime = 2.00;                                                              %Sliding window time should be in between 1-2 seconds (page 31 of Book)
% % M = (sWTime / wT);                                                          %Sliding Window Size (Span of the sliding window calculated by dividing the sliding time by the time of each window)
% %In the paper the authors took M = 96 and tracing back the equations they
% %have wT = 0.032 seconds giving sWTime approximately equal to 3 seconds but
% %we chose as mentioned in slides between 1-2 seconds

% nPSD(:,1:M-1) = yPSD_expS(:,1:M-1);

% for k = 1:size(yPSD_expS,1)                                                 %For all the frequency bins
    % for l = M:size(yPSD_expS,2)                                             %For time frames in window M
        % nPSD(k,l) = min(yPSD_expS(k,l-M+1:l));                              %Noise PSD estimation using min of PSD of Yk(l)                
    % end    
% end

% nPSD = max(nPSD,0);

% %-------------------------------------------------------------------------- 
% %Calculating the bias compensation
% %--------------------------------------------------------------------------

% %Optimizing alpha for each frame and frequency bin

% yPSD_expSnew(:,1) = (1-alpha).*yPSD(:,1);
% yPSD_expS2(:,1) = (1-0.8^2)*yPSD_expSnew(:,1);
% yPSD_expSq(:,1) = (1-0.8^2).*yPSD_expSnew(:,1).^2;

% alpha_max = 0.96;
% alpha_min = 0.30;
% alpha_Opt(1:size(yPSD,1),1) = 0.50;
% alpha_Corr(1) = 0.5;
% alpha_CorrRes(1) = 0.5;

% for l = 2:size(yPSD,2)                                                      %For all the frames
    % average_STPSD_PrevFr = mean(yPSD_expSnew(:,l-1));
    % average_Periodogram = mean(yPSD(:,l));
    % alpha_Corr(l) = 1 / (1 + (average_STPSD_PrevFr / average_Periodogram - 1).^2);
    % alpha_CorrRes(l) = 0.7*alpha_CorrRes(l-1)+ 0.3*(max(alpha_Corr(l),0.7));
    % for k = 1:size(yPSD,1)                                                  %For all the frequency bins
        % alpha_Opt(k,l) = (alpha_CorrRes(l)).*min(alpha_max,1/(1+(yPSD_expS(k,l-1)/nPSD(k,l-1)-1).^2));                  %Recalculating smoothing factor using the formula (7) given in paper
        % alpha_Opt(k,l) = max(alpha_Opt(k,l),alpha_min);
        % beta = min((alpha_Opt(k,l).^2),0.8);
        % yPSD_expSnew(k,l) = alpha_Opt(k,l).*yPSD_expSnew(k,l-1) + (1-alpha_Opt(k,l)).*yPSD(k,l);
        % yPSD_expS2(k,l) = beta.*yPSD_expS2(k,l-1) + (1-beta).*yPSD_expS(k,l);
        % yPSD_expSq(k,l) = beta.*yPSD_expSq(k,l-1) + (1-beta).*(yPSD_expS(k,l).^2);
    % end    
% end

% D = M; %For us D = 80 that is close to 80 so using the table given in paper we find MD and HD
% MD = 0.865;
% HD = 3.25;

% var_NS = yPSD_expSq - yPSD_expS2.^2;
% Qeq(:,1) = 2;
% QeqE(:,1) = 2;
% nNewPSD(:,1) = D.*nPSD(:,1);    %In paper it says Bmin = D for Qeq = 2
 
% for l = 2:size(yPSD,2)                                                      %For all the frames
    % for k = 1:size(yPSD,1)                                                  %For all the frequency bins
        % Qeq(k,l) = min((2*(nPSD(k,l-1)).^2) / (var_NS(k,l)), 2);            %Because mentioned in paper
% %         Qeq(k,l) = min((2*(nNewPSD(k,l-1)).^2) / abs(var_NS(k,l)),50);        
        % QeqE(k,l) = (Qeq(k,l) - 2*MD) / (1 - MD);
        % Bmin(k,l) = 1 + ((D-1)*2/(QeqE(k,l)));%*(gamma((1 + 2/Qeq(k,l)).^HD));
        % nNewPSD(k,l) = max(Bmin(k,l)*nPSD(k,l),0);
    % end
% end

% %% Gain Function

% %Either choose Noise as estimation without bias or with bias compensation
% %by uncommenting or commenting nNewPSD = nPSD below, respectively.

% %--------------------------------------------------------------------------
% %Speech PSD Estimation using Wiener Filter 
% %--------------------------------------------------------------------------

% PSD_Speech = max(yPSD_expSnew - nNewPSD,0);
        
% for r = 1:size(yPSD,1)
    % for c = 1:size(yPSD,2)
        % filt_Res(r,c) = PSD_Speech(r,c) / (yPSD_expSnew(r,c));
        % speech_est_Wiener(r,c) = filt_Res(r,c)*outWFFT(r,c);
    % end
% end

% %--------------------------------------------------------------------------
% %Simple Spectral Subtraction
% %--------------------------------------------------------------------------

% a = 2;
% b = 3;

% sEstSS = ((max((1 - (b.*abs(nNewPSD).^a)./((abs(outWFFT)).^a)),0)).^(1/a)).*outWFFT;

% %--------------------------------------------------------------------------
% %MMSE 
% %--------------------------------------------------------------------------

% kappa = yPSD./nPSD;
% sigmaSNR(:,1) = kappa(:,1);
% alphaSNR = 0.98;    %Typical Value
% sigmaSNR(:,1) = (1-alphaSNR).*max(kappa(:,1)-1,0);
% for i = 2:size(yPSD,2)
    % sigmaSNR(:,i) = alphaSNR*(4/pi).*(PSD_Speech(:,i-1)./nPSD(:,i)) + (1-alphaSNR).*max(kappa(:,i)-1,0);
% end

% % for i = 2:size(yPSD,2)
% %     sigmaSNR(:,i) = (yPSD(:,i)./nPSD(:,i-1)) - 1;                        %Slide 18 Lecture 4, Using Maximum Likelihood Estimator we get noise only at the output
% % end

% u = (sigmaSNR./(sigmaSNR + 1)).*kappa;
% HMMSE = ((sqrt(pi*u))./(2*kappa)).*(exp(-u/2)).*((1+u).*besseli(0,u./2) + u.*besseli(1,u./2));
% sEMMSE = max((HMMSE.*abs(outWFFT)).*exp(complex(0,angle(outWFFT))),0);
% % sEMMSE = (outWFFT - max(HMMSE.*abs(outWFFT),0)).*exp(complex(0,angle(outWFFT)));

% %% IFFT

% outWI(:,:,1) = real(ifft(sEstSS));                                                    
% outWI(:,:,2) = real(ifft(speech_est_Wiener));
% outWI(:,:,3) = real(ifft(sEMMSE));                                                    

% %% Overlap and Add

% orig = stoi(speech, noisySS(1:length(speech)), fs);
% % fprintf('STOI score unprocessed: %.4f\n', orig);
% figure(4)
% subplot(4,1,1)
% plot(noisySS)
% title(sprintf('Noisy Speech Signal with SNR = %.4g dB and STOI Score = %.4g', SNR, orig))
% ylim([-0.5 0.5])
% set(gca,'XTick',[]);

% for filt = 1:3
    % procSS = zeros(size(noisySS));                                              %Initializing the size of processed speech signal
    % procSS(1:wS+1) = outWI(:,1,filt);                                                %Getting back the first frame

    % for i = 1:numOfWins-1
        % procSS(i*wS/2:i*wS/2+wS) = outWI(:,i+1,filt) + procSS(i*wS/2:i*wS/2+wS);     %Overlapping and Adding the received frames to recover the processed speech in procSS
    % end

    % e = (speech - procSS(1:length(speech))).^2;
    % eNan = isnan(e);
    % eI = find(eNan == 0);
    % er = (sum(e(eI)));
    % ee = (e(eI));
    
    % error = sqrt(sum(e(eI)));
% %     fprintf('Sum of squared errors: %.4f\n', error);
    % intelligibility = stoi(speech, procSS(1:length(speech)), fs);
% %     fprintf('STOI score   processed: %.4f\n', intelligibility);
    % figure(4)
    % subplot(4,1,filt+1)
    % plot(procSS)
    % ylim([-0.5 0.5])
    % if filt < 4
      % set(gca,'XTick',[]);
    % end   
    % switch filt
        % case 1
            % title(sprintf('Enhanced Speech Signal using Spectral Subtraction with STOI Score = %.4g and SSE = %.4g', intelligibility, error));
        % case 2
            % title(sprintf('Enhanced Speech Signal using Wiener Filtering with STOI Score = %.4g and SSE = %.4g', intelligibility, error));
        % case 3
            % title(sprintf('Enhanced Speech Signal using MMSE Filtering with STOI Score = %.4g and SSE = %.4g', intelligibility, error));
    % end    
    % audiowrite(strcat('procSS', num2str(filt),'.wav'),procSS,fs)
% end
% figure()
% spectrogram(procSS,win,'yaxis')
% title('Spectrogram of Enhanced Speech Signal')


end