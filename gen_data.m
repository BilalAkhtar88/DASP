function X = gen_data(S,M,N,dist,theta,speed,numFBin,fs)
%S is the signals from two angles, it is (numel(theta) x N) matrix
%M is the number of antennas in array
%N is the number of samples in S
%Delta is spacing between the antennas in meters, to be converted to wavelengths
%Assuming channel response is impulse for both the channels and that the
%channels are noise free
%Assuming Noise is added to already added to S and noise free channel
%X is after receiving the signal from beamforming antennas

fs
fq = ((0:numFBin-1)*fs/numFBin);
fq(numel(fq))
Delta = (dist.*((0:numFBin-1)*fs/numFBin))./(speed);


for t = 1:numel(theta)
    a(:,((t-1)*length(Delta)+1):t*length(Delta)) = gen_a(M,Delta,theta(t));
end

U = abs(ctranspose(a)*a);
% ssA = size(a)
% ssS = size(S)
X = a*S;
% ssX = size(a)


end