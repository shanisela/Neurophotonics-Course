clc;clear;close all;
% Loading the intensity meadata
dataFileList = {'FN_031_V2_Postdose2_Nback.mat','FN_032_V1_Postdose1_Nback.mat'};
SDS = 3; %Sourse-Detector Separation distance [cm]
relDPF =  "RelativeDPFCoefficients.csv"; 
DPFperTissue  = "DPFperTissue.txt"; %Relative DPF according to wavelength
extinctionCoefficients = "ExtinctionCoefficientsData.csv"; %.csv file with the following columns : wavelength, Water, HbO2, HHb, FatSoybean
tissueType = 'adult_head'; %Options: 'adult_forearm' \ 'baby_head' \ 'adult_head' \ 'adult_leg' 
plotChannelIdx = [1,2];   %Vector indicating channels to plot. 

for i = 1:length(dataFileList)
    [ dHbR , dHbO, fig ] = CalcNIRS(dataFileList{i}, SDS, tissueType, plotChannelIdx, extinctionCoefficients , DPFperTissue, relDPF );
    
    % Set unique figure identifier to the file used

end 


%% Question 2
close all;
FirstFile = load(dataFileList{1});
FirstChannel = FirstFile.d(:,1);
L = length(FirstChannel);
Fs = 7.8; % Sampling frequency in Hz
Y = abs(fft(FirstChannel)/L);% Normalized Fourier Transform
f = (0:L-1) * (Fs/L); % Frequency axis

% Only take the first half of the spectrum
Y = Y(1:floor(L/2));
f = f(1:floor(L/2));

% Define frequency range for signal and noise
noiseFrequencyStart = 2.5; % Noise frequency starts at 2.5 Hz
HeartBeatFreq = [1 2];
fMaxVal = 3.9;

% Find the indices for the relevant frequency ranges
pulseIndices = find(f >= HeartBeatFreq(1) & f <= HeartBeatFreq(2));
noiseIndices = find(f >= noiseFrequencyStart & f <= fMaxVal);

% Calculate the signal strength at the pulse frequency
[~, HRpulsePeakIDX] = max(Y(pulseIndices));
HRpulsePeakfft = Y(pulseIndices(HRpulsePeakIDX));

% Calculate the noise as the mean FFT values in the noise range
noiseStrength = mean(Y(noiseIndices));

% Calculate the Signal-to-Noise Ratio (SNR)
SNR = HRpulsePeakfft / noiseStrength;

% Plot the Fourier Transform
figure;
semilogy(f, Y.');
hold on
plot(f(pulseIndices(HRpulsePeakIDX)),HRpulsePeakfft,'o')
hold off
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Fourier Transform of First Channel');
legend('fft signal', 'HR peak')
% Display the SNR
fprintf('The Signal-to-Noise Ratio (SNR) is: %.2f\n', SNR);