clc;clear;close all;
% Loading the intensity meadata
dataFileList = {'FN_032_V1_Postdose1_Nback.mat','FN_031_V2_Postdose2_Nback.mat'};

SDS = 3; %Sourse-Detector Separation distance [cm]
relDPF =  "RelativeDPFCoefficients.csv"; 
DPFperTissue  = "DPFperTissue.txt"; %Relative DPF according to wavelength
extinctionCoefficients = "ExtinctionCoefficientsData.csv"; %.csv file with the following columns : wavelength, Water, HbO2, HHb, FatSoybean
tissueType = 'adult_head'; %Options: 'adult_forearm' \ 'baby_head' \ 'adult_head' \ 'adult_leg' 
plotChannelIdx = [1,2];   %Vector indicating channels to plot. 
for i = length(dataFileList)
    [ dHbR , dHbO, fig ] = CalcNIRS(dataFileList{i}, SDS, tissueType, plotChannelIdx, extinctionCoefficients , DPFperTissue, relDPF );
    
end 


%% Question 2

firstChannel = load(dataFile1).d (:,1);
L = length(firstChannel);
Fs = 7.8; % Sampling frequency in Hz
% Fourier Transform
Y = fft(firstChannel) / L;
% Frequency axis
f = (0:L-1) * Fs / L; % Frequency range
% Only take the first half of the spectrum
%Y = Y(1:floor(L/2));
%f = f(1:floor(L/2));
% Plot the Fourier Transform
figure;
plot(f, abs(Y));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Fourier Transform of First Channel');