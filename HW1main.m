clc;clear;close all;
% Loading the data
data1 = load('FN_032_V1_Postdose1_Nback.mat');
data2 = load('FN_032_V1_Postdose1_Nback.mat');

SDS = 3; %Sourse-Detector Separation distance [cm]
relDPF =  "RelativeDPFCoefficients.csv"; 
DPFperTissue  = "DPFperTissue.txt"; %relative DPF according to wavelength
extinctionCoefficients = "ExtinctionCoefficientsData.csv"; %.csv file with the following columns : wavelength, Water, HbO2, HHb, FatSoybean
tissueType = 'adult_forearm'; %Options: 'adult_forearm' \ 'baby_head' \ 'adult_head' \ 'adult_leg' 
plotChannelIdx = [1,2];   %vector indicating channels to plot. 
[ dHbR , dHbO, fig ] = CalcNIRS(data1, SDS, tissueType, plotChannelIdx, extinctionCoefficients , DPFperTissue, relDPF );


