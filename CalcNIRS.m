function [ dHbR , dHbO, fig ] = CalcNIRS(dataFile, SDS, tissueType, plotChannelIdx, extinctionCoefficientsFile, DPFperTissueFile, relDPFfile )
%% Extract relevant data from the input dataFile
wavelengths =dataFile.SD.Lambda; %Two wavelengths [nm]
time = dataFile.t; %Time vector
intensities = dataFile.d;
%intensityLow = dataFile.d(:,1:20); %intensity levels at low WL
%intensityHigh = dataFile.d(:,21:end); %intensity levels at high WL

%% Set default values 
if nargin < 7 || isempty(relDPFfile)
    relDPFfile = '.\RelativeDPFCoefficients.csv';
end
if nargin < 6 || isempty(DPFperTissueFile)
    DPFperTissueFile = '.\DPFperTissue.txt';
end
if nargin < 5 || isempty(extinctionCoefficientsFile)
    extinctionCoefficientsFile = '.\ExtinctionCoefficientsData.csv';
end
if nargin < 4 || isempty(plotChannelIdx)
    plotChannelIdx = [];
end

%% Load extinction coefficients and DPF data
extinctionCoefficients = readtable(extinctionCoefficientsFile);
DPFperTissue = readtable(DPFperTissueFile);
relDPF = readtable(relDPFfile);


%% Calculate DPF for the given tissue type and wavelengths
DPF_807nm = DPFperTissue.DPF(strcmp(DPFperTissue.Tissue, tissueType));
relDPF_factors = interp1(relDPF.wavelength, relDPF{:, 2:end}, wavelengths);
DPF = DPF_807nm .* relDPF_factors;

%% Calculate optical densities
I0 = intensities(1, :);
OD = log10(I0 ./ intensities);

%% Calculate extinction coefficients for the given wavelengths
epsilon_HbR = interp1(extinctionCoefficients.wavelength, extinctionCoefficients.HHb, wavelengths);
epsilon_HbO = interp1(extinctionCoefficients.wavelength, extinctionCoefficients.HbO2, wavelengths);

%% Output
dHbR = OD(:,1:20) ./ (epsilon_HbR(1) * DPF(1) * SDS);
dHbO = OD(:,20:end) ./ (epsilon_HbO (2)* DPF(2) * SDS);

%% Plot the specified channels
figs = [];
if ~isempty(plotChannelIdx) && isvector(plotChannelIdx) && all(plotChannelIdx >= 1 & plotChannelIdx <= 20 & rem(plotChannelIdx, 1) == 0)
    for ch = plotChannelIdx
        fig = figure;
        plot(time, dHbR(:, ch), 'b', 'LineWidth', 2);
        hold on;
        plot(time, dHbO(:, ch), 'r', 'LineWidth', 2);
        hold off;
        title(sprintf('Channel %d', ch));
        xlabel('Time (s)');
        ylabel('Concentration Change');
        legend('dHbR', 'dHbO', 'Location', 'best');
        figs = [figs, fig]; % Store the figure handle
    end
elseif isempty(plotChannelIdx)
    disp('No channels specified for plotting.');
else
    disp('Invalid input for plotChannelIdx. Please provide a vector with values in the range [1-20].');
end


