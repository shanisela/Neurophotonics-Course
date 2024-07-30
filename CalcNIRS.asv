function [ dHbR , dHbO, fig ] = CalcNIRS(dataFile, SDS, tissueType, plotChannelIdx, extinctionCoefficientsFile, DPFperTissueFile, relDPFfile )
%% Verifying input validity
% Check if dataFile is a .mat file and exists
if ~ischar(dataFile) || ~exist(dataFile, 'file') || ~strcmp(dataFile(end-3:end), '.mat')
    error('dataFile must be a valid .mat file.');
end
% Check if SDS is a numeric value
if ~isnumeric(SDS) || ~isscalar(SDS)
    error('SDS must be a numeric scalar.');
end
% Check if tissueType is a string
if ~ischar(tissueType)
    error('tissueType must be a string.');
end
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

%% Extract relevant data from the input dataFile
% Load data from .mat file and validate its contents
data = load(dataFile); 
if ~isfield(data, 'SD') || ~isfield(data.SD, 'Lambda') || ~isfield(data, 't') || ~isfield(data, 'd')
    error('dataFile must contain the fields: DS.Lambda, t, and d.');
end
if size(data.d, 2) ~= 40
    error('Intensity data (d) must have 40 rows.');
end

wavelengths =data.SD.Lambda; %Two wavelengths [nm]
time = data.t; %Time vector [sec]
intensities = data.d; %Intensity measurements 

%% Load extinction coefficients and DPF data
extinctionCoefficients = readtable(extinctionCoefficientsFile); %Absorption coefficient (epsilon) [L/mol*cm]
DPFperTissue = readtable(DPFperTissueFile);
relDPF = readtable(relDPFfile);


%% Calculate DPF for the given tissue type and wavelengths
DPF807nm = DPFperTissue.DPF(strcmp(DPFperTissue.Tissue, tissueType)); %Finding the DPF at 807 for the required tissue type
relDPFCoeff = interp1(relDPF.wavelength, relDPF.relDPFcoeff, wavelengths); %Find the relative DPF for the required wl
DPF = DPF807nm .* relDPFCoeff; %calculate final DPF

%% Calculate optical densities
I0 = intensities(1,:);
OD = log10(I0 ./ intensities);

%% Calculate extinction coefficients for the given wavelengths
epsilonHbR = interp1(extinctionCoefficients.wavelength, extinctionCoefficients.HHb, wavelengths);
epsilonHbO = interp1(extinctionCoefficients.wavelength, extinctionCoefficients.HbO2, wavelengths);

%% Output
dHbR = OD(:,1:20) ./ (epsilonHbR(1) * DPF(1) * SDS);
dHbO = OD(:,21:end) ./ (epsilonHbO (2)* DPF(2) * SDS);

%% Plot the specified channels
figs = [];
if ~isempty(plotChannelIdx) && isvector(plotChannelIdx) 
    for ch = plotChannelIdx
        fig = figure;
        plot(time, dHbR(:, ch), 'b');
        hold on;
        plot(time, dHbO(:, ch), 'r');
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


