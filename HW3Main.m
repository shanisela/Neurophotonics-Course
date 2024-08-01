clc; close all; clear; 
% Neurophotonics HW3

%% part 1+2
% 1. Ask User For Recording Paths per recording type

% Select the folder containing the short recording
%With the laser and trigger ON, blackLevel =0 , exposureTime=15ms short- 30sec 
shortFolderPath = uigetdir('', 'Select the Folder Containing short Recordings');
if shortFolderPath == 0
    disp('User canceled the folder selection.');
    return;
else
    disp(['User selected folder: ', shortFolderPath]);
end

% Select the folder containing the Read Noise recording
%(no laser, cover the input, minimum exposure time(0.021ms), same game, blackLevel 20DU – 200 frames
readNoiseFolderPath = uigetdir('', 'Select the Folder Containing Read Noise Recordings');
if readNoiseFolderPath == 0
    disp('User canceled the folder selection.');
    return;
else
    disp(['User selected folder: ', readNoiseFolderPath]);
end

% Select the folder containing the Dark Background recording
%(no laser, blackLevel =0, exposureTime=15ms, everything else is the same – 400 frames
backgroundFolderPath = uigetdir('', 'Select the Folder Containing Background Recordings');
if backgroundFolderPath == 0
    disp('User canceled the folder selection.');
    return;
else
    disp(['User selected folder: ', backgroundFolderPath]);
end

% Select the folder containing the Long recording
%With the laser and trigger ON, blackLevel =0 , exposureTime=15ms 
% long- breath holding x 4 times , 30 seconds each time
longFolderPath = uigetdir('', 'Select the Folder Containing long Recordings');
if longFolderPath == 0
    disp('User canceled the folder selection.');
    return;
else
    disp(['User selected folder: ', longFolderPath]);
end


% Get a list of TIFF files in the selected folders
shortFilePattern = fullfile(shortFolderPath, '*.tiff');
readNoiseFilePattern = fullfile(readNoiseFolderPath, '*.tiff');
backgroundFilePattern = fullfile(backgroundFolderPath, '*.tiff');
longFilePattern = fullfile(longFolderPath, '*.tiff');

shortTiffFiles = dir(shortFilePattern);
readNoiseTiffFiles = dir(readNoiseFilePattern);
backgroundTiffFiles = dir(backgroundFilePattern);
longTiffFiles = dir(longFilePattern);


% 2. Determine an ROI (plot one frame then ask the user to draw a circle in it)
firstFileName = shortTiffFiles(1).name;
firstFilePath = fullfile(shortFolderPath, firstFileName);
info = imfinfo(firstFilePath);
frame1 = imread(firstFilePath, 1)./16;

% Display the first frame
figure;
imshow(frame1, []);
title('Draw a circular ROI and double-click to confirm');
h = imellipse(gca);
position = wait(h); % Wait for the user to double-click to finalize the position
mask = createMask(h);
close;

disp('Finished part 1+2')
%% Part 3 
clc;
% 3.Calculate Read Noise Matrix (read noise in each window)
%opening readnois nois frames
numFramesReadNoise = length(readNoiseTiffFiles);
readNoiseFrames = zeros([size(frame1), numFramesReadNoise]);

for i = 1:numFramesReadNoise
    readNoiseFileName = readNoiseTiffFiles(i).name;
    readNoiseFilePath = fullfile(readNoiseFolderPath, readNoiseFileName);
    readNoiseFrames(:,:,i) = imread(readNoiseFilePath)./16;
end

readNois = std(readNoiseFrames,0,3);

% Set window size
window = 7; 

% Apply the filter to the input matrix A
reasNoisMatrix = imboxfilt(readNois, window);

disp('Finished part 3')

%% part 4 + 5
clc;
%4. Calculate DarkBackground Image
%opening background nois frames
numFramesBackgroungNois = length(backgroundTiffFiles);
backgroundNoiseFrames = zeros([size(frame1), numFramesBackgroungNois]);

for i = 1:numFramesBackgroungNois
    backgroundFileName = backgroundTiffFiles(i).name;
    backgroundNoiseFilePath = fullfile(backgroundFolderPath, backgroundFileName);
    backgroundNoiseFrames(:,:,i) = imread(backgroundNoiseFilePath)./16;
end
%DarkBackground Image
backgroundNoisMatrix = mean(backgroundNoiseFrames,3);


% 5. Calculate Pixels-Non-Uniformity (σ_sp)
% Opening short recording frames
numFramesShortRecording = 500; 
shortFrames = zeros([size(frame1), numFramesShortRecording]);
for i = 1:numFramesShortRecording
    laserFileName = shortTiffFiles(i).name; % Adjust if necessary
    laserFilePath = fullfile(shortFolderPath, laserFileName);
    shortFrames(:,:,i) = double(imread(laserFilePath))./16;
end

%Opening long recording frames
numFramesLongRecording = 500; 
longFrames500 = zeros([size(frame1), numFramesLongRecording]);
for i = 1:numFramesLongRecording
    laserFileName = longTiffFiles(i).name; % Adjust if necessary
    laserFilePath = fullfile(longFolderPath, laserFileName);
    longFrames500(:,:,i) = double(imread(laserFilePath))./16;
end


shortAverage = mean(shortFrames,3);
longAverage = mean(longFrames500,3);

minusBackgroundShort = shortAverage - backgroundNoisMatrix;
minusBackgroundLong = longAverage - backgroundNoisMatrix;

nonUnifShort = stdfilt(minusBackgroundShort,true(window)).^2; %calc variance in 7X7 window
nonUnifLong = stdfilt(minusBackgroundLong,true(window)).^2; %calc variance in 7X7 window

disp('Finished part 4+5')

%% Part 6
clc;
% 6. Calculate G[DU/e]
gain = 24; %[dB]
maxCapacity = 10500;  %[e]
nBits = 12; %bits

Gin = 10^(gain/20);
Gbase = (2^nBits)/maxCapacity;
G = Gbase * Gin; % [DU/e]

disp('Finished part 6')

%% Part 7
clc;
% 7. For every frame, subtract the background, calc contrast per window
KShortList = [];
numFramesShortRecording = length(shortTiffFiles);
for i = 1:numFramesShortRecording
    laserFileName = shortTiffFiles(i).name; % Adjust if necessary
    laserFilePath = fullfile(shortFolderPath, laserFileName);
    frame = double(imread(laserFilePath))./16;
    %Subtract background from the frame
    noBackground = frame - backgroundNoisMatrix;
    %Raw contrast per window
    kRawPerWindow = (stdfilt(noBackground,true(window)) ./ imboxfilt(noBackground,window) ).^2;
    % Fixed contrast per window
    Ks = G ./ imboxfilt(noBackground, window);
    Kr = (reasNoisMatrix ./ imboxfilt(noBackground, window)).^2;
    Ksp = (nonUnifShort ./ imboxfilt(noBackground, window)).^2;
    Kq = 1/12 .* ((imboxfilt(noBackground, window)).^2);
    KfShort = kRawPerWindow - Ks- Kr - Ksp - Kq ;
    KinROI = mean2(KfShort .* mask);
    KShortList = [KShortList,KinROI];
end

KLongList = [];
numFramesLongRecording = length(longTiffFiles);
for i = 1:numFramesLongRecording
    laserFileName = longTiffFiles(i).name; % Adjust if necessary
    laserFilePath = fullfile(longFolderPath, laserFileName);
    frame = double(imread(laserFilePath))./16;
    %Subtract background from the frame
    noBackground = frame - backgroundNoisMatrix;
    %Raw contrast per window
    kRawPerWindow = (stdfilt(noBackground,true(window)) ./ imboxfilt(noBackground,window) ).^2;
    % Fixed contrast per window
    Ks = G ./ imboxfilt(noBackground, window);  %Shot nois
    Kr = (reasNoisMatrix ./ imboxfilt(noBackground, window)).^2; %Readout nois 
    Ksp = (nonUnifLong ./ imboxfilt(noBackground, window)).^2; %Spatial non-uniformity
    Kq = 1/12 .* ((imboxfilt(noBackground, window)).^2);  %Quantisation 
    KfLong = kRawPerWindow - Ks- Kr - Ksp - Kq ;
    KinROI = mean2(KfLong .* mask);
    KLongList = [KLongList,KinROI];
end

disp('Finished part 7')


%% Part 8+9
clc;
%Plot <Kf> vs time for short and long recordings
numShortFrames = size(shortTiffFiles);
numLongFrames = size(longTiffFiles);
samplingRate = 20; % in Hz
dt = 1 / samplingRate; % Time between frames in seconds

% Create the time vector
timeVectorShort = (0:numShortFrames-1) * dt;
timeVectorLong = (0:numLongFrames-1) * dt;

figure; 
plot(timeVectorShort,KShortList);  
xlabel('Time[sec]')
ylabel('<K_f>')
title('Short Recording');

figure;
plot(timeVectorLong,KLongList)
xlabel('Time[sec]')
ylabel('<K_f>')
disp('Finished part 8')
title('Long Recording');

%Adding marks of breath hold start and end for the long recording

taskStart = [30, 90, 150, 210]; % Start times in seconds (adjust according to your data)
taskDuration = 30; % Duration of each breath hold in seconds
taskColor = [0.8, 0.8, 0.8]; % Light grey color for the patch

% Get the current Y-axis limits
ylims = get(gca, 'YLim');

% Loop through each breath hold event and add the markup
for iter_i = 1:numel(taskStart)
    patch([taskStart(iter_i) taskStart(iter_i) + taskDuration, taskStart(iter_i) + taskDuration taskStart(iter_i)], ...
          [ylims(1) ylims(1) ylims(2) ylims(2)], taskColor, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
end

% Ensure the Y-axis limits remain unchanged
set(gca, 'YLim', ylims);

disp('Finished part 9')
