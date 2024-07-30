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
backgroungdFilePattern = fullfile(backgroundFolderPath, '*.tiff');
longFilePattern = fullfile(longFolderPath, '*.tiff');

shortTiffFiles = dir(shortFilePattern);
readNoiseTiffFiles = dir(readNoiseFilePattern);
backgroundTiffFiles = dir(backgroungdFilePattern);
longTiffFiles = dir(longFilePattern);


% 2. Determine an ROI (plot one frame then ask the user to draw a circle in it)
firstFileName = shortTiffFiles(1).name;
firstFilePath = fullfile(shortFolderPath, firstFileName);
info = imfinfo(firstFilePath);
frame1 = imread(firstFilePath, 1);

% Display the first frame
figure;
imshow(frame1, []);
title('Draw a circular ROI and double-click to confirm');
h = imellipse(gca);
position = wait(h); % Wait for the user to double-click to finalize the position
mask = createMask(h);
close;
%% Part 3 
% 3.Calculate Read Noise Matrix (read noise in each window)
numFramesReadNoise = length(readNoiseTiffFiles); % Number of frames in Read Noise recordings
readNoiseFrames = zeros([size(frame1), numFramesReadNoise]);

% Read the frames from the Read Noise recordings and apply the mask
maskedReadNoiseFrames = [];
for i = 1:numFramesReadNoise
    readNoiseFileName = readNoiseTiffFiles(i).name;
    readNoiseFilePath = fullfile(readNoiseFolderPath, readNoiseFileName);
    currentFrame = imread(readNoiseFilePath);
    % Cast the mask to the same type as currentFrame
    maskCasted = cast(mask, class(currentFrame));
    % Apply the mask to the current frame
    maskedFrame = currentFrame .* maskCasted; % Apply mask
    maskedReadNoiseFrames = cat(3, maskedReadNoiseFrames, maskedFrame);
  
end

% Calculate the mean of the masked Read Noise frames
meanReadNoise = mean(maskedReadNoiseFrames, 3);
imshow(meanReadNoise)


%% part 4
% 4. Calculate Pixels-Non-Uniformity (σ_sp)
numFramesShortRecording = 500; % Number of frames with laser (adjust as necessary)
laserFrames = zeros([size(frame1), numFramesShortRecording]);

% Assuming the laser recordings start from the appropriate frame
for i = 1:numFramesShortRecording
    laserFileName = mainTiffFiles(400 + i).name; % Adjust if necessary
    laserFilePath = fullfile(mainFolderPath, laserFileName);
    laserFrames(:,:,i) = imread(laserFilePath);
end

% Average the frames
meanLaserFrames = mean(laserFrames, 3);

% Load Dark Background frame (assuming the first 400 frames are Dark Background)
numFramesDarkBackground = 400;
darkBackgroundFrames = zeros([size(frame1), numFramesDarkBackground]);
for i = 1:numFramesDarkBackground
    darkBackgroundFileName = mainTiffFiles(i).name; % Adjust if necessary
    darkBackgroundFilePath = fullfile(mainFolderPath, darkBackgroundFileName);
    darkBackgroundFrames(:,:,i) = imread(darkBackgroundFilePath);
end
meanDarkBackground = mean(darkBackgroundFrames, 3);

% Subtract the background
correctedFrames = meanLaserFrames - meanDarkBackground;

% Calculate variance in each window
varianceMatrix = stdfilt(correctedFrames, ones(3)).^2; % Standard deviation filter

%% Part 5
% 5. Calculate G[DU/e]
% Placeholder: Replace with actual calculation depending on your setup
G = 1; % Dummy value

% 6. For every frame, subtract the background
for i = 1:numFramesShortRecording
    laserFrames(:,:,i) = laserFrames(:,:,i) - meanDarkBackground;
end

% Display results
disp('Read Noise Matrix:');
disp(readNoiseMatrix);
disp('Variance Matrix:');
disp(varianceMatrix);
disp('G [DU/e]:');
disp(G);