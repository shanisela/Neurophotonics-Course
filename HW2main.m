clc; clear; close all;
% Neurophotonics HW 2
%Reading the recoedings
addpath("auxiliary_code\")
addpath("records\")

numOfRecords = 4;
path = 'C:\Users\shani\Documents\Masters Degree\Neurophotonics\HW2\records';
recordingsNames = dir(path);
noiseCell = repmat(struct('RecordingName',[],'GlobalTemporalNoise',[], 'GlobalSpatialNoise', [],...
    'LocalSpatialNoise_AvginTime', [],...
    'LocalSpatialNoise_AvgFrames',[]),1,numOfRecords);
% Filter out '.' and '..' and keep only folders
recordingsNames = recordingsNames([recordingsNames.isdir]);
recordingsNames = recordingsNames(~ismember({recordingsNames.name}, {'.', '..'}));
window = 7;
maxCapacity = 10500;  %[e]
nBits = 12; %bits

for i = 1:numOfRecords
    subFolderPath = fullfile(path,recordingsNames(i).name);
    noiseCell(i).RecordingName= recordingsNames(i).name;
    [rec, info] = ReadRecord(subFolderPath);
    %Calculating Temporal Noise
    noiseCell(i).GlobalTemporalNois = mean2(std(rec,0,3)); 
    %Global Spatial Noise (after averaging in time)
    noiseCell(i).GlobalSpatialNoise = std(mean(rec,3),0,"all"); 
    %Local Spatial Noise with window size of 7 (after averaging in time)
    noiseCell(i).LocalSpatialNoise_AvgInTime = mean2(stdfilt(mean(rec,3), true(window))); 
    %Local Spatial Noise with window size of 7 per frame – then average over all frames.
    numFrames = size(rec, 3);
    localNoisePerFrame = zeros(numFrames, 1);
    for k = 1:numFrames
        % Extract the current frame
        currentFrame = rec(:, :, k);
             
        % Compute standard deviation filter for the current frame
        stdDevFilteredFrame = stdfilt(currentFrame, true(7));
        
        % Compute mean of the filtered standard deviation
        localNoisePerFrame(k) = mean(stdDevFilteredFrame(:));
    end
    noiseCell(i).LocalSpatialNoise_AvgFrames = mean2(localNoisePerFrame);
    disp(['Processing subfolder: ', noiseCell(i).RecordingName ]);
    
    totalNoise = ((noiseCell(i).GlobalTemporalNois).^2+(noiseCell(i).LocalSpatialNoise_AvgInTime).^2)^0.5;
    fprintf('Total noise is %.2f\n', totalNoise)
% Answer to question 5 
    if noiseCell(i).RecordingName == "WhitePaper_Gain24dB_expT0.5ms_BlackLevel0DU"
        q5ans = ((noiseCell(i).GlobalTemporalNois)^2)/ (mean2(rec));
        compValue = ((2^nBits)/maxCapacity)*(10^(info.name.Gain/20));
        disp ([q5ans,compValue])
        imagesc(mean(rec,3))
        colorbar;
    end 
    
end 



