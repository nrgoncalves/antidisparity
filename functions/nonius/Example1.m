% Example 1: analyze calibration, validation and task data.
% 
% In this script we assume that we have collected two sets of measurements 
% for which we know the "ground truth" positions, and one set of 
% measurements during our actual experiment. We will use them
% to calibrate, validate and analyze task data. 

clear all;
close all;

eyetrackingFile = 'sample-data/ex1/NG.mat';
eyetrackingCalFile = 'sample-data/ex1/Cal1.mat';
eyetrackingValFile = 'sample-data/ex1/Cal2.mat';
D = EyeTrackingDataset(eyetrackingFile, ...
                       eyetrackingCalFile, ...
                       eyetrackingValFile);


% Calibration
% - calibration design
design.xcoords = [2 0 -2 0 2 0 -4 0]';        % Horizontal target positions (deg)
design.ycoords = [0 2 0 -2 0 -2 0 4]';        % Vertical target positions (deg)

design.onsetTimes = [3, 7, 11, 15, 19, 23, 27, 31];    % seconds
design.offsetTimes = [5, 9, 13, 17, 21, 25, 29, 33];   % seconds

D.setCalibrationDesign(design);

% - estimate calibration parameters
D.estimateCalibrationParams();

% - validate calibration parameters
D.setValidationDesign(design);
D.validateCalibrationParams();

% - apply calibration
D.applyCalibrationParams();


% Epoch data
D.epochAround('Stim1', [-100, 500], 'Same');
D.epochAround('Stim2', [-100, 500], 'Opposite');

% Plot epoched data (e.g.)
% D.plotEpoch('Same', 'R', 'H')
% D.plotEpoch('Opposite', 'L', 'V')


% Group trials per condition
D.groupByCondition();


% Reject outliers
D.rejectOutliers();


% Plot comparison between conditions
D.plotConditions({'Same', 'Opposite'}, 'L', 'H');


% Compute vergence
D.computeCondVergence({'Same', 'Opposite'})


% Plot heatmap of eye position
bins = {linspace(-5,5,51), linspace(-5,5,51)};
map = D.plotHeatmap('Same', 'L', bins);