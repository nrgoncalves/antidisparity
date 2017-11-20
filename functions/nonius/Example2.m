% Example 2: analyze calibration and validation data.
% 
% In this script we assume that we have collected two sets of measurements 
% for which we know the "ground truth" positions. We will use one of them
% to calibrate the measurements, and the other one for validation purposes.

clear all;
close all;

eyetrackingFile = 'sample-data/ex2/Cal1.mat';
eyetrackingCalFile = 'sample-data/ex2/Cal1.mat';
eyetrackingValFile = 'sample-data/ex2/Cal3.mat';
D = EyeTrackingDataset(eyetrackingFile, ...
                       eyetrackingCalFile, ...
                       eyetrackingValFile);

% Calibration
% - calibration design
[Xoffset, Yoffset] = meshgrid(-4:2:4, -4:2:4);
design.xcoords = Xoffset(:);
design.ycoords = Yoffset(:);

design.onsetTimes = 3:4:(numel(design.xcoords)*4);
design.offsetTimes = design.onsetTimes + 2;

D.setCalibrationDesign(design);

% - estimate calibration parameters
D.estimateCalibrationParams();

% Validation
D.setValidationDesign(design);
D.validateCalibrationParams();