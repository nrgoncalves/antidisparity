addpath(genpath([pwd filesep 'functions']))

resultsLabel = 'prestim_sem';    

subjects = {'AW', 'NG', 'RR', 'LS', 'EM'};
baseDir = [pwd '\data\eye-tracking-data'];

% calibration design
design.xcoords = [2 0 -2 0 2 0 -4 0]';        % Horizontal target positions (deg)
design.ycoords = [0 2 0 -2 0 -2 0 4]';        % Vertical target positions (deg)

design.onsetTimes = [3, 7, 11, 15, 19, 23, 27, 31];    % seconds
design.offsetTimes = [5, 9, 13, 17, 21, 25, 29, 33];   % seconds

for s = 1:numel(subjects)
    
    eyetrackingFiles = {[baseDir '/' subjects{s} '/exp2_run1.mat'],...
                        [baseDir '/' subjects{s} '/exp2_run2.mat']};

    eyetrackingCalFiles = {[baseDir '/' subjects{s} '/exp2_run1_cal.mat'],...
                           [baseDir '/' subjects{s} '/exp2_run2_cal.mat']};
    
    D = EyeTrackingDataset(eyetrackingFiles, ...
                           eyetrackingCalFiles);
    
    
	% set calibration design
	D.setCalibrationDesign(design);
                       
    % estimate calibration parameters
    D.estimateCalibrationParams();
    
    % apply calibration
    D.applyCalibrationParams();
    
    % pause;
    close all;
    
    % Epoch data
    D.epochAround('Stim1', [-100, 400], 'Same');
    D.epochAround('Stim2', [-100, 400], 'Opposite');
        
    
    % Group trials per condition
    D.groupByCondition();
    
    
    % Reject outliers
    D.rejectOutliers();
    
        
    % Compute vergence
    D.computeCondVergence({'Same', 'Opposite'}, 'prestim')
    ylim([-1, 1])
    plot2svg(['eyetracking/' subjects{s} '_vergence_' resultsLabel '.svg'], gcf);
    
    % D.computeSaccadeStatistics({'Same', 'Opposite'})
    
    % Plot heatmap of eye position
    bins = {linspace(-5,5,51), linspace(-5,5,51)};
    map = D.plotHeatmap('Same', 'L', bins);
    
    % pause;
    close all;
    save(['eyetracking/' subjects{s} '_' resultsLabel '.mat'], 'D');
    clear D;               
end
