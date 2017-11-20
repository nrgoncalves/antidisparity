clc; clear all;
addpath(genpath([pwd filesep 'functions']))
fitPsychFlag = 0; 

% directory where results are saved
resultsDir = [pwd '\data\signed-corr-mws-levels\subjects'];
subjects = {'AW', 'NG', 'LS'}; 

% experiment identifiers
expId = [4 1 2 6];
pltColor = [0.7, 0.7, 0.7];

for ss = 1:numel(subjects)
    subplot(1,numel(subjects),ss);
    hold on;
    for ee = expId
        
        runFiles = wildcardsearch(fullfile(resultsDir, subjects{ss}, 'results'),...
            ['exp' num2str(ee)]);
        
        if ~isempty(runFiles)
            plotResultsLevelsBootstrapped(subjects{ss}, ee, runFiles, pltColor, fitPsychFlag, 'none', 0);
        end
        
    end
            
end

left = [];
right = [];

for ss = 1:numel(subjects)
    for ee = expId
        fname = ['lr-correlation/' subjects{ss} '_exp' num2str(ee) 'LevelsTuningCurveLeftRightCorrelation.mat'];
        if exist(fname,'file')==2
            res = load(fname);
            left = [left; res.left];
            right = [right; res.right];
            clear res;
        end
    end        
end    

[r, p] = corr(left, right);
figure;
plot(left, right, 'o')