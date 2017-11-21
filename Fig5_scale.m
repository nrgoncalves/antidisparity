clc; clear all;
addpath(genpath([pwd filesep 'functions']))

% directory where results are saved
resultsDir = fullfile(pwd, 'data', 'signed-corr-mws-scale', 'subjects');
subjects = {'AW', 'NG', 'LS'}; 

% experiment identifiers
expId = 1;

figure('pos', [100,100,800,200],...
       'DefaultAxesFontName', 'Helvetica',...
       'DefaultAxesFontSize', 8);

for ss = 1:numel(subjects)
        subplot(1,numel(subjects), ss);
        hold on;
        for ee = expId
            
            runFiles = wildcardsearch(fullfile(resultsDir, subjects{ss}, 'results'),...
                                               ['exp' num2str(ee)]);
            
            if ~isempty(runFiles)
                plotResultsScale(subjects{ss}, ee, runFiles, 0);
            end
            
        end
        
end

plot2svg('figures/Fig5b.svg')