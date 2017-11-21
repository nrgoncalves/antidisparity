clc; close all; clear all;
addpath(genpath([pwd filesep 'functions']))

% directory where results are saved
resultsDir = fullfile(pwd,'data', 'signed-corr-mws-time', 'subjects');
subjects = {'AW', 'NG', 'LS'}; 
expId = 1;  % for panel b and c
% expId = 3;  % for panel d (control)

% uncomment the following two lines to plot results obtained for
% experiment with longer stimulus duration 
% subjects = {'NG'};
% expId = 4;

if numel(subjects) == 1
    figure('pos', [100,100,250,200],...
           'DefaultAxesFontName', 'Helvetica',...
            'DefaultAxesFontSize', 8);
else
    figure('pos', [100,100,800,200],...
           'DefaultAxesFontName', 'Helvetica',...
           'DefaultAxesFontSize', 8);
end

for ss = 1:numel(subjects)
    
    subplot(1,numel(subjects),ss);
    hold on;
    for ee = expId
        
        runFiles = wildcardsearch(fullfile(resultsDir, subjects{ss}, 'results'),...
            ['exp' num2str(ee)]);
        
        if ~isempty(runFiles)
            plotResultsTime(subjects{ss}, ee, runFiles);
        end
        
    end
    
end

if expId == 1
    plot2svg('figures/Fig6b.svg')
elseif expId == 4
    plot2svg('figures/Fig6_longerPresentationTime.svg')
else
    plot2svg('figures/Fig6_control.svg')
end