clc; clear all;
addpath(genpath([pwd filesep 'functions']))
pathToPsignifit = 'D:\NRG\tools\psignifit';
addpath(genpath(pathToPsignifit))

% change this to 1 if we want to fit a psychometric function using 
% psignifit. pathToPsignifit must point to the psignifit toolbox
fitPsychFlag = 1; 

% directory where results are saved
resultsDir = fullfile(pwd,'data', 'signed-corr-mws-levels', 'subjects');
subjects = {'AW', 'NG', 'LS'}; 

% experiment identifiers
expId = 1;
pltColor = [0.7, 0.7, 0.7];

figure('pos', [100,100,800,200],...
       'DefaultAxesFontName', 'Helvetica',...
       'DefaultAxesFontSize', 8);

for ss = 1:numel(subjects)
    subplot(1,numel(subjects),ss);
    hold on;
    for ee = expId
        
        runFiles = wildcardsearch(fullfile(resultsDir, subjects{ss}, 'results'),...
            ['exp' num2str(ee)]);
        
        if ~isempty(runFiles)
            plotResultsLevelsBootstrapped(subjects{ss}, ee, runFiles, pltColor, fitPsychFlag, 'abs', 1);
        end
        
    end
            
end

plot2svg('figures/Fig3b.svg')