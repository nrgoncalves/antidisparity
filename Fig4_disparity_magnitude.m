clc; clear all;
addpath(genpath([pwd filesep 'functions']))
pathToPsignifit = 'D:\NRG\tools\psignifit';
addpath(genpath(pathToPsignifit))

% directory where results are saved
resultsDir = [pwd '\data\signed-corr-mws-levels\subjects'];
subjects = {'AW', 'NG', 'LS'};

% experiment identifiers per subject
expId{1} = [4 1 2 6];
expId{2} = [4 1 2 6];
expId{3} = [1 2];

disparityMagnitude = {[1.5, 3, 9, 15], [1.5, 3, 9, 15], [3, 9]};

% colors for plotting
pltColor = [224, 141, 50;
            55, 126, 184;
            0.2, 0.2, 0.2;
            228, 26, 28;
            0.2, 0.2, 0.2;
            152,78,163;]/255;

f1 = figure('pos', [100,100,200,200],...
       'DefaultAxesFontName', 'Helvetica',...
       'DefaultAxesFontSize', 8);

subjectsToPlot = 2;
levelsToPlot = [4, 2, 6];
thisSubplot = 1;

for ss = 1:numel(subjects)
    
    halfWidth = [];

    for ee = expId{ss}
        
        runFiles = wildcardsearch(fullfile(resultsDir, subjects{ss}, 'results'),...
            ['exp' num2str(ee)]);
        
        if ~isempty(runFiles)
            
            if (any(subjectsToPlot==ss)) && (any(levelsToPlot==ee))
                hold on;
                if numel(subjectsToPlot)>1
                    subplot(1, numel(subjectsToPlot), thisSubplot);
                end
                plotFlag = 1;
            else
                plotFlag = 0;
            end
            results = plotResultsLevelsBootstrapped(subjects{ss}, ee, runFiles, pltColor(ee,:), 1, 'abs', plotFlag);
            halfWidth = [halfWidth results.halfWidth];
        end
        
    end
    if (any(subjectsToPlot==ss))
        thisSubplot = thisSubplot + 1;
    end
    thisDispMag = disparityMagnitude{ss};
    save(['hw-magnitude/' subjects{ss} '_halfWidth_magnitude.mat'], 'thisDispMag', 'halfWidth');

end

xlim([0 10])
plot2svg('figures/Fig4b.svg')

plotHalfWidthMagnitude;

