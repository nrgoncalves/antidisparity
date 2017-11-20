function results = plotResultsLevelsBootstrapped(subjId, expId, runFiles, pltColor, fitFlag, mode, plotFlag)
% plots the proportion of correct responses for several
% levels of difference between correlated and anticorrelated disparity

if nargin < 5
    fitFlag = 0;
end

if nargin < 6
    mode = 'abs';
end

if nargin < 7
    plotFlag = 0;
end

% load data from each run
design = [];
correctness = [];
for rr = 1:numel(runFiles)
    d = load(runFiles{rr});
    design = [design d.design];
    correctness = [correctness d.correctness];
end

% fetch some experiment parameters
sparam = d.sparam;
disparity = sparam.disparity;
dparam = d.dparam;
cDisp = d.cDisp;
aDisp = d.aDisp;

diffAxis = [disparity-(sparam.disparity(sparam.cDispInd(1))); disparity-sparam.disparity(sparam.cDispInd(2))];
if strcmp(mode, 'abs')
    diffAxis = abs(diffAxis);
end
diffAxisCond=reshape((1:26), [13,2])';

uniqueDiffs = unique(diffAxis);
pCorrect = zeros(numel(uniqueDiffs), 3);
nTrialsPerCondition = zeros(numel(uniqueDiffs),1);
nCorrectTrials = zeros(numel(uniqueDiffs),1);

for dd = 1:numel(uniqueDiffs)

    theseTrials = ismember(design, diffAxisCond(diffAxis==uniqueDiffs(dd)));
    nTrialsPerCondition(dd) = sum(theseTrials);
    nCorrectTrials(dd) = sum(correctness(theseTrials));
    bootstat = bootstrp(5000,@mean,correctness(theseTrials));
    ci = prctile(bootstat, [50 84 16]);
    pCorrect(dd,:) = [ci(1), abs(ci(2:3)-ci(1))];
    
end

if strcmp(mode, 'none')
    zeroPoint = find(uniqueDiffs==0);
    left = pCorrect(1:(zeroPoint-1), 1);
    right = pCorrect(end:-1:(zeroPoint+1), 1);

    [r, p] = corr(left, right);

    save(['lr-correlation/' subjId '_exp' num2str(expId) 'LevelsTuningCurveLeftRightCorrelation.mat'], 'left', 'right', 'r', 'p');
end

maxHeight = max(pCorrect(:, 1));
minHeight = min(pCorrect(:, 1));
halfHeight = (maxHeight + minHeight)/2;

diffs = uniqueDiffs;

hold on;

if fitFlag
    data = [diffs, nCorrectTrials, nTrialsPerCondition];
    
    options = struct;   % initialize as an empty struct
    options.sigmoidName = 'norm';
    options.priors = cell(5,1);
    options.borders = nan(5,2);
    options.expType = 'YesNo';
    
    options.priors{4} = @(x) betapdf(x,2,2);
    options.borders(3,:) = [0,.1];
    options.borders(4,:) = [.11,.89];
    options.fixedPars = nan(5,1);
    options.fixedPars(5) = 0;
    options.stepN   = [40,40,40,40,1];
    options.mbStepN = [30,30,20,20,1];
    disp('Fitting psychometric function...');

    res = psignifit(data, options);
    
    plotOptions.dataColor = pltColor;
    if plotFlag
        plotPsych(res, plotOptions);
    end
    [halfWidth, CI] = getThreshold(res, halfHeight, 0);
	
else
    if plotFlag
        h1 = shadedErrorBar(diffs, pCorrect(:,1), pCorrect(:,2:3), {'ko', 'markerfacecolor', pltColor,'markeredgecolor','none'}, 1);
        set(h1.patch, 'faceColor', pltColor);
        set(h1.edge(1), 'color', pltColor);
        set(h1.edge(2), 'color', pltColor);
    end
end

if plotFlag
    plot([diffs(1)-1, diffs(end)+1], 0.5*ones(1,2), 'k--');
    ylim([0,1.1])
    xlim([diffs(1)-1, diffs(end)+1])
    set(gca, 'xtick', diffs(1:2:end))
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.025 0.025])
    ylabel('Proportion correct')
    xlabel('Disparity difference (arcmin)')
    title(subjId)
end
if fitFlag
    results.res = res;
    results.halfWidth = halfWidth;
else
    results = [];
end

