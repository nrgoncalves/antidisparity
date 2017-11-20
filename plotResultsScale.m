function plotResultsScale(subjId, expId, runFiles, fitFlag)
% plots the proportion of correct responses for several
% levels of visuotopic separation between correlated and anticorrelated dots

colorSame = [127,201,127]/255;
colorOpp = [190,174,212]/255;

% load data from each run
design = [];
correctness = [];
reactionTime = [];
for rr = 1:numel(runFiles)
    d = load(runFiles{rr});
    design = [design d.design];
    correctness = [correctness d.correctness];
    reactionTime = [reactionTime d.reactionTime];
end

spatialFactor = 2.^(1:7);
nSpatialFactors = numel(spatialFactor);
spatialFactorCond = repmat(1:nSpatialFactors, [1,2]);
dparam.numConds = numel(spatialFactorCond);
invert_anti = [zeros(1,nSpatialFactors), ones(1,nSpatialFactors)];

% plot results
pCorrectSame = zeros(3,nSpatialFactors);
pCorrectOpp = zeros(3,nSpatialFactors);

nCorrectTrials = zeros(nSpatialFactors,1);
nTrialsPerCondition = zeros(nSpatialFactors,1);

rTimeSame = zeros(1,nSpatialFactors);
rTimeOpp = zeros(1,nSpatialFactors);

iSame = 1;
iOpp = 1;
for jj = 1:numel(spatialFactorCond)
    
    pCorrect = mean(correctness(design==jj));
    rTime = mean(reactionTime(design==jj));
    
    bootstat = bootstrp(5000,@mean,correctness(design==jj));
    pCorrect = prctile(bootstat, [50 84 16]); % mean, upper and lower limits
    
    pCorrect(2:3) = abs(pCorrect(2:3) - pCorrect(1));
    
    if invert_anti(jj) == 0
        nTrialsPerCondition(iSame) = sum(design==jj);
        nCorrectTrials(iSame) = sum(correctness(design==jj));
        
        pCorrectSame(:, iSame) = pCorrect';
        rTimeSame(iSame) = rTime;
        iSame = iSame + 1;
    elseif invert_anti(jj) == 1
        pCorrectOpp(:, iOpp) = pCorrect';
        rTimeOpp(iOpp) = rTime;
        iOpp = iOpp + 1;
    else
        sca;
        keyboard;
    end
end

spatialPeriod = 60*spatialFactor/d.sparam.pix_per_deg;

plot(spatialPeriod, 0.5*ones(size(spatialFactor)), 'k--');
hold on;
h1 = shadedErrorBar(spatialPeriod, pCorrectSame(1,:), pCorrectSame(2:3,:), {'ko', 'markerfacecolor', colorSame,'markeredgecolor','none'}, 1);
set(h1.patch, 'faceColor', colorSame);
set(h1.edge(1), 'color', colorSame)
set(h1.edge(2), 'color', colorSame)
h2 = shadedErrorBar(spatialPeriod, pCorrectOpp(1,:), pCorrectOpp(2:3,:), {'ko', 'markerfacecolor', colorOpp,'markeredgecolor','none'}, 1);
set(h2.patch, 'faceColor', colorOpp);
set(h2.edge(1), 'color', colorOpp)
set(h2.edge(2), 'color', colorOpp)

if fitFlag
    data = [spatialPeriod', nCorrectTrials, nTrialsPerCondition];

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

    plotOptions.dataColor = colorSame;
    plotPsych(res, plotOptions);    
end


set(gca,'TickDir','out')
set(gca,'TickLength',[0.025 0.025])

ylim([0,1.05])
xlim([spatialPeriod(1)-1, spatialPeriod(end)+1])
ylabel('Proportion correct')
xlabel('Spatial period (arcmin)')
title(subjId)
