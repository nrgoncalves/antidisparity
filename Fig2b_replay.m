clc, clear all, close all; 
addpath(genpath([pwd filesep 'functions']))

% directory where results are saved
resultsDir = [pwd '\data\signed-corr-tms\subjects'];
subjects = {'AW', 'NG', 'LS', 'RR', 'WH', 'VK',...
            'EM', 'MG', 'ZR', 'PF', 'LL', 'SW'}; 

% experiment identifiers
expId = 3;

thr = zeros(numel(subjects), 2);

figure('pos', [200 700 150 200],...
       'DefaultAxesFontName', 'Helvetica',...
       'DefaultAxesFontSize', 8);
hold on;

% loop over subjects
for ss = 1:numel(subjects)
            
        % loop over experiments
        for ee = expId
            
            % get final run
            runFiles = wildcardsearch(fullfile(resultsDir, subjects{ss}, 'results'),...
                                               ['exp' num2str(ee) '*run_01.mat']);
            
            if numel(runFiles)==1
                runOutput = load(runFiles{1});
                
                % compute proportion correct for correlated and
                % anticorrelated RDS
                thr(ss,:) = [mean(runOutput.response(runOutput.design==1)==runOutput.disparities(runOutput.design==1)'),...
                             mean(runOutput.response(runOutput.design==2)==runOutput.disparities(runOutput.design==2)')];
                                
                plot([1,2], thr(ss,:), 'o-', 'Color',[0.6 0.6 0.6],...
                                       'markerfacecolor', 'w');
                
            else
                error('More than one results file found!')
            end
            
        end
            
end

xlim([0.5 2.5])
ylim([0 1])
ylabel('Proportion correct')
set(gca, 'xtick', [1, 2],'TickLength',[0.025 0.025],'TickDir', 'out')
set(gca, 'xticklabel', {'Corr', 'Anti'})

plot2svg('figures/Fig2b.svg')

% difference between correlated and anticorrelated performance
[hThr, pThr, CI, STATS] = ttest(thr(:,1), thr(:,2));

% difference between anticorrelated and chance
[hThrA, pThrA, CI_A, STATS_A] = ttest(thr(:,2), 0.5);



