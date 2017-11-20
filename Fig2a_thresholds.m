clc, clear all, close all; 
addpath(genpath([pwd filesep 'functions']))

% directory where results are saved
resultsDir = [pwd '\data\signed-corr-tms\subjects'];
subjects = {'AW', 'NG', 'LS', 'RR', 'WH', 'VK',...
            'EM', 'MG', 'ZR', 'PF', 'LL', 'SW'}; 

% experiment identifier
expId = 2;

thr = zeros(numel(subjects), 2);

figure('pos', [200 700 150 200],...
       'DefaultAxesFontName', 'Helvetica',...
       'DefaultAxesFontSize', 8);

% loop over subjects
for ss = 1:numel(subjects)
        
        % loop over experiments
        for ee = expId
            
            % find last run of QUEST
            runFiles = wildcardsearch(fullfile(resultsDir, subjects{ss}, 'results'),...
                                               ['exp' num2str(ee) '*run_02.mat']);
            
            % plot threshold value
            if numel(runFiles)==1
                runOutput = load(runFiles{1});
                thr(ss,:) = 10 .^ runOutput.thr_mean;
                            
                hold on;
                plot([1,2], thr(ss,:), 'o-', 'Color',[0.6 0.6 0.6],...
                                       'markerfacecolor', 'w');                
            else
                error('More than one results file found!')
            end
            
        end
            
end

xlim([0.5 2.5])
ylim([0 1])
ylabel('Proportion of correlated dots')
set(gca, 'xtick', [1, 2],'TickLength',[0.025 0.025],'TickDir', 'out')
set(gca, 'xticklabel', {'Same', 'Opposite'})

% two-tailed, paired t-test on thresholds for the two conditions
[hThr,pThr,CI,STATS] = ttest(thr(:,1), thr(:,2));

% save plot
plot2svg('figures/Fig2a.svg')


