function plotResultsTime(subjId, expId, runFiles)
% plots the proportion of correct responses for several
% levels of difference between correlated and anticorrelated dots

colorPre = [190,174,212]/255;
colorPost = [127,201,127]/255;

% load data from each run
design = [];
correctness = [];
for rr = 1:numel(runFiles)
    d = load(runFiles{rr});
    design = [design d.design];
    correctness = [correctness d.correctness];
end

pcorrect = zeros(d.dparam.numConds, 3);
for c = 1:d.dparam.numConds
    bootstat = bootstrp(5000,@mean,correctness(design==c));
    pcorrect(c,:) = prctile(bootstat, [50 84 16]);
    pcorrect(c, 2:3) = abs(pcorrect(c,2:3) - pcorrect(c,1));
end

caDelay = [-4:-1, 1:4];
frameRate = 120;
tAxis = caDelay/frameRate*1000;

h1 = shadedErrorBar(tAxis(1:4), pcorrect(1:4,1), pcorrect(1:4, 2:3), {'ko', 'markerfacecolor', colorPre,'markeredgecolor','none'}, 1);
set(h1.patch, 'faceColor', colorPre);
set(h1.edge(1), 'color', colorPre)
set(h1.edge(2), 'color', colorPre)

h2 = shadedErrorBar(tAxis(5:8), pcorrect(5:8,1), pcorrect(5:8, 2:3), {'ko', 'markerfacecolor', colorPost,'markeredgecolor','none'}, 1);
set(h2.patch, 'faceColor', colorPost);
set(h2.edge(1), 'color', colorPost)
set(h2.edge(2), 'color', colorPost)

set(gca,'TickDir','out')
set(gca,'TickLength',[0.025 0.025])

% plot(caDelay(1:4)/frameRate*1000, pcorrect(1:4), 'o-');
% hold on;
% plot(caDelay(5:end)/frameRate*1000, pcorrect(5:end), 'o-');
xlabel('Correlated - anticorrelated dots onset (ms)')
ylabel('Proportion correct')
ylim([0.4,1.05]);
title(subjId);


% pooled by positive or negative onset
bootstat = bootstrp(5000,@mean,correctness(design<=4));
pcorrectPre = prctile(bootstat, [50 84 16]);

bootstat = bootstrp(5000,@mean,correctness(design>5));
pcorrectPost = prctile(bootstat, [50 84 16]);

fig2 = figure;
hold on;
bar(1, pcorrectPre(1),'facecolor', colorPre); 
plot([1, 1], pcorrectPre(2:3), 'color', 'k')

bar(2, pcorrectPost(1),'facecolor', colorPost); 
plot([2, 2], pcorrectPost(2:3), 'color', 'k')
set(gca,'TickDir','out')
set(gca,'TickLength',[0.025 0.025])

ylabel('Proportion correct')
ylim([0.4,1.05]);
title(subjId);
set(gca, 'xtick', [1, 2])
set(gca, 'xticklabel', {'Anti lead', 'Corr lead'})

if expId == 1
    plot2svg(['figures/Fig6c_' subjId '.svg'])
elseif expId == 4
    plot2svg(['figures/Fig6_' subjId '_longerPresentationTime.svg'])
else
    plot2svg(['figures/Fig6d_control_' subjId '.svg'])
end
close(fig2);

