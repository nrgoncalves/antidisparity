figure;
subplot(1,1,1);

% first we plot the data collected by Stevenson and colleagues (1992)
load('hw-magnitude/stevData_converted.mat')

xx = scx;
yy = scy/2; % convert from full-width to half-width

% bootstrap linear fit
nSamples = 1000;
m = numel(xx);

xi_max = 16.5;
samples = randi(m, [nSamples, m]);
xi = linspace(0,xi_max,50);
yi = zeros(nSamples, size(xi, 2));
for i = 1:nSamples

    thisx = xx(samples(i,:));
    thisy = yy(samples(i,:));

    p = polyfit(thisx, thisy, 1);

    yi(i, :) = p(1)*xi + p(2);
    hold on;
end

ciBounds = [2.5 97.5];
ci = zeros(size(xi, 2), 2);
for j = 1:size(xi, 2)
    prc = prctile(yi(:, j), ciBounds);
    ci(j, :) = prc;
end

plot(xi, ci(:, 1), 'r--')
plot(xi, ci(:, 2), 'r--')
plot(xi, mean(yi, 1), 'r')

% then we plot fMRI estimates from V3A and V3B 
% (Goncalves et al., 2015, J Neurosci)
load('hw-magnitude/V3A_fwhm_fits.mat')
p = p/2; % convert from fwhm to hwhm
plot(xi, p(1)*xi + p(2), 'g')
disp(p(1))

load('hw-magnitude/V3B_fwhm_fits.mat')
p = p/2; % convert from fwhm to hwhm
plot(xi, p(1)*xi + p(2), 'c')
disp(p(1))

% finally we plot our latest behavioural data
subjects = {'AW', 'NG', 'LS'};
lspec = {'^', 'o', 's'};
cspec = [145, 118, 182; 71, 169, 71; 252 147 48] / 255;
hold on;

disparityMagnitude = [];
halfWidth = [];
for ss = 1:numel(subjects)
    
    res = load(['hw-magnitude/' subjects{ss} '_halfWidth_magnitude.mat']);
    
    disparityMagnitude = [disparityMagnitude res.thisDispMag];
    halfWidth = [halfWidth res.halfWidth];
    
    plot(res.thisDispMag, res.halfWidth, 'LineStyle', 'none', 'Marker', lspec{ss}, 'color', cspec(ss,:), 'linewidth', 2)
    
end

[r, p] = corr(disparityMagnitude', halfWidth');

set(gca,'TickDir','out')
set(gca,'TickLength',[0.025 0.025])
ylabel('Half width (arcmin)')
xlabel('Disparity magnitude (arcmin)')
xlim([0, 16.5]);
ylim([0, 8]);
axis square;

plot2svg('figures/Fig4c.svg')
