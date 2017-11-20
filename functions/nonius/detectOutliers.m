function outliers = detectOutliers(x)
% Detect outliers using k-NN with 2nd order derivative and range features.
%
% INPUT
%
% x         (arr[double])       epoched data (M epochs-by-T timepoints)
%
%
% OUTPUT
%
% outliers  (arr[logical])      M-by-1 vector containing 1's if epoch is
%                               outlier
%
% nrg, '02-Sep-2017 00:55:38'

dx = diff(x, 1, 2);
x1 = zscore(mean(dx.^2, 2));
x2 = zscore(max(dx,[],2)- min(dx,[],2));

X = [x1, x2];
ds = sort(squareform(pdist(X, 'euclidean')), 2, 'ascend');

k = 8;
knnd = mean(ds(:, 2:2+k-1), 2);

dcutoff = 1;
outliers = knnd > dcutoff;

%{
figure;
plot(this.epochTime.Stim1, x(outliers == 0, :), 'color', 'k')
hold on;
plot(this.epochTime.Stim1, x(outliers == 1, :), 'color', 'r')
pause;
%}
