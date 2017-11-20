function [trf, sourceTransformed] = coherentPointDrift(source, target)
% Register two point sets (source to target) using the Coherent Point Drift
% algorithm (Myronenko & Song, 2009, https://arxiv.org/abs/0905.2635).
% 
% !! IMPORTANT NOTE !!: This implementation is experimental and is still
% unstable.
%
% INPUT
%
% source            (arr[double])   coordinates for the points to be 
%                                   transformed (M-by-2)
%
% target            (arr[double])   coordinates for the target points (N-by-2)
%
% OUTPUT
% 
% trf               (struct)        transformation parameters
%  .T               (arr[double])   affine matrix (2-by-2)
%  .t               (arr[double])   x and y translation parameters (2-by-1)
% 
% sourceTransformed (arr[double])   transformed source points
%
% nrg, '02-Sep-2017 00:17:49'

% exclude nan points
source(isnan(source), :) = [];
target(isnan(target), :) = [];

% initial transform (identity)
T = eye(2);
t = [0, 0];

D = size(source, 2);
if size(target, 2) ~= D
    warning('Number of source and target features do not match.')
end

M = size(source, 1);
N = size(target, 1);

% initialize parameters
w = 0.001;  % [nrg] needs quite a bit of adjustment depending on how much noise we have
s2 = 1/(D*N*M) * sum(sum(pdist2(source, target)));  % variance

fh = figure;
hold on;
plot(target(:,1), target(:,2), 'ro')
plot(source(:,1), source(:,2), 'bx')
drawnow;

maxIter = 1000;
it = 1;
done = 0;
while (done == 0)
    
    % transformed data
    sourceTransformed = applyTransform(source, T, t);
    
    % E-step, compute P
    distanceMatrix = pdist2(sourceTransformed, target);
    
    pij = exp(-1 * distanceMatrix .^ 2 / (2*s2));
    pij = pij ./ (repmat(sum(pij, 1), [M, 1]) + (2*pi()*s2)^(D/2)*(w/(1-w))*(M/N));
    
    % M-step, solve
    [tmpT, tmpt, tmps2, success] = solve(target, source, pij);
    
    if success
        T = tmpT;
        t = tmpt';
        s2 = tmps2;
    else
        break;
    end
    
    figure(fh);
    cla;
    plot(target(:,1), target(:,2), 'ro')
    plot(sourceTransformed(:,1), sourceTransformed(:,2), 'bx')
    drawnow;
    % [nrg] could maybe work a bit more on the vizualization to plot the 
    % variance too
    
    % [nrg] this early stopping criteria does not work very well. I think
    % I should probably just run this inside a try block, catch the warning
    % and then return the previous iteration
    if it > maxIter
        done = 1;
    else
        if it ~= 1
           if max(max((pij_prev - pij).^2)) < eps
               done = 1;
               fprintf('Early stopping criteria met at iteration %i.\n', it)
           end
        end
        pij_prev = pij;      
        it = it + 1;
    end
        
end

sourceTransformed = applyTransform(source, T, t);
trf.T = T;
trf.t = t;


function [B, t, s2, success] = solve(S, M, P)
% S: matrix containing target points (N-by-D)
% T: matrix containing source points (the ones to adjust) (M-by-D)
% P: probabilities matrix

D = size(S, 2);
Np = ones(1, size(M, 1)) * P * ones(size(S, 1), 1);

mu_s = 1/Np * S' * P' * ones(size(S, 1), 1);
mu_m = 1/Np * M' * P' * ones(size(M, 1), 1);

Se = S - mu_s';
Me = M - mu_m';

mpm = Me' * diag(P*ones(size(P,2),1)) * Me;
rc = rcond(mpm);

if isnan(rc) || (rc < 10^(-5))
    success = 0;
    [B, t, s2] = deal([]);
else
    % B = (Se' * P' * Me) * inv(mpm);
    B = (Se' * P' * Me) / mpm;
    t = mu_s - B * mu_m;
    s2 = 1/(Np * D) * (trace(Se' * diag(P'*ones(size(P,1), 1)) * Se) - ...
                   trace(Se' * P' * Me * B'));
    success = 1;
end







