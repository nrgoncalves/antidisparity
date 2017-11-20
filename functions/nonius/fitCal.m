function res = fitCal(calFile, design, distCorr)
% Fit the calibration parameters on data with known reference points.
%
% INPUT
%
% calFile           (string)        path to mat file containing raw and
%                                   event data
%
% design            (struct)        structure containing the position and
%  .xcoords                         onset/offset times for reference data
%  .ycoords                         (see Documentation for more info)
%  .onsetTimes
%  .offsetTimes
%
% distCorr          (string)        method for distortion correction
%                                   { 'none' | 'manual' | 'auto_cpd' }
%
%
% OUTPUT
%
% res               (struct)        structure containing fit results
%
%  .(eye).r2        (double)        proportion of variance explained by fit
%  .(eye).vb        (arr[double])   array of parameters ([slope, intersect])
%                                   of the fitted linear model to vertical
%                                   coordinates
%  .(eye).hb        (arr[double])   array of parameters ([slope, intersect])
%                                   of the fitted linear model to
%                                   horizontal coordinates
%  .(eye).movdelay  (double)        estimated delay between stim onset and
%                                   eye movement
%
%
% nrg, '01-Sep-2017 23:48:30'


if nargin < 3
    distCorr = 'none';
end

dstruct = load(calFile);
eyes = fieldnames(dstruct.E);

for e = 1:numel(eyes)

    % calibration data
    gaze = dstruct.E.(eyes{e});

    % exclude blink periods +- 25 ms
    blinks = zeros(numel(gaze.T), 1);
    for i = 1:dstruct.EVT.blink.n
        blinks((dstruct.E.L.T > dstruct.EVT.blink.Tstart(i)-25) & (dstruct.E.L.T < dstruct.EVT.blink.Tend(i)+25)) = 1;
    end

    gaze.H(blinks==1) = nan;
    gaze.V(blinks==1) = nan;

    movementDelay = linspace(0.100, 0.400, 16);
    fprintf('Estimating %s eye position\n', eyes{e})
    r2 = 0;
    for d = 1:numel(movementDelay)

        fprintf('... [delay = %.3f s', movementDelay(d))

        time = (gaze.T - gaze.T(1))/10^3;
        [hX, vX] = deal(zeros(size(time)));

        % add the delay
        theseOnsets = design.onsetTimes + movementDelay(d);
        theseOffsets = design.offsetTimes + movementDelay(d);

        % build regressor
        for i = 1:numel(design.xcoords)
            hX((time>theseOnsets(i)) & (time<theseOffsets(i))) = design.xcoords(i);
            vX((time>theseOnsets(i)) & (time<theseOffsets(i))) = design.ycoords(i);
        end

        flipHorizontal = 1;
        if flipHorizontal == 1
            hX = -1 * hX;
        end

        % exclude nan samples
        exclude = isnan(gaze.V) | isnan(gaze.H);

        vY = gaze.V;
        hY = gaze.H;

        vY(exclude) = [];
        vX(exclude) = [];

        hY(exclude) = [];
        hX(exclude) = [];

        time(exclude) = [];

        % add bias term
        vX(:, 2) = ones(size(vX));
        hX(:, 2) = ones(size(hX));

        % [nrg] later we might include a first order term

        % least-squares regression
        vB = regress(vY, vX);
        vYe = vX*vB;
        vR2 = 1 - (std(vY-vYe)^2)/(std(vY)^2);

        hB = regress(hY, hX);
        hYe = hX * hB;
        hR2 = 1 - (std(hY-hYe)^2)/(std(hY)^2);

        R2 = (vR2 + hR2)/2;
        fprintf(', R = %.3f]\n', R2)

        % check if this is the best model so far
        if R2 > r2
            r2 = R2;        % average r2 statistic (horizontal and vertical)
            vb = vB;        % parameters for vertical position
            hb = hB;        % parameters for horizontal poisition
            v = vY;         % measured vertical position (pixel)
            h = hY;         % measured horizontal position (pixel)
            ve = vYe;       % estimated vertical position (pixel)
            he = hYe;       % estimated horizontal position (pixel)
            t = time;       % time
            movDelay = movementDelay(d);    % movement delay (seconds)
        end
    end

    res.(eyes{e}).r2 = r2;
    res.(eyes{e}).vb = vb;
    res.(eyes{e}).hb = hb;
    res.(eyes{e}).movdelay = movDelay;

    hye = (hY-hb(2))/hb(1);
    vye = (vY-vb(2))/vb(1);

    % distortion correction
    source = [-1*design.xcoords, design.ycoords];
    target = zeros(numel(design.xcoords), 2);
    for i = 1:numel(design.xcoords)
        target(i,:) = [nanmedian(hye(hX==design.xcoords(i) & vX==design.ycoords(i))), ...
                       nanmedian(vye(hX==design.xcoords(i) & vX==design.ycoords(i)))];
    end

    switch distCorr
        case 'none'
            trf.T = eye(2);
            trf.t = zeros(1, 2);
        case 'manual'
            trf = manualRegistration(source, target);
        case 'auto_cpd'
            trf = psreg(source, target);
    end

    res.(eyes{e}).trf = trf;

    % plotting
    figure('pos', [100 700 800 200]);
    subplot(1,3,1)
    plot(t, vX(:,1))
    hold on; plot(t, vye)
    ylabel('Gaze position')
    xlabel('Time (s)')
    title('Vertical')
    xlim([t(1), t(end)]);

    subplot(1,3,2)
    plot(t, hX(:, 1))
    hold on; plot(t, hye)
    xlabel('Time (s)')
    title('Horizontal')
    xlim([t(1), t(end)]);

    subplot(1,3,3)
    plot(hye, vye, 'k.')
    hold on; plot(-1*design.xcoords, design.ycoords, 'ro')
    ylabel('Vertical position (deg)')
    xlabel('Horizontal position (deg)')
    xlim([-7,7])
    ylim([-7,7])
    
    
    
end
