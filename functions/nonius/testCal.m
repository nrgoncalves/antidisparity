function res = testCal(testFile, design, params)
% Test the calibration parameters on data with known reference points.
% 
% INPUT
%
% testFile          (string)        path to mat file containing raw and
%                                   event data
% 
% design            (struct)        structure containing the position and
%  .xcoords                         onset/offset times for reference data
%  .ycoords                         (see Documentation for more info)
%  .onsetTimes
%  .offsetTimes
%
% params            (struct)        structure containing calibration
%                                   parameters for the left and right eyes
%
% OUTPUT
%
% res               (struct)        structure containing test results
%  .(eye).V.r2                      variance for vertical (V.r2) and 
%  .(eye).H.r2                      horizontal (H.r2) coordinates, as well as
%  .(eye).r2                        the average explained variance (.r2)
% 
% nrg, '01-Sep-2017 19:13:30'


dstruct = load(testFile);
eyes = fieldnames(dstruct.E);

for e = 1:numel(eyes)
    
    % test data
    gaze = dstruct.E.(eyes{e});
    vY_test = gaze.V;
    hY_test = gaze.H;

    % design
    time = (gaze.T - gaze.T(1))/10^3;
    [hX, vX] = deal(zeros(size(time)));
    
    theseOnsets = design.onsetTimes + params.(eyes{e}).movdelay;
    theseOffsets = design.offsetTimes + params.(eyes{e}).movdelay;
    
    for i = 1:numel(design.xcoords)
        hX((time>theseOnsets(i)) & (time<theseOffsets(i))) = design.xcoords(i);
        vX((time>theseOnsets(i)) & (time<theseOffsets(i))) = design.ycoords(i);
    end
    
    flipHorizontal = 1;
    if flipHorizontal == 1
        hX = -1 * hX;
    end
    
    % apply calibration parameters
    vye = (vY_test - params.(eyes{e}).vb(2)) / params.(eyes{e}).vb(1);
    hye = (hY_test - params.(eyes{e}).hb(2)) / params.(eyes{e}).hb(1);
    
    % apply transform
    applyAffineCorrection = 1;
    if applyAffineCorrection == 1
        xyTrf = applyTransform([hye, vye], params.(eyes{e}).trf.T, params.(eyes{e}).trf.t, 'inverted');
        hye = xyTrf(:, 1);
        vye = xyTrf(:, 2);
    end
    
    % plotting
    figure('pos', [100 200 800 200]);
    subplot(1,3,1)
    plot(time, vye)
    hold on; plot(time, vX)
    ylabel('Gaze position')
    xlabel('Time (s)')
    
    subplot(1,3,2)
    plot(time, hye)
    hold on; plot(time, hX)
    xlabel('Time (s)')
    
    subplot(1,3,3)
    plot(hye, vye, 'k.')
    hold on; plot(-1*design.xcoords, design.ycoords, 'ro')
    ylabel('Vertical position (deg)')
    xlabel('Horizontal position (deg)')
    xlim([-7, 7])
    ylim([-7, 7])
    
    res.(eyes{e}).V.ye = vye;
    res.(eyes{e}).V.xe = vX;
    res.(eyes{e}).H.ye = hye;
    res.(eyes{e}).H.xe = hX;
    
    % compute explained variance (with respect to ground truth design)
    exclude = isnan(vye) | isnan(vX);
    vye(exclude) = [];
    vX(exclude) = [];
    vr2 = 1 - (std(vye-vX)^2) / (std(vye)^2);
    
    exclude = isnan(hye) | isnan(hX);
    hye(exclude) = [];
    hX(exclude) = [];
    hr2 = 1 - (std(hye-hX)^2) / (std(hye)^2);
    
    res.(eyes{e}).V.r2 = vr2;
    res.(eyes{e}).H.r2 = hr2;
    res.(eyes{e}).r2 = mean([vr2, hr2]);
    
    fprintf('Proportion of explained variance (R^2) = %0.2f\n', res.(eyes{e}).r2)
    
end