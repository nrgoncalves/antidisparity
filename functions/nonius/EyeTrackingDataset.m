 classdef EyeTrackingDataset < handle
    % Handle eye-tracking datasets.
    
    % base properties
    properties
        files
        meas
        evnt
        calibration
        validation
        epochs
        epochTime
        conditions
        samplingRate
    end
    properties (Access = private)
        state
    end
    
    
    % methods that handle the data
    methods
        
        % Constructor
        function obj = EyeTrackingDataset(files, calfiles, valfiles)
            
            if nargin < 2
                calfiles = {};
            end
            if nargin < 3
                valfiles = {};
            end
            if ~iscell(files)
                files = {files};
            end
            if ~iscell(calfiles)
                calfiles = {calfiles};
            end
            if ~iscell(valfiles)
                valfiles = {valfiles};
            end
            
            obj.files = files;
            for i = 1:numel(files)
                runData = load(files{i});
                obj.evnt{i} = runData.EVT;
                obj.meas{i} = runData.E;
                obj.samplingRate{i} = 10^3 / (runData.E.L.T(2) - runData.E.L.T(1));
            end
            
            % calibration and validation data
            obj.calibration.files = calfiles;
            obj.validation.files = valfiles;
            obj.state = 'raw';  % this really should be private
        end
        
        
        % Set calibration files
        function setCalibrationFiles(this, calfiles)
            if ~iscell(calfiles)
                calfiles = {calfiles};
            end
            this.calibration.files = calfiles;
        end
        
        
        % Set validation files
        function setValidationFiles(this, valfiles)
            if ~iscell(valfiles)
                valfiles = {valfiles};
            end
            this.validation.files = valfiles;
        end
        
        
        % Set calibration design
        function setCalibrationDesign(this, design)
            % [nrg] TODO: check if design structure is valid
            this.calibration.design = design;
        end
        
        
        % Set validation design
        function setValidationDesign(this, design)
            % [nrg] TODO: check if design structure is valid
            this.validation.design = design;
        end
        
        
        % Estimate calibration parameters
        function estimateCalibrationParams(this, distCorr)
            
            if nargin < 2
                distCorr = 'none';
            end
                        
            if isempty(this.calibration.files)
                error(['No calibration files have been defined. ',...
                       'You can specify calibration files in EyeTrackingDataset() ',...
                       'or using the method setCalibrationFiles(...)'])
            else
                calFiles = this.calibration.files;
            end
            
            nCalFiles = numel(calFiles);
            
            if ~any(strcmp(fieldnames(this.calibration), 'design'))
                error(['No calibration design has been specified. ',...
                       'You can specify the calibration design using ',...
                       'the method setCalibrationDesign(...)'])
            end
            
            switch numel(this.calibration.design)
                case 1
                    fprintf('One calibration design found, using it for all calibration files...');
                    design(1:nCalFiles) = this.calibration.design;
                case numel(calFiles)
                    design = this.calibration.design;
                otherwise
                    error('The number of specified calibration designs does not match the number or calibration files')
            end
            
            
            this.calibration.fit(1:nCalFiles) = struct('L', [], 'R', []);
            
            for i = 1:nCalFiles
                
                % estimate
                res = fitCal(calFiles{i}, design(i), distCorr);
                                
                % store params
                this.calibration.fit(i) = res;
                
            end
            
        end
        
        
        % Validate calibration parameters
        function validateCalibrationParams(this) 
            if isempty(this.validation.files)
                error(['No validation files have been defined. ',...
                       'You can specify calibration files in EyeTrackingDataset() ',...
                       'or using the method setCalibrationFiles(...)'])
            else
                valFiles = this.validation.files;
                nValFiles = numel(valFiles);
                calFiles = this.calibration.files;
                nCalFiles = numel(calFiles);
            end
            
            if nValFiles ~= nCalFiles
                error(['The number of calibration files (%i) does not match ', ...
                       'the number of validation files (%i)'], ...
                       nCalFiles, nValFiles);
            end
            
            switch numel(this.validation.design)
                case 1
                    fprintf('One validation design found, using it for all validation files...\n');
                    design(1:nValFiles) = this.validation.design;
                case nValFiles
                    design = this.validation.design;
                otherwise
                    error('The number of specified calibration designs does not match the number or calibration files')
            end
            
            this.validation.res(1:nValFiles) = struct('L', [], 'R', []);
            
            for i = 1 : nValFiles
                
                % test using the calibration parameters
                res = testCal(valFiles{i}, design(i), this.calibration.fit(i));
                
                % store results in validation structure
                this.validation.res(i) = res;
                
            end
        end
        
        
        % Apply parameters to measurement data
        function applyCalibrationParams(this)
            
            if strcmp(this.state, 'calibrated')
                error(['The dataset appears to have been calibrated previously. ', ...
                       'Perhaps you called applyCalibrationParams twice on this object?'])
            end
            
            nFiles = numel(this.meas);
            for i = 1 : nFiles
                params = this.calibration.fit(i);
                eyes = fieldnames(params);
                for e = 1 : numel(eyes)
                    this.meas{i}.(eyes{e}).H = (this.meas{i}.(eyes{e}).H - params.(eyes{e}).hb(2)) / params.(eyes{e}).hb(1);
                    this.meas{i}.(eyes{e}).V = (this.meas{i}.(eyes{e}).V - params.(eyes{e}).vb(2)) / params.(eyes{e}).vb(1);
                end
            end
            this.state = 'calibrated';
        end
                
        
        % Epoch data around marker
        function epochAround(this, markerName, timeWindow, epochName)
                        
            if nargin < 4
                epochName = markerName;
            end
            
            ep = cell(numel(this.meas), 1);
            
            for r = 1:numel(this.meas)
            
                timeWindow = timeWindow / 10^3 * this.samplingRate{r};
                t = this.evnt{r}.msg.time(cellfun(@(x) ~isempty(strfind(x, markerName)), this.evnt{r}.msg.text));
                ep{r} = cell(numel(t), 1);
                
                for i = 1:numel(t)
                    idx = find(this.meas{r}.L.T == t(i));
                    frames = (idx + timeWindow(1)):(idx + timeWindow(2));
                    ep{r}{i} = applyToFields(this.meas{r}, @(x) x(frames));
                
                end
                
            end
            
            this.epochs.(epochName) = ep;
            this.epochTime.(epochName) = timeWindow(1):timeWindow(2);
        end
        
        
        % Group all trials per condition
        function groupByCondition(this)
            
            conds = fieldnames(this.epochs);
            
            for i = 1:numel(conds)
                eyes = fieldnames(this.epochs.(conds{i}){1}{1});
                for e = 1:numel(eyes)
                    vars = fieldnames(this.epochs.(conds{i}){1}{1}.(eyes{e}));
                    for v = 1:numel(vars)
                        this.conditions.(conds{i}).(eyes{e}).(vars{v}) = this.stackTrials(conds{i}, eyes{e}, vars{v});
                    end
                end
            end
            
        end
        
        
        % Reject outliers
        function rejectOutliers(this, conds)
            if nargin < 2
                conds = fieldnames(this.conditions);
            end
            if isempty(conds)
                conds = fieldnames(this.conditions);
            end
            if ~iscell(conds)
                conds = {conds};
            end
            
            for c = 1:numel(conds)
                % remove blinks
                T = this.conditions.(conds{c}).L.T;

                blinks = zeros(size(T));
                for r = 1:numel(this.evnt)  % loop through runs
                    for b = 1:this.evnt{r}.blink.n  % loop through blinks
                        blinks(T>=this.evnt{r}.blink.Tstart(b)-100 & T<=this.evnt{r}.blink.Tend(b)+100) = 1;
                    end
                end
                
                blinkTrials = any(blinks, 2);
                this.conditions.(conds{c}).L = applyToFields(this.conditions.(conds{c}).L, @(x) x(~blinkTrials, :));
                this.conditions.(conds{c}).R = applyToFields(this.conditions.(conds{c}).R, @(x) x(~blinkTrials, :));
                fprintf('[BLINK REJECTION] Removing %i trials (%.2f %%) for condition %s\n', ...
                    sum(blinkTrials), mean(blinkTrials)*100, conds{c})
                
                % run outlier detection
                outliersL = detectOutliers(this.conditions.(conds{c}).L.H);
                outliersR = detectOutliers(this.conditions.(conds{c}).R.H);
                reject = outliersL | outliersR;
                fprintf('[OUTLIER REJECTION] Rejecting %i trials (%.2f %%) for condition %s\n', ...
                        sum(reject), mean(reject)*100, conds{c})
                this.conditions.(conds{c}).L = applyToFields(this.conditions.(conds{c}).L, @(x) x(~reject, :));
                this.conditions.(conds{c}).R = applyToFields(this.conditions.(conds{c}).R, @(x) x(~reject, :));
            
                this.conditions.propRejected = mean(reject);
            
            end
        end
        
        
        % Stack trials into 2d array
        function S = stackTrials(this, epochName, whichEye, whichVar)
            
            data = this.epochs.(epochName);
            timeAx = this.epochTime.(epochName);
            
            nRuns = numel(data);
            nTrials = 0;
            for r = 1 : nRuns
                nTrials = nTrials + numel(data{r});
            end

            S = zeros(nTrials, numel(timeAx));
            trialNum = 1;
            for r = 1 : nRuns
                nRunTrials = numel(data{r});
                for t = 1 : nRunTrials
                    S(trialNum, :) = data{r}{t}.(whichEye).(whichVar);
                    trialNum = trialNum + 1;
                end
            end
            
        end
                
        
        % PLOTTING
        % Plot mean and sd for one or more conditions
        function plotConditions(this, conds, whichEye, field)
            
            figure;
            pcolors = get(gca, 'ColorOrder');
            
            cHandles = zeros(1, numel(conds));
            mver = version();
            
            for c = 1:numel(conds)

                if mver(1)=='9'
                    h = shadedErrorBar(this.epochTime.(conds{c}), ...
                                   this.conditions.(conds{c}).(whichEye).(field), ...
                                   {@(x) mean(x,1), @(x) std(x,[],1)/sqrt(size(x,1))}, ...
                                   'lineProps', {'color', pcolors(c,:)});
                    hold on;
                else
                    hold on;
                    warning('The shaded error bar visualization is only available for recent Matlab 2017. Default to simple line plots')
                    my = mean(this.conditions.(conds{c}).(whichEye).(field), 1);
                    sdy = std(this.conditions.(conds{c}).(whichEye).(field), [], 1);
                    h = plot(this.epochTime.(conds{c}), ...
                             my,...
                             'color', pcolors(c,:));
                    hl = plot(this.epochTime.(conds{c}), ...
                             my - sdy, '--',...
                             'color', pcolors(c,:));
                    hu = plot(this.epochTime.(conds{c}), ...
                             my + sdy, '--',...
                             'color', pcolors(c,:));
                    
                    cHandles(c) = h;
                end

                cHandles(c) = h.mainLine;
            end
            xlabel('Time (ms)')
            % add ylabel
            ylabel(['Eye: ' whichEye ', Field :' field])
            
            % add legend for conditions
            legend(cHandles, conds, 'Box', 'off')
        end
        
        
        % Plot raw epoched data (useful for quality control)
        function plotEpoch(this, epochName, whichEye, whichVar)
            
            if ~iscell(whichEye)
                whichEye = {whichEye};
            end
            
            if ~iscell(whichVar)
                whichVar = {whichVar};
            end
            
            data = this.epochs.(epochName);
            timeAx = this.epochTime.(epochName);
            
            nRuns = numel(data);
            
            nEyes = numel(whichEye);
            for e = 1 : nEyes
                hf = figure;
                subplot(1, nEyes, e);
                hold on;
                
                for v = 1 : numel(whichVar)
                    for r = 1 : nRuns
                        nTrials = numel(data{r});
                        for t = 1 : nTrials
                            plot(timeAx, data{r}{t}.(whichEye{e}).(whichVar{v}))
                        end
                    end
                end
                xlabel('Time (ms)')
            end
            
        end
        
        
        % Plot heatmap of eye position
        function map = plotHeatmap(this, cond, whichEye, histBins)
            
            if nargin < 4
                histBins = [];
            end
            
            xy = [this.conditions.(cond).(whichEye).H(:), ...
                  this.conditions.(cond).(whichEye).V(:)];
            
            if isempty(histBins)
                bounds = [min(xy, [], 1); max(xy, [], 1)];
                xbins = linspace(bounds(1,1), bounds(2,1), 50);
                ybins = linspace(bounds(1,2), bounds(2,2), 50);
                histBins = {xbins, ybins};
            else
                xbins = histBins{1};
                ybins = histBins{2};
            end

            xbinCenters = xbins(1:end-1) + diff(xbins)/2;
            ybinCenters = ybins(1:end-1) + diff(ybins)/2;
            
            % NOTE: we are transposing the result so that X corresponds to
            % horizontal coordinates
            rawMap = hist3(xy, 'Edges', histBins)';
            map = imfilter(rawMap(1:end-1, 1:end-1), ones(6,6)/36);
            
            figure;
            imagesc(map);
            
            % plot reference visual angles
            xCenter = numel(xbinCenters) / 2 + 0.5;
            yCenter = numel(ybinCenters) / 2 + 0.5;
            theta = 0 : 0.01 : 2*pi;
            
            xScale = abs(xbinCenters(end)-xbinCenters(1))/(numel(xbinCenters)-1);
            yScale = abs(ybinCenters(end)-ybinCenters(1))/(numel(ybinCenters)-1);
                        
            diameterDeg = [4, 8];
            diameter = diameterDeg ./ [xScale, yScale];
            for dd = 1:numel(diameter)
                radius = diameter(dd)/2;
                x = radius * cos(theta) + xCenter;
                y = radius * sin(theta) + yCenter;
                hold on;
                if mod(dd,2) == 0
                    lstyle = '-';
                else 
                    lstyle = '--';
                end
                plot(x, y, ['w' lstyle], 'LineWidth', 3);
                ht = text((radius + 0.2/xScale) * cos(-1*pi()/4) + xCenter, ...
                          (radius + 0.2/yScale) * sin(-1*pi()/4) + yCenter, ...
                          [num2str(diameterDeg(dd)/2) ' deg'], 'color', 'w');
                set(ht, 'Rotation', -45)
            end
            
            % tick locations, in degrees
            xTickLoc = linspace(-5, 5, 5);
            yTickLoc = linspace(-5, 5, 5);
            axis square;
            set(gca, 'xtick', xTickLoc / xScale + xCenter)
            set(gca, 'xticklabel', xTickLoc)
            set(gca, 'ytick', yTickLoc / yScale + yCenter)
            set(gca, 'yticklabel', yTickLoc)
            ylabel('Vertical position (deg)') 
            xlabel('Horizontal position (deg)') 
        end
        
        
        % STATISTICS
        
        % Compute vergence and plot vergence trace
        % [nrg] some things here could be paramaterized and I should
        % include the possibility of plotting data without shadedErrorBar
        function computeCondVergence(this, conds, reference)
            
            matlabVersion = version();
            
            if nargin < 3
                reference = 'none';
            end
            figure; hold on;
            pcolors = get(gca, 'ColorOrder');
            v = cell(numel(conds),1);
            
            mver = version();
            for c = 1 : numel(conds)
                v{c} = this.conditions.(conds{c}).L.H - this.conditions.(conds{c}).R.H;
                
                if strcmp(reference, 'prestim')
                    v{c} = v{c} - nanmean(nanmean(v{c}(:,this.epochTime.(conds{c})<0)));
                end
                
                if matlabVersion(1)=='9'
                    h = shadedErrorBar(this.epochTime.(conds{c}), ...
                                   v{c}, ...
                                   {@(x) nanmean(x,1), @(x) nanstd(x,[],1)/sqrt(size(x, 1))}, ...
                                   'lineProps', {'color', pcolors(c,:)});
                else
                    hold on;
                    warning('The shaded error bar visualization is only available for recent Matlab 2017. Default to simple line plots')
                    mv = mean(v{c}, 1);
                    sdv = std(v{c}, [], 1);
                    semv = sdv/sqrt(size(v{c},1));
                    h = plot(this.epochTime.(conds{c}), ...
                             mv,...
                             'color', pcolors(c,:));
                    hl = plot(this.epochTime.(conds{c}), ...
                             mv - semv, '--',...
                             'color', pcolors(c,:));
                    hu = plot(this.epochTime.(conds{c}), ...
                             mv + semv, '--',...
                             'color', pcolors(c,:));
                end
            
            end
            
            xlabel('Time (ms)')
            ylabel('Vergence (deg)')
            
            % [nrg] so far we are only comparing pairs of conditions
            % comparison of more conditions will be implemented later
            if numel(conds) == 2
                pval = ones(size(v{1}, 2), 1); 
                for i=1:size(v{1}, 2)
                    [p,~] = ranksum(v{1}(:,i), v{2}(:,i));
                    pval(i) = p;
                end
                
            end
            
            alpha = 0.05;
            
            if any(pval<alpha)
                ylimits = get(gca, 'YLim');
                xc = ylimits(2) - (ylimits(2)-ylimits(1))*0.1;
                hold on;
                signiDots = this.epochTime.(conds{c})(pval<alpha);
                plot(signiDots, ones(size(signiDots)) * xc, 'ko', ...
                     'markerfacecolor', 'k');
            end     
            
        end
        
        % Compute saccade statistics
        % [nrg] this needs work. Ideally we would look at number of
        % saccades, saccade velocity/duration, and so on...
        function computeSaccadeStatistics(this, conds)
            
            nSac = zeros(1,numel(conds));
            for c = 1:numel(conds)
                % time
                T = this.conditions.(conds{c}).L.T;

                sac = zeros(size(T));
                for r = 1:numel(this.evnt)  % loop through runs
                    for b = 1:this.evnt{r}.sac.n  % loop through saccades
                        sac(T>=this.evnt{r}.sac.Tstart(b) & T<=this.evnt{r}.sac.Tend(b)) = 1;
                    end
                end
                
                nSac(c) = sum(any(sac, 2))/size(sac, 1);
                
            end
            
            figure;
            bar(1:numel(conds), nSac); 
            xlim([0,numel(conds)+1]); 
            ylim([0,0.2])
            set(gca, 'XTickLabel', conds)
            ylabel('Number of saccades per trial')

        end
        
        
        function v = expectedVergence(ipd, vd)
            % expectedVergence = @(ipd, vd) 2 * 180/pi() * atan(ipd/vd)
            v = 2 * 180/pi() * atan(ipd/vd);
        end
        
        
    end
    
end