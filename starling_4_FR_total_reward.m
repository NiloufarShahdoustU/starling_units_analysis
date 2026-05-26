clc;
clear;
close all;

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};

inputFolder = '\\155.100.91.44\d\Data\Nill\starling\spikes\spike_data\';
OutputFolder = '\\155.100.91.44\d\Code\Nill\Starling_units_analysis\starling_4_FR_total_reward\';

if ~exist(OutputFolder, 'dir')
    mkdir(OutputFolder);
end

% Extract -1 to +1 so baseline normalization is possible.
% Visualization and statistics are ONLY 0 to +1 sec.
windowBeforeSec = 1;
windowAfterSec  = 1;
binSizeSec = 0.05;
smoothBins = 7;

minMeanFRHz = 0.5;
minFracTrialsWithSpikes = 0.10;

nPerm = 1000;
alphaPerm = 0.05;
clusterFormingAlpha = 0.05;

analysisWindowSec = windowBeforeSec + windowAfterSec;

timeEdges = -windowBeforeSec:binSizeSec:windowAfterSec;
timeCenters = timeEdges(1:end-1) + binSizeSec/2;

baselineIdx = timeCenters >= -1 & timeCenters < -0.25;
testIdx = timeCenters >= 0 & timeCenters <= windowAfterSec;

testTimeCenters = timeCenters(testIdx);

winColor  = [0 0.6 0];
loseColor = [0.85 0 0];

for pt = 1:length(microPts)

    ptID = microPts{pt};
    fprintf('\nprocessing patient %s\n', ptID);

    data = fullfile(inputFolder, sprintf('%s_spikeData.mat', ptID));
    load(data);

    ChanUnitTimestamp = spikeData.ChanUnitTimestamp;
    eventTime = spikeData.eventTimes;
    bhvData = eventTime.bhvData;
    inclChans = spikeData.inclChans;
    microLabels = spikeData.microLabels;

    if isfield(eventTime, 'totalRewardTime')
        totalRewardTime = double(eventTime.totalRewardTime);
    elseif isfield(eventTime, 'totalReward')
        totalRewardTime = double(eventTime.totalReward);
    elseif isfield(eventTime, 'totalRewardShowTime')
        totalRewardTime = double(eventTime.totalRewardShowTime);
    elseif isfield(eventTime, 'totalRewardOnsetTime')
        totalRewardTime = double(eventTime.totalRewardOnsetTime);
    else
        error('Could not find total reward time field in spikeData.eventTimes');
    end

    SampleRes = double(spikeData.SampleRes);

    outcome = lower(string(bhvData.outcome));

    winTrials  = find(outcome == "win");
    loseTrials = find(outcome == "lose");

    allChanUnits = unique(ChanUnitTimestamp(:,1:2), 'rows');
    allChanUnits = double(allChanUnits);

    allChanUnits(allChanUnits(:,2) == 255, :) = [];

    maxInclChan = max(double(inclChans(:)));
    allChanUnits(allChanUnits(:,1) > maxInclChan, :) = [];

    for cu = 1:size(allChanUnits, 1)

        chanNum = allChanUnits(cu, 1);
        unitNum = allChanUnits(cu, 2);

        chanIdx = find(double(inclChans) == chanNum, 1);

        if isempty(chanIdx)
            areaName = 'unknown_area';
        else
            areaName = microLabels{chanIdx};
        end

        areaNameClean = regexprep(char(areaName), '[^\w]', '_');

        unitSpikeTimes = ChanUnitTimestamp( ...
            double(ChanUnitTimestamp(:,1)) == chanNum & ...
            double(ChanUnitTimestamp(:,2)) == unitNum, 3);

        unitSpikeTimes = double(unitSpikeTimes);

        if isempty(unitSpikeTimes)
            continue
        end

        winRasterX = [];
        winRasterY = [];
        loseRasterX = [];
        loseRasterY = [];

        winTrialFR = [];
        loseTrialFR = [];

        winTrialSpikeCounts = [];
        loseTrialSpikeCounts = [];

        validWinTrials = 0;
        validLoseTrials = 0;

        for i = 1:length(winTrials)

            tr = winTrials(i);

            if tr > length(totalRewardTime) || isnan(totalRewardTime(tr))
                continue
            end

            validWinTrials = validWinTrials + 1;

            windowStart = totalRewardTime(tr) - windowBeforeSec * SampleRes;
            windowEnd   = totalRewardTime(tr) + windowAfterSec  * SampleRes;

            spikesInWindow = unitSpikeTimes( ...
                unitSpikeTimes >= windowStart & ...
                unitSpikeTimes <= windowEnd);

            relSpikesSec = (spikesInWindow - totalRewardTime(tr)) ./ SampleRes;

            winRasterX = [winRasterX; relSpikesSec(:)];
            winRasterY = [winRasterY; validWinTrials .* ones(length(relSpikesSec), 1)];

            counts = histcounts(relSpikesSec, timeEdges);

            winTrialFR(validWinTrials, :) = counts ./ binSizeSec;
            winTrialSpikeCounts(validWinTrials, 1) = sum(counts);

        end

        for i = 1:length(loseTrials)

            tr = loseTrials(i);

            if tr > length(totalRewardTime) || isnan(totalRewardTime(tr))
                continue
            end

            validLoseTrials = validLoseTrials + 1;

            windowStart = totalRewardTime(tr) - windowBeforeSec * SampleRes;
            windowEnd   = totalRewardTime(tr) + windowAfterSec  * SampleRes;

            spikesInWindow = unitSpikeTimes( ...
                unitSpikeTimes >= windowStart & ...
                unitSpikeTimes <= windowEnd);

            relSpikesSec = (spikesInWindow - totalRewardTime(tr)) ./ SampleRes;

            loseRasterX = [loseRasterX; relSpikesSec(:)];
            loseRasterY = [loseRasterY; validLoseTrials .* ones(length(relSpikesSec), 1)];

            counts = histcounts(relSpikesSec, timeEdges);

            loseTrialFR(validLoseTrials, :) = counts ./ binSizeSec;
            loseTrialSpikeCounts(validLoseTrials, 1) = sum(counts);

        end

        if isempty(winTrialFR) || isempty(loseTrialFR)
            continue
        end

        allTrialSpikeCounts = [winTrialSpikeCounts; loseTrialSpikeCounts];

        nValidTrialsTotal = length(allTrialSpikeCounts);
        totalSpikesInWindow = sum(allTrialSpikeCounts);

        meanFRHz = totalSpikesInWindow / (nValidTrialsTotal * analysisWindowSec);

        fracTrialsWithSpikes = sum(allTrialSpikeCounts > 0) / nValidTrialsTotal;

        if meanFRHz < minMeanFRHz || fracTrialsWithSpikes < minFracTrialsWithSpikes
            continue
        end

        % ================================================================
        % Shared baseline normalization across win and lose
        % baseline = -1 to -0.25 sec before total reward
        % ================================================================
        allTrialFR = [winTrialFR; loseTrialFR];

        baselineVals = allTrialFR(:, baselineIdx);
        baselineMean = mean(baselineVals(:), 'omitnan');
        baselineStd  = std(baselineVals(:), 'omitnan');

        if baselineStd == 0 || isnan(baselineStd)
            baselineStd = 1;
        end

        winZ  = (winTrialFR  - baselineMean) ./ baselineStd;
        loseZ = (loseTrialFR - baselineMean) ./ baselineStd;

        winMean  = mean(winZ, 1, 'omitnan');
        loseMean = mean(loseZ, 1, 'omitnan');

        winSEM  = std(winZ, 0, 1, 'omitnan') ./ sqrt(size(winZ, 1));
        loseSEM = std(loseZ, 0, 1, 'omitnan') ./ sqrt(size(loseZ, 1));

        winMeanSmooth  = smoothdata(winMean,  'gaussian', smoothBins);
        loseMeanSmooth = smoothdata(loseMean, 'gaussian', smoothBins);

        winSEMSmooth  = smoothdata(winSEM,  'gaussian', smoothBins);
        loseSEMSmooth = smoothdata(loseSEM, 'gaussian', smoothBins);

        % ================================================================
        % Non-parametric cluster-based permutation test
        % win vs lose, ONLY 0 to +1 sec
        % Uses ranksum per time bin.
        % Cluster mass = sum(abs(z-values)) inside contiguous significant bins.
        % Correction = compare observed cluster mass to max permuted cluster mass.
        % ================================================================

        sigSegments = [];

        winTestZ  = winZ(:, testIdx);
        loseTestZ = loseZ(:, testIdx);

        nTestBins = size(winTestZ, 2);

        observedP = nan(1, nTestBins);
        observedZ = nan(1, nTestBins);

        for bb = 1:nTestBins

            x = winTestZ(:, bb);
            y = loseTestZ(:, bb);

            x = x(~isnan(x));
            y = y(~isnan(y));

            if length(x) < 2 || length(y) < 2
                continue
            end

            [pVal, ~, stats] = ranksum(x, y);

            observedP(bb) = pVal;

            if isfield(stats, 'zval')
                observedZ(bb) = stats.zval;
            else
                observedZ(bb) = 0;
            end

        end

        observedSigBins = observedP < clusterFormingAlpha;

        observedClusterMasses = [];
        observedClusterStarts = [];
        observedClusterEnds = [];

        bb = 1;
        while bb <= nTestBins

            if observedSigBins(bb)

                clusterStartIdx = bb;

                while bb <= nTestBins && observedSigBins(bb)
                    bb = bb + 1;
                end

                clusterEndIdx = bb - 1;

                clusterMass = sum(abs(observedZ(clusterStartIdx:clusterEndIdx)), 'omitnan');

                observedClusterMasses = [observedClusterMasses; clusterMass];
                observedClusterStarts = [observedClusterStarts; testTimeCenters(clusterStartIdx) - binSizeSec/2];
                observedClusterEnds   = [observedClusterEnds;   testTimeCenters(clusterEndIdx)   + binSizeSec/2];

            else
                bb = bb + 1;
            end

        end

        allTestZ = [winTestZ; loseTestZ];

        nWinHere = size(winTestZ, 1);
        nLoseHere = size(loseTestZ, 1);
        nTotalHere = nWinHere + nLoseHere;

        maxPermClusterMass = zeros(nPerm, 1);

        for pp = 1:nPerm

            shuffledIdx = randperm(nTotalHere);

            permWinIdx  = shuffledIdx(1:nWinHere);
            permLoseIdx = shuffledIdx(nWinHere+1:end);

            permWinZ  = allTestZ(permWinIdx, :);
            permLoseZ = allTestZ(permLoseIdx, :);

            permP = nan(1, nTestBins);
            permZ = nan(1, nTestBins);

            for bb = 1:nTestBins

                x = permWinZ(:, bb);
                y = permLoseZ(:, bb);

                x = x(~isnan(x));
                y = y(~isnan(y));

                if length(x) < 2 || length(y) < 2
                    continue
                end

                [pVal, ~, stats] = ranksum(x, y);

                permP(bb) = pVal;

                if isfield(stats, 'zval')
                    permZ(bb) = stats.zval;
                else
                    permZ(bb) = 0;
                end

            end

            permSigBins = permP < clusterFormingAlpha;

            permClusterMasses = [];

            bb = 1;
            while bb <= nTestBins

                if permSigBins(bb)

                    clusterStartIdx = bb;

                    while bb <= nTestBins && permSigBins(bb)
                        bb = bb + 1;
                    end

                    clusterEndIdx = bb - 1;

                    clusterMass = sum(abs(permZ(clusterStartIdx:clusterEndIdx)), 'omitnan');

                    permClusterMasses = [permClusterMasses; clusterMass];

                else
                    bb = bb + 1;
                end

            end

            if ~isempty(permClusterMasses)
                maxPermClusterMass(pp) = max(permClusterMasses);
            else
                maxPermClusterMass(pp) = 0;
            end

        end

        clusterPvals = [];

        for cc = 1:length(observedClusterMasses)

            clusterP = (sum(maxPermClusterMass >= observedClusterMasses(cc)) + 1) / ...
                       (nPerm + 1);

            clusterPvals = [clusterPvals; clusterP];

            if clusterP < alphaPerm
                sigSegments = [sigSegments; observedClusterStarts(cc) observedClusterEnds(cc)];
            end

        end

        if ~isempty(sigSegments)
            sigSegments(sigSegments(:,2) <= 0, :) = [];
            sigSegments(:,1) = max(sigSegments(:,1), 0);
            sigSegments(:,2) = min(sigSegments(:,2), windowAfterSec);
        end

        % ================================================================
        % Plot
        % ================================================================

        fig = figure('Visible', 'off', 'Color', 'w');
        set(fig, 'Position', [100 100 900 700]);

        subplot(2,1,1)
        hold on

        scatter(winRasterX, winRasterY, 6, ...
            'filled', ...
            'MarkerFaceColor', winColor, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.4)

        scatter(loseRasterX, loseRasterY + validWinTrials, 6, ...
            'filled', ...
            'MarkerFaceColor', loseColor, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.4)

        xline(0, '--k')
        yline(validWinTrials + 0.5, '--k')

        xlim([0 windowAfterSec])
        xlabel('time from total reward onset (s)')
        ylabel('trials')

        title(sprintf('%s | %s | Ch %d Unit %d | Raster | FR %.2f Hz', ...
            ptID, areaNameClean, chanNum, unitNum, meanFRHz), ...
            'Interpreter', 'none')

        legend({'win', 'lose'}, 'Location', 'best')
        box off

        subplot(2,1,2)
        hold on

        plotIdx = timeCenters >= 0 & timeCenters <= windowAfterSec;

        fill([timeCenters(plotIdx) fliplr(timeCenters(plotIdx))], ...
            [winMeanSmooth(plotIdx) + winSEMSmooth(plotIdx), ...
             fliplr(winMeanSmooth(plotIdx) - winSEMSmooth(plotIdx))], ...
            winColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

        fill([timeCenters(plotIdx) fliplr(timeCenters(plotIdx))], ...
            [loseMeanSmooth(plotIdx) + loseSEMSmooth(plotIdx), ...
             fliplr(loseMeanSmooth(plotIdx) - loseSEMSmooth(plotIdx))], ...
            loseColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

        plot(timeCenters(plotIdx), winMeanSmooth(plotIdx), ...
            'Color', winColor, 'LineWidth', 1.8)

        plot(timeCenters(plotIdx), loseMeanSmooth(plotIdx), ...
            'Color', loseColor, 'LineWidth', 1.8)

        xline(0, '--k')
        yline(0, ':k')

        upperVals = [winMeanSmooth(plotIdx) + winSEMSmooth(plotIdx), ...
                     loseMeanSmooth(plotIdx) + loseSEMSmooth(plotIdx)];

        lowerVals = [winMeanSmooth(plotIdx) - winSEMSmooth(plotIdx), ...
                     loseMeanSmooth(plotIdx) - loseSEMSmooth(plotIdx)];

        curveHigh = max(upperVals, [], 'omitnan');
        curveLow  = min(lowerVals, [], 'omitnan');

        if isempty(curveHigh) || isnan(curveHigh)
            curveHigh = 1;
        end

        if isempty(curveLow) || isnan(curveLow)
            curveLow = -1;
        end

        yRange = curveHigh - curveLow;

        if yRange == 0 || isnan(yRange)
            yRange = 1;
        end

        sigY = curveHigh + 0.15 * yRange;

        for ss = 1:size(sigSegments, 1)

            sigStart = sigSegments(ss, 1);
            sigEnd   = sigSegments(ss, 2);

            plot([sigStart sigEnd], [sigY sigY], ...
                '-', ...
                'Color', [0.5 0.5 0.5], ...
                'LineWidth', 4)

            labelY = sigY + 0.04 * yRange;

            text(sigStart, labelY, sprintf('%.0f ms', sigStart * 1000), ...
                'Color', [0.35 0.35 0.35], ...
                'FontSize', 8, ...
                'Rotation', 90, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle')

            text(sigEnd, labelY, sprintf('%.0f ms', sigEnd * 1000), ...
                'Color', [0.35 0.35 0.35], ...
                'FontSize', 8, ...
                'Rotation', 90, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle')

        end

        if ~isempty(sigSegments)
            ylim([curveLow - 0.15 * yRange, sigY + 0.55 * yRange])
        else
            ylim([curveLow - 0.15 * yRange, curveHigh + 0.25 * yRange])
        end

        xlim([0 windowAfterSec])
        xlabel('time from total reward onset (s)')
        ylabel('baseline z-scored firing rate')
        title('PSTH: mean ± SEM, cluster-corrected ranksum test');
        box off

        pdfName = sprintf('%s_%s_ch%d_unit%d_totalReward_clusterPerm.pdf', ...
            ptID, areaNameClean, chanNum, unitNum);

        pdfPath = fullfile(OutputFolder, pdfName);

        exportgraphics(fig, pdfPath, 'ContentType', 'vector');
        close(fig);

        fprintf('saved: %s \n', pdfName);

    end
end