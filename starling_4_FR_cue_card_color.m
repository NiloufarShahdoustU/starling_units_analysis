clc;
clear;
close all;

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};

inputFolder = '\\155.100.91.44\d\Data\Nill\starling\spikes\spike_data\';
OutputFolder = '\\155.100.91.44\d\Code\Nill\Starling_units_analysis\starling_4_FR_cue_card_color\';

if ~exist(OutputFolder, 'dir')
    mkdir(OutputFolder);
end

windowBeforeSec = 3;      % needed only for baseline
windowAfterSec  = 1;      % show/test cardShowTime to +1 sec
binSizeSec = 0.05;
smoothBins = 7;

minMeanFRHz = 0.5;
minFracTrialsWithSpikes = 0.10;

permWindowMs = 200;
permStrideMs = 50;
nPerm = 1000;
alphaPerm = 0.05;

permWindowSec = permWindowMs / 1000;
permStrideSec = permStrideMs / 1000;

analysisWindowSec = windowBeforeSec + windowAfterSec;

timeEdges = -windowBeforeSec:binSizeSec:windowAfterSec;
timeCenters = timeEdges(1:end-1) + binSizeSec/2;

baselineIdx = timeCenters >= -1 & timeCenters < -0.25;
plotIdx = timeCenters >= 0 & timeCenters <= 1;

uniformColor =  [0.5 0.5 0.5];
lowColor     = [0.85 0.35 0];
highColor    = [0.25 0.65 0.25];

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

    cardShowTime = double(eventTime.cardShowTime);
    SampleRes = double(spikeData.SampleRes);

    deckType = lower(string(bhvData.distribution));

    uniformTrials = find(deckType == "uniform");
    lowTrials     = find(deckType == "low");
    highTrials    = find(deckType == "high");
    
    
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

        [uniformRasterX, uniformRasterY, uniformTrialFR, uniformTrialSpikeCounts, validUniformTrials] = ...
            getTrialFR(uniformTrials, cardShowTime, unitSpikeTimes, SampleRes, ...
            windowBeforeSec, windowAfterSec, timeEdges, binSizeSec);

        [lowRasterX, lowRasterY, lowTrialFR, lowTrialSpikeCounts, validLowTrials] = ...
            getTrialFR(lowTrials, cardShowTime, unitSpikeTimes, SampleRes, ...
            windowBeforeSec, windowAfterSec, timeEdges, binSizeSec);

        [highRasterX, highRasterY, highTrialFR, highTrialSpikeCounts, validHighTrials] = ...
            getTrialFR(highTrials, cardShowTime, unitSpikeTimes, SampleRes, ...
            windowBeforeSec, windowAfterSec, timeEdges, binSizeSec);

        if isempty(uniformTrialFR) && isempty(lowTrialFR) && isempty(highTrialFR)
            continue
        end

        allTrialSpikeCounts = [uniformTrialSpikeCounts; lowTrialSpikeCounts; highTrialSpikeCounts];

        nValidTrialsTotal = length(allTrialSpikeCounts);
        totalSpikesInWindow = sum(allTrialSpikeCounts);

        meanFRHz = totalSpikesInWindow / (nValidTrialsTotal * analysisWindowSec);
        fracTrialsWithSpikes = sum(allTrialSpikeCounts > 0) / nValidTrialsTotal;

        if meanFRHz < minMeanFRHz || fracTrialsWithSpikes < minFracTrialsWithSpikes
            continue
        end

        allTrialFR = [uniformTrialFR; lowTrialFR; highTrialFR];

        baselineVals = allTrialFR(:, baselineIdx);
        baselineMean = mean(baselineVals(:), 'omitnan');
        baselineStd  = std(baselineVals(:), 'omitnan');

        if baselineStd == 0 || isnan(baselineStd)
            baselineStd = 1;
        end

        uniformZ = (uniformTrialFR - baselineMean) ./ baselineStd;
        lowZ     = (lowTrialFR     - baselineMean) ./ baselineStd;
        highZ    = (highTrialFR    - baselineMean) ./ baselineStd;

        uniformMean = mean(uniformZ, 1, 'omitnan');
        lowMean     = mean(lowZ, 1, 'omitnan');
        highMean    = mean(highZ, 1, 'omitnan');

        uniformSEM = std(uniformZ, 0, 1, 'omitnan') ./ sqrt(size(uniformZ, 1));
        lowSEM     = std(lowZ, 0, 1, 'omitnan') ./ sqrt(size(lowZ, 1));
        highSEM    = std(highZ, 0, 1, 'omitnan') ./ sqrt(size(highZ, 1));

        uniformMeanSmooth = smoothdata(uniformMean, 'gaussian', smoothBins);
        lowMeanSmooth     = smoothdata(lowMean,     'gaussian', smoothBins);
        highMeanSmooth    = smoothdata(highMean,    'gaussian', smoothBins);

        uniformSEMSmooth = smoothdata(uniformSEM, 'gaussian', smoothBins);
        lowSEMSmooth     = smoothdata(lowSEM,     'gaussian', smoothBins);
        highSEMSmooth    = smoothdata(highSEM,    'gaussian', smoothBins);

        sigSegments_uniform_low = runClusterPermutationTwoGroups( ...
            uniformZ, lowZ, timeCenters, ...
            permWindowSec, permStrideSec, nPerm, alphaPerm);

        sigSegments_uniform_high = runClusterPermutationTwoGroups( ...
            uniformZ, highZ, timeCenters, ...
            permWindowSec, permStrideSec, nPerm, alphaPerm);

        sigSegments_low_high = runClusterPermutationTwoGroups( ...
            lowZ, highZ, timeCenters, ...
            permWindowSec, permStrideSec, nPerm, alphaPerm);

        fig = figure('Visible', 'off', 'Color', 'w');
        set(fig, 'Position', [100 100 950 750]);

        subplot(2,1,1)
        hold on

        scatter(uniformRasterX, uniformRasterY, 6, ...
            'filled', ...
            'MarkerFaceColor', uniformColor, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.4)

        scatter(lowRasterX, lowRasterY + validUniformTrials, 6, ...
            'filled', ...
            'MarkerFaceColor', lowColor, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.4)

        scatter(highRasterX, highRasterY + validUniformTrials + validLowTrials, 6, ...
            'filled', ...
            'MarkerFaceColor', highColor, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.4)

        xline(0, '--k')
        yline(validUniformTrials + 0.5, '--k')
        yline(validUniformTrials + validLowTrials + 0.5, '--k')

        xlim([0 1])
        xlabel('time from card show onset (s)')
        ylabel('trials')

        title(sprintf('%s | %s | Ch %d Unit %d | Raster | FR %.2f Hz', ...
            ptID, areaNameClean, chanNum, unitNum, meanFRHz), ...
            'Interpreter', 'none')

        legend({'uniform', 'low', 'high'}, 'Location', 'best')
        box off

        subplot(2,1,2)
        hold on

        fill([timeCenters(plotIdx) fliplr(timeCenters(plotIdx))], ...
            [uniformMeanSmooth(plotIdx) + uniformSEMSmooth(plotIdx), ...
            fliplr(uniformMeanSmooth(plotIdx) - uniformSEMSmooth(plotIdx))], ...
            uniformColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

        fill([timeCenters(plotIdx) fliplr(timeCenters(plotIdx))], ...
            [lowMeanSmooth(plotIdx) + lowSEMSmooth(plotIdx), ...
            fliplr(lowMeanSmooth(plotIdx) - lowSEMSmooth(plotIdx))], ...
            lowColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

        fill([timeCenters(plotIdx) fliplr(timeCenters(plotIdx))], ...
            [highMeanSmooth(plotIdx) + highSEMSmooth(plotIdx), ...
            fliplr(highMeanSmooth(plotIdx) - highSEMSmooth(plotIdx))], ...
            highColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

        plot(timeCenters(plotIdx), uniformMeanSmooth(plotIdx), ...
            'Color', uniformColor, 'LineWidth', 1.8)

        plot(timeCenters(plotIdx), lowMeanSmooth(plotIdx), ...
            'Color', lowColor, 'LineWidth', 1.8)

        plot(timeCenters(plotIdx), highMeanSmooth(plotIdx), ...
            'Color', highColor, 'LineWidth', 1.8)

        xline(0, '--k')
        yline(0, ':k')

        upperVals = [uniformMeanSmooth(plotIdx) + uniformSEMSmooth(plotIdx), ...
                     lowMeanSmooth(plotIdx) + lowSEMSmooth(plotIdx), ...
                     highMeanSmooth(plotIdx) + highSEMSmooth(plotIdx)];

        lowerVals = [uniformMeanSmooth(plotIdx) - uniformSEMSmooth(plotIdx), ...
                     lowMeanSmooth(plotIdx) - lowSEMSmooth(plotIdx), ...
                     highMeanSmooth(plotIdx) - highSEMSmooth(plotIdx)];

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

        sigY1 = curveHigh + 0.15 * yRange;
        sigY2 = curveHigh + 0.30 * yRange;
        sigY3 = curveHigh + 0.45 * yRange;

        for ss = 1:size(sigSegments_uniform_low, 1)
            plot([sigSegments_uniform_low(ss,1) sigSegments_uniform_low(ss,2)], ...
                 [sigY1 sigY1], '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 4)
            text(sigSegments_uniform_low(ss,1), sigY1, 'U-L', ...
                'Color', [0.35 0.35 0.35], 'FontSize', 8, ...
                'VerticalAlignment', 'bottom')
        end

        for ss = 1:size(sigSegments_uniform_high, 1)
            plot([sigSegments_uniform_high(ss,1) sigSegments_uniform_high(ss,2)], ...
                 [sigY2 sigY2], '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 4)
            text(sigSegments_uniform_high(ss,1), sigY2, 'U-H', ...
                'Color', [0.15 0.15 0.15], 'FontSize', 8, ...
                'VerticalAlignment', 'bottom')
        end

        for ss = 1:size(sigSegments_low_high, 1)
            plot([sigSegments_low_high(ss,1) sigSegments_low_high(ss,2)], ...
                 [sigY3 sigY3], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 4)
            text(sigSegments_low_high(ss,1), sigY3, 'L-H', ...
                'Color', [0.45 0.45 0.45], 'FontSize', 8, ...
                'VerticalAlignment', 'bottom')
        end

        if ~isempty(sigSegments_uniform_low) || ~isempty(sigSegments_uniform_high) || ~isempty(sigSegments_low_high)
            ylim([curveLow - 0.15 * yRange, sigY3 + 0.35 * yRange])
        end

        xlim([0 1])
        xlabel('time from card show onset (s)')
        ylabel('baseline z-scored firing rate')
        title('Card show PSTH: 0 to 1 sec | cluster-corrected permutation tests');

        box off

        pdfName = sprintf('%s_%s_ch%d_unit%d_cardShow_uniform_low_high.pdf', ...
            ptID, areaNameClean, chanNum, unitNum);

        pdfPath = fullfile(OutputFolder, pdfName);

        exportgraphics(fig, pdfPath, 'ContentType', 'vector');
        close(fig);

        fprintf('saved: %s \n', pdfName);

    end
end


function [rasterX, rasterY, trialFR, trialSpikeCounts, validTrials] = getTrialFR( ...
    trials, alignTime, unitSpikeTimes, SampleRes, ...
    windowBeforeSec, windowAfterSec, timeEdges, binSizeSec)

    rasterX = [];
    rasterY = [];
    trialFR = [];
    trialSpikeCounts = [];
    validTrials = 0;

    for i = 1:length(trials)

        tr = trials(i);

        if tr > length(alignTime) || isnan(alignTime(tr))
            continue
        end

        validTrials = validTrials + 1;

        windowStart = alignTime(tr) - (windowBeforeSec * SampleRes);
        windowEnd   = alignTime(tr) + (windowAfterSec  * SampleRes);

        spikesInWindow = unitSpikeTimes( ...
            unitSpikeTimes >= windowStart & ...
            unitSpikeTimes <= windowEnd);

        relSpikesSec = (spikesInWindow - alignTime(tr)) ./ SampleRes;

        rasterX = [rasterX; relSpikesSec(:)];
        rasterY = [rasterY; validTrials .* ones(length(relSpikesSec), 1)];

        counts = histcounts(relSpikesSec, timeEdges);

        trialFR(validTrials, :) = counts ./ binSizeSec;
        trialSpikeCounts(validTrials, 1) = sum(counts);

    end
end


function sigSegments = runClusterPermutationTwoGroups( ...
    group1Z, group2Z, timeCenters, ...
    permWindowSec, permStrideSec, nPerm, alphaPerm)

    sigSegments = [];

    if isempty(group1Z) || isempty(group2Z)
        return
    end

    if size(group1Z, 1) < 2 || size(group2Z, 1) < 2
        return
    end

    testStartSec = 0;
    testEndSec = 1;

    permStartTimes = testStartSec:permStrideSec:(testEndSec - permWindowSec);
    permEndTimes = permStartTimes + permWindowSec;

    nWindows = length(permStartTimes);

    group1WindowMat = nan(size(group1Z, 1), nWindows);
    group2WindowMat = nan(size(group2Z, 1), nWindows);

    for ww = 1:nWindows

        thisStart = permStartTimes(ww);
        thisEnd   = permEndTimes(ww);

        thisIdx = timeCenters >= thisStart & timeCenters < thisEnd;

        if sum(thisIdx) < 1
            continue
        end

        group1WindowMat(:, ww) = mean(group1Z(:, thisIdx), 2, 'omitnan');
        group2WindowMat(:, ww) = mean(group2Z(:, thisIdx), 2, 'omitnan');

    end

    observedDiffs = nan(1, nWindows);

    for ww = 1:nWindows
        observedDiffs(ww) = mean(group1WindowMat(:, ww), 'omitnan') - ...
                            mean(group2WindowMat(:, ww), 'omitnan');
    end

    allWindowMat = [group1WindowMat; group2WindowMat];

    nGroup1Here = size(group1WindowMat, 1);
    nTotalHere = size(allWindowMat, 1);

    permDiffs = nan(nPerm, nWindows);

    for pp = 1:nPerm

        shuffledIdx = randperm(nTotalHere);

        permGroup1Idx = shuffledIdx(1:nGroup1Here);
        permGroup2Idx = shuffledIdx(nGroup1Here+1:end);

        permGroup1Mat = allWindowMat(permGroup1Idx, :);
        permGroup2Mat = allWindowMat(permGroup2Idx, :);

        for ww = 1:nWindows
            permDiffs(pp, ww) = mean(permGroup1Mat(:, ww), 'omitnan') - ...
                                mean(permGroup2Mat(:, ww), 'omitnan');
        end

    end

    permThresholds = prctile(abs(permDiffs), 100 * (1 - alphaPerm), 1);

    sigWindowIdx = abs(observedDiffs) > permThresholds;

    observedClusterMasses = [];
    observedClusterStarts = [];
    observedClusterEnds = [];

    ww = 1;
    while ww <= nWindows

        if sigWindowIdx(ww)

            clusterStartIdx = ww;

            while ww <= nWindows && sigWindowIdx(ww)
                ww = ww + 1;
            end

            clusterEndIdx = ww - 1;

            clusterMass = sum(abs(observedDiffs(clusterStartIdx:clusterEndIdx)), 'omitnan');

            observedClusterMasses = [observedClusterMasses; clusterMass];
            observedClusterStarts = [observedClusterStarts; permStartTimes(clusterStartIdx)];
            observedClusterEnds   = [observedClusterEnds; permEndTimes(clusterEndIdx)];

        else
            ww = ww + 1;
        end
    end

    maxPermClusterMass = zeros(nPerm, 1);

    for pp = 1:nPerm

        permSigWindowIdx = abs(permDiffs(pp, :)) > permThresholds;

        permClusterMasses = [];

        ww = 1;
        while ww <= nWindows

            if permSigWindowIdx(ww)

                clusterStartIdx = ww;

                while ww <= nWindows && permSigWindowIdx(ww)
                    ww = ww + 1;
                end

                clusterEndIdx = ww - 1;

                clusterMass = sum(abs(permDiffs(pp, clusterStartIdx:clusterEndIdx)), 'omitnan');

                permClusterMasses = [permClusterMasses; clusterMass];

            else
                ww = ww + 1;
            end
        end

        if ~isempty(permClusterMasses)
            maxPermClusterMass(pp) = max(permClusterMasses);
        end
    end

    for cc = 1:length(observedClusterMasses)

        clusterP = (sum(maxPermClusterMass >= observedClusterMasses(cc)) + 1) / ...
                   (nPerm + 1);

        if clusterP < alphaPerm
            sigSegments = [sigSegments; observedClusterStarts(cc) observedClusterEnds(cc)];
        end
    end

    if ~isempty(sigSegments)
        sigSegments(sigSegments(:,2) <= 0, :) = [];
        sigSegments(sigSegments(:,1) >= 1, :) = [];
        sigSegments(:,1) = max(sigSegments(:,1), 0);
        sigSegments(:,2) = min(sigSegments(:,2), 1);
    end
end