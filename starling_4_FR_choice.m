clc;
clear;
close all;

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};

inputFolder = '\\155.100.91.44\d\Data\Nill\starling\spikes\spike_data\';
OutputFolder = '\\155.100.91.44\d\Code\Nill\Starling_units_analysis\starling_4_FR_choice\';

if ~exist(OutputFolder, 'dir')
    mkdir(OutputFolder);
end

windowBeforeSec = 3;
windowAfterSec  = 1;
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
baselineIdx = timeCenters >= -1 & timeCenters < -.25;

arrowUpColor   = [0 0.45 0.85];
arrowDownColor = [0.85 0.35 0];

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

    choiceTime = double(spikeData.eventTimes.choiceAndFeedbackTime);
    SampleRes  = double(spikeData.SampleRes);

    choice = lower(string(bhvData.choice));

    arrowUpTrials   = find(choice == "arrowup");
    arrowDownTrials = find(choice == "arrowdown");

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

        arrowUpRasterX = [];
        arrowUpRasterY = [];
        arrowDownRasterX = [];
        arrowDownRasterY = [];

        arrowUpTrialFR = [];
        arrowDownTrialFR = [];

        arrowUpTrialSpikeCounts = [];
        arrowDownTrialSpikeCounts = [];

        validArrowUpTrials = 0;
        validArrowDownTrials = 0;

        for i = 1:length(arrowUpTrials)

            tr = arrowUpTrials(i);

            if tr > length(choiceTime) || isnan(choiceTime(tr))
                continue
            end

            validArrowUpTrials = validArrowUpTrials + 1;

            windowStart = choiceTime(tr) - (windowBeforeSec * SampleRes);
            windowEnd   = choiceTime(tr) + (windowAfterSec  * SampleRes);

            spikesInWindow = unitSpikeTimes( ...
                unitSpikeTimes >= windowStart & ...
                unitSpikeTimes <= windowEnd);

            relSpikesSec = (spikesInWindow - choiceTime(tr)) ./ SampleRes;

            arrowUpRasterX = [arrowUpRasterX; relSpikesSec(:)];
            arrowUpRasterY = [arrowUpRasterY; validArrowUpTrials .* ones(length(relSpikesSec), 1)];

            counts = histcounts(relSpikesSec, timeEdges);

            arrowUpTrialFR(validArrowUpTrials, :) = counts ./ binSizeSec;
            arrowUpTrialSpikeCounts(validArrowUpTrials, 1) = sum(counts);

        end

        for i = 1:length(arrowDownTrials)

            tr = arrowDownTrials(i);

            if tr > length(choiceTime) || isnan(choiceTime(tr))
                continue
            end

            validArrowDownTrials = validArrowDownTrials + 1;

            windowStart = choiceTime(tr) - (windowBeforeSec * SampleRes);
            windowEnd   = choiceTime(tr) + (windowAfterSec  * SampleRes);

            spikesInWindow = unitSpikeTimes( ...
                unitSpikeTimes >= windowStart & ...
                unitSpikeTimes <= windowEnd);

            relSpikesSec = (spikesInWindow - choiceTime(tr)) ./ SampleRes;

            arrowDownRasterX = [arrowDownRasterX; relSpikesSec(:)];
            arrowDownRasterY = [arrowDownRasterY; validArrowDownTrials .* ones(length(relSpikesSec), 1)];

            counts = histcounts(relSpikesSec, timeEdges);

            arrowDownTrialFR(validArrowDownTrials, :) = counts ./ binSizeSec;
            arrowDownTrialSpikeCounts(validArrowDownTrials, 1) = sum(counts);

        end

        if isempty(arrowUpTrialFR) && isempty(arrowDownTrialFR)
            continue
        end

        allTrialSpikeCounts = [arrowUpTrialSpikeCounts; arrowDownTrialSpikeCounts];

        nValidTrialsTotal = length(allTrialSpikeCounts);
        totalSpikesInWindow = sum(allTrialSpikeCounts);

        meanFRHz = totalSpikesInWindow / (nValidTrialsTotal * analysisWindowSec);

        fracTrialsWithSpikes = sum(allTrialSpikeCounts > 0) / nValidTrialsTotal;

        if meanFRHz < minMeanFRHz || fracTrialsWithSpikes < minFracTrialsWithSpikes
            continue
        end

        allTrialFR = [arrowUpTrialFR; arrowDownTrialFR];

        baselineVals = allTrialFR(:, baselineIdx);
        baselineMean = mean(baselineVals(:), 'omitnan');
        baselineStd  = std(baselineVals(:), 'omitnan');

        if baselineStd == 0 || isnan(baselineStd)
            baselineStd = 1;
        end

        arrowUpZ   = (arrowUpTrialFR   - baselineMean) ./ baselineStd;
        arrowDownZ = (arrowDownTrialFR - baselineMean) ./ baselineStd;

        arrowUpMean   = mean(arrowUpZ, 1, 'omitnan');
        arrowDownMean = mean(arrowDownZ, 1, 'omitnan');

        arrowUpSEM   = std(arrowUpZ, 0, 1, 'omitnan') ./ sqrt(size(arrowUpZ, 1));
        arrowDownSEM = std(arrowDownZ, 0, 1, 'omitnan') ./ sqrt(size(arrowDownZ, 1));

        arrowUpMeanSmooth   = smoothdata(arrowUpMean,   'gaussian', smoothBins);
        arrowDownMeanSmooth = smoothdata(arrowDownMean, 'gaussian', smoothBins);

        arrowUpSEMSmooth   = smoothdata(arrowUpSEM,   'gaussian', smoothBins);
        arrowDownSEMSmooth = smoothdata(arrowDownSEM, 'gaussian', smoothBins);

        sigSegments = [];

        if ~isempty(arrowUpZ) && ~isempty(arrowDownZ)

            % ONLY test the -3 to 0 sec pre-choice window
            permStartTimes = -windowBeforeSec:permStrideSec:(0 - permWindowSec);
            permEndTimes = permStartTimes + permWindowSec;

            nWindows = length(permStartTimes);

            arrowUpWindowMat   = nan(size(arrowUpZ, 1), nWindows);
            arrowDownWindowMat = nan(size(arrowDownZ, 1), nWindows);

            for ww = 1:nWindows

                thisStart = permStartTimes(ww);
                thisEnd   = permEndTimes(ww);

                thisIdx = timeCenters >= thisStart & timeCenters < thisEnd;

                if sum(thisIdx) < 1
                    continue
                end

                arrowUpWindowMat(:, ww)   = mean(arrowUpZ(:, thisIdx), 2, 'omitnan');
                arrowDownWindowMat(:, ww) = mean(arrowDownZ(:, thisIdx), 2, 'omitnan');

            end

            observedDiffs = nan(1, nWindows);

            for ww = 1:nWindows
                observedDiffs(ww) = mean(arrowUpWindowMat(:, ww), 'omitnan') - ...
                                    mean(arrowDownWindowMat(:, ww), 'omitnan');
            end

            allWindowMat = [arrowUpWindowMat; arrowDownWindowMat];

            nArrowUpHere   = size(arrowUpWindowMat, 1);
            nArrowDownHere = size(arrowDownWindowMat, 1);
            nTotalHere = nArrowUpHere + nArrowDownHere;

            permDiffs = nan(nPerm, nWindows);

            for pp = 1:nPerm

                shuffledIdx = randperm(nTotalHere);

                permArrowUpIdx   = shuffledIdx(1:nArrowUpHere);
                permArrowDownIdx = shuffledIdx(nArrowUpHere+1:end);

                permArrowUpMat   = allWindowMat(permArrowUpIdx, :);
                permArrowDownMat = allWindowMat(permArrowDownIdx, :);

                for ww = 1:nWindows
                    permDiffs(pp, ww) = mean(permArrowUpMat(:, ww), 'omitnan') - ...
                                        mean(permArrowDownMat(:, ww), 'omitnan');
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
            % Keep ONLY significant clusters fully/partly in -3 to 0 sec
            if ~isempty(sigSegments)
                sigSegments(sigSegments(:,1) >= 0, :) = [];
                sigSegments(:,2) = min(sigSegments(:,2), 0);
                sigSegments(:,1) = max(sigSegments(:,1), -windowBeforeSec);
            end
        end

        fig = figure('Visible', 'off', 'Color', 'w');
        set(fig, 'Position', [100 100 900 700]);

        subplot(2,1,1)
        hold on

        scatter(arrowUpRasterX, arrowUpRasterY, 6, ...
            'filled', ...
            'MarkerFaceColor', arrowUpColor, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.4)

        scatter(arrowDownRasterX, arrowDownRasterY + validArrowUpTrials, 6, ...
            'filled', ...
            'MarkerFaceColor', arrowDownColor, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.4)

        xline(0, '--k')
        yline(validArrowUpTrials + 0.5, '--k')

        xlim([-windowBeforeSec windowAfterSec])
        xlabel('time from choice onset (s)')
        ylabel('trials')

        title(sprintf('%s | %s | Ch %d Unit %d | Raster | FR %.2f Hz', ...
            ptID, areaNameClean, chanNum, unitNum, meanFRHz), ...
            'Interpreter', 'none')

        legend({'arrowup', 'arrowdown'}, 'Location', 'best')
        box off

        subplot(2,1,2)
        hold on

        fill([timeCenters fliplr(timeCenters)], ...
            [arrowUpMeanSmooth + arrowUpSEMSmooth fliplr(arrowUpMeanSmooth - arrowUpSEMSmooth)], ...
            arrowUpColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

        fill([timeCenters fliplr(timeCenters)], ...
            [arrowDownMeanSmooth + arrowDownSEMSmooth fliplr(arrowDownMeanSmooth - arrowDownSEMSmooth)], ...
            arrowDownColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

        plot(timeCenters, arrowUpMeanSmooth, ...
            'Color', arrowUpColor, 'LineWidth', 1.8)

        plot(timeCenters, arrowDownMeanSmooth, ...
            'Color', arrowDownColor, 'LineWidth', 1.8)

        xline(0, '--k')
        yline(0, ':k')

        upperVals = [arrowUpMeanSmooth + arrowUpSEMSmooth, arrowDownMeanSmooth + arrowDownSEMSmooth];
        lowerVals = [arrowUpMeanSmooth - arrowUpSEMSmooth, arrowDownMeanSmooth - arrowDownSEMSmooth];

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
        end

        xlim([-windowBeforeSec windowAfterSec])
        xlabel('time from choice onset (s)')
        ylabel('baseline z-scored firing rate')
        title('PSTH: mean ± sem');
        box off

        pdfName = sprintf('%s_%s_ch%d_unit%d_choice.pdf', ...
            ptID, areaNameClean, chanNum, unitNum);

        pdfPath = fullfile(OutputFolder, pdfName);

        exportgraphics(fig, pdfPath, 'ContentType', 'vector');
        close(fig);

        fprintf('saved: %s \n', pdfName);

    end
end