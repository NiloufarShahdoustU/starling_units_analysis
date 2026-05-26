clc;
clear;
close all;

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};

inputFolder = '\\155.100.91.44\d\Data\Nill\starling\spikes\spike_data\';
OutputFolder = '\\155.100.91.44\d\Code\Nill\Starling_units_analysis\starling_4_FR_flip_card_color\';

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
baselineIdx = timeCenters >= -1 & timeCenters < -0.25;

uniformColor = [0.5 0.5 0.5];
lowColor     = [0.85 0.35 0];
highColor    = [0.25 0.65 0.25];

allStats = table();
allClusterStats = table();

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

    distribution = lower(string(bhvData.distribution));

    uniformTrials = find(distribution == "uniform");
    lowTrials     = find(distribution == "low");
    highTrials    = find(distribution == "high");

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

        uniformRasterX = [];
        uniformRasterY = [];
        lowRasterX = [];
        lowRasterY = [];
        highRasterX = [];
        highRasterY = [];

        uniformTrialFR = [];
        lowTrialFR = [];
        highTrialFR = [];

        uniformTrialSpikeCounts = [];
        lowTrialSpikeCounts = [];
        highTrialSpikeCounts = [];

        validUniformTrials = 0;
        validLowTrials = 0;
        validHighTrials = 0;

        for i = 1:length(uniformTrials)

            tr = uniformTrials(i);

            if tr > length(choiceTime) || isnan(choiceTime(tr))
                continue
            end

            validUniformTrials = validUniformTrials + 1;

            windowStart = choiceTime(tr) - windowBeforeSec * SampleRes;
            windowEnd   = choiceTime(tr) + windowAfterSec  * SampleRes;

            spikesInWindow = unitSpikeTimes(unitSpikeTimes >= windowStart & unitSpikeTimes <= windowEnd);
            relSpikesSec = (spikesInWindow - choiceTime(tr)) ./ SampleRes;

            uniformRasterX = [uniformRasterX; relSpikesSec(:)];
            uniformRasterY = [uniformRasterY; validUniformTrials .* ones(length(relSpikesSec), 1)];

            counts = histcounts(relSpikesSec, timeEdges);

            uniformTrialFR(validUniformTrials, :) = counts ./ binSizeSec;
            uniformTrialSpikeCounts(validUniformTrials, 1) = sum(counts);

        end

        for i = 1:length(lowTrials)

            tr = lowTrials(i);

            if tr > length(choiceTime) || isnan(choiceTime(tr))
                continue
            end

            validLowTrials = validLowTrials + 1;

            windowStart = choiceTime(tr) - windowBeforeSec * SampleRes;
            windowEnd   = choiceTime(tr) + windowAfterSec  * SampleRes;

            spikesInWindow = unitSpikeTimes(unitSpikeTimes >= windowStart & unitSpikeTimes <= windowEnd);
            relSpikesSec = (spikesInWindow - choiceTime(tr)) ./ SampleRes;

            lowRasterX = [lowRasterX; relSpikesSec(:)];
            lowRasterY = [lowRasterY; validLowTrials .* ones(length(relSpikesSec), 1)];

            counts = histcounts(relSpikesSec, timeEdges);

            lowTrialFR(validLowTrials, :) = counts ./ binSizeSec;
            lowTrialSpikeCounts(validLowTrials, 1) = sum(counts);

        end

        for i = 1:length(highTrials)

            tr = highTrials(i);

            if tr > length(choiceTime) || isnan(choiceTime(tr))
                continue
            end

            validHighTrials = validHighTrials + 1;

            windowStart = choiceTime(tr) - windowBeforeSec * SampleRes;
            windowEnd   = choiceTime(tr) + windowAfterSec  * SampleRes;

            spikesInWindow = unitSpikeTimes(unitSpikeTimes >= windowStart & unitSpikeTimes <= windowEnd);
            relSpikesSec = (spikesInWindow - choiceTime(tr)) ./ SampleRes;

            highRasterX = [highRasterX; relSpikesSec(:)];
            highRasterY = [highRasterY; validHighTrials .* ones(length(relSpikesSec), 1)];

            counts = histcounts(relSpikesSec, timeEdges);

            highTrialFR(validHighTrials, :) = counts ./ binSizeSec;
            highTrialSpikeCounts(validHighTrials, 1) = sum(counts);

        end

        if isempty(uniformTrialFR) || isempty(lowTrialFR) || isempty(highTrialFR)
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

        permStartTimes = -windowBeforeSec:permStrideSec:(0 - permWindowSec);
        permEndTimes = permStartTimes + permWindowSec;
        nWindows = length(permStartTimes);

        uniformWindowMat = nan(size(uniformZ, 1), nWindows);
        lowWindowMat     = nan(size(lowZ, 1), nWindows);
        highWindowMat    = nan(size(highZ, 1), nWindows);

        for ww = 1:nWindows

            thisStart = permStartTimes(ww);
            thisEnd   = permEndTimes(ww);

            thisIdx = timeCenters >= thisStart & timeCenters < thisEnd;

            if sum(thisIdx) < 1
                continue
            end

            uniformWindowMat(:, ww) = mean(uniformZ(:, thisIdx), 2, 'omitnan');
            lowWindowMat(:, ww)     = mean(lowZ(:, thisIdx), 2, 'omitnan');
            highWindowMat(:, ww)    = mean(highZ(:, thisIdx), 2, 'omitnan');

        end

        observedF = nan(1, nWindows);
        observedP = nan(1, nWindows);

        for ww = 1:nWindows

            y = [uniformWindowMat(:, ww); lowWindowMat(:, ww); highWindowMat(:, ww)];
            g = [repmat({'uniform'}, size(uniformWindowMat, 1), 1); ...
                 repmat({'low'},     size(lowWindowMat, 1), 1); ...
                 repmat({'high'},    size(highWindowMat, 1), 1)];

            validIdx = ~isnan(y);
            y = y(validIdx);
            g = g(validIdx);

            if numel(unique(g)) < 3
                continue
            end

            [p, tbl] = anova1(y, g, 'off');

            observedP(ww) = p;
            observedF(ww) = tbl{2,5};

        end

        fdrP = bh_fdr(observedP);
        fdrSig = fdrP < alphaPerm;

        allWindowMat = [uniformWindowMat; lowWindowMat; highWindowMat];

        nUniformHere = size(uniformWindowMat, 1);
        nLowHere     = size(lowWindowMat, 1);
        nHighHere    = size(highWindowMat, 1);
        nTotalHere   = nUniformHere + nLowHere + nHighHere;

        groupLabels = [ones(nUniformHere, 1); ...
                       2 .* ones(nLowHere, 1); ...
                       3 .* ones(nHighHere, 1)];

        permF = nan(nPerm, nWindows);

        for pp = 1:nPerm

            shuffledLabels = groupLabels(randperm(nTotalHere));

            for ww = 1:nWindows

                y = allWindowMat(:, ww);
                g = shuffledLabels;

                validIdx = ~isnan(y);
                y = y(validIdx);
                g = g(validIdx);

                if numel(unique(g)) < 3
                    continue
                end

                [~, tbl] = anova1(y, g, 'off');
                permF(pp, ww) = tbl{2,5};

            end
        end

        permThresholds = prctile(permF, 100 * (1 - alphaPerm), 1);

        sigWindowIdx = observedF > permThresholds;

        observedClusterMasses = [];
        observedClusterStarts = [];
        observedClusterEnds = [];
        observedClusterStartIdx = [];
        observedClusterEndIdx = [];

        ww = 1;
        while ww <= nWindows

            if sigWindowIdx(ww)

                clusterStartIdx = ww;

                while ww <= nWindows && sigWindowIdx(ww)
                    ww = ww + 1;
                end

                clusterEndIdx = ww - 1;

                clusterMass = sum(observedF(clusterStartIdx:clusterEndIdx), 'omitnan');

                observedClusterMasses = [observedClusterMasses; clusterMass];
                observedClusterStarts = [observedClusterStarts; permStartTimes(clusterStartIdx)];
                observedClusterEnds   = [observedClusterEnds; permEndTimes(clusterEndIdx)];
                observedClusterStartIdx = [observedClusterStartIdx; clusterStartIdx];
                observedClusterEndIdx   = [observedClusterEndIdx; clusterEndIdx];

            else
                ww = ww + 1;
            end
        end

        maxPermClusterMass = zeros(nPerm, 1);

        for pp = 1:nPerm

            permSigWindowIdx = permF(pp, :) > permThresholds;
            permClusterMasses = [];

            ww = 1;
            while ww <= nWindows

                if permSigWindowIdx(ww)

                    clusterStartIdx = ww;

                    while ww <= nWindows && permSigWindowIdx(ww)
                        ww = ww + 1;
                    end

                    clusterEndIdx = ww - 1;

                    clusterMass = sum(permF(pp, clusterStartIdx:clusterEndIdx), 'omitnan');
                    permClusterMasses = [permClusterMasses; clusterMass];

                else
                    ww = ww + 1;
                end
            end

            if ~isempty(permClusterMasses)
                maxPermClusterMass(pp) = max(permClusterMasses);
            end
        end

        sigSegments = [];
        sigLabels = {};

        for cc = 1:length(observedClusterMasses)

            clusterP = (sum(maxPermClusterMass >= observedClusterMasses(cc)) + 1) / (nPerm + 1);

            if clusterP < alphaPerm

                sigStart = observedClusterStarts(cc);
                sigEnd   = observedClusterEnds(cc);

                sigStart = max(sigStart, -windowBeforeSec);
                sigEnd   = min(sigEnd, 0);

                if sigStart < 0

                    idx1 = observedClusterStartIdx(cc);
                    idx2 = observedClusterEndIdx(cc);

                    clusterUniformVals = mean(uniformWindowMat(:, idx1:idx2), 2, 'omitnan');
                    clusterLowVals     = mean(lowWindowMat(:, idx1:idx2), 2, 'omitnan');
                    clusterHighVals    = mean(highWindowMat(:, idx1:idx2), 2, 'omitnan');

                    pairNames = {'U-L', 'U-H', 'L-H'};
                    pairP = nan(3,1);

                    pairP(1) = permutation_pair_p(clusterUniformVals, clusterLowVals, nPerm);
                    pairP(2) = permutation_pair_p(clusterUniformVals, clusterHighVals, nPerm);
                    pairP(3) = permutation_pair_p(clusterLowVals, clusterHighVals, nPerm);

                    pairPFDR = bh_fdr(pairP);

                    sigPairIdx = find(pairPFDR < alphaPerm);

                    if isempty(sigPairIdx)
                        labelTxt = 'ANOVA';
                    else
                        labelTxt = strjoin(pairNames(sigPairIdx), ', ');
                    end

                    sigSegments = [sigSegments; sigStart sigEnd];
                    sigLabels{end+1,1} = labelTxt;

                    tmpCluster = table();
                    tmpCluster.Patient = string(ptID);
                    tmpCluster.Area = string(areaNameClean);
                    tmpCluster.Channel = chanNum;
                    tmpCluster.Unit = unitNum;
                    tmpCluster.ClusterStartSec = sigStart;
                    tmpCluster.ClusterEndSec = sigEnd;
                    tmpCluster.ClusterP = clusterP;
                    tmpCluster.Pair_U_L_p = pairP(1);
                    tmpCluster.Pair_U_H_p = pairP(2);
                    tmpCluster.Pair_L_H_p = pairP(3);
                    tmpCluster.Pair_U_L_FDR = pairPFDR(1);
                    tmpCluster.Pair_U_H_FDR = pairPFDR(2);
                    tmpCluster.Pair_L_H_FDR = pairPFDR(3);
                    tmpCluster.SignificantPairs = string(labelTxt);

                    allClusterStats = [allClusterStats; tmpCluster];

                end
            end
        end

        tmpStats = table();

        tmpStats.Patient = repmat(string(ptID), nWindows, 1);
        tmpStats.Area = repmat(string(areaNameClean), nWindows, 1);
        tmpStats.Channel = repmat(chanNum, nWindows, 1);
        tmpStats.Unit = repmat(unitNum, nWindows, 1);
        tmpStats.WindowStartSec = permStartTimes(:);
        tmpStats.WindowEndSec = permEndTimes(:);
        tmpStats.ANOVA_F = observedF(:);
        tmpStats.ANOVA_p = observedP(:);
        tmpStats.FDR_p = fdrP(:);
        tmpStats.FDR_sig = fdrSig(:);

        allStats = [allStats; tmpStats];

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

        xlim([-windowBeforeSec windowAfterSec])
        xlabel('time from choice onset (s)')
        ylabel('trials')

        title(sprintf('%s | %s | Ch %d Unit %d | Raster | FR %.2f Hz', ...
            ptID, areaNameClean, chanNum, unitNum, meanFRHz), ...
            'Interpreter', 'none')

        legend({'uniform', 'low', 'high'}, 'Location', 'best')
        box off

        subplot(2,1,2)
        hold on

        fill([timeCenters fliplr(timeCenters)], ...
            [uniformMeanSmooth + uniformSEMSmooth fliplr(uniformMeanSmooth - uniformSEMSmooth)], ...
            uniformColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

        fill([timeCenters fliplr(timeCenters)], ...
            [lowMeanSmooth + lowSEMSmooth fliplr(lowMeanSmooth - lowSEMSmooth)], ...
            lowColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

        fill([timeCenters fliplr(timeCenters)], ...
            [highMeanSmooth + highSEMSmooth fliplr(highMeanSmooth - highSEMSmooth)], ...
            highColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

        plot(timeCenters, uniformMeanSmooth, ...
            'Color', uniformColor, 'LineWidth', 1.8)

        plot(timeCenters, lowMeanSmooth, ...
            'Color', lowColor, 'LineWidth', 1.8)

        plot(timeCenters, highMeanSmooth, ...
            'Color', highColor, 'LineWidth', 1.8)

        xline(0, '--k')
        yline(0, ':k')

        upperVals = [uniformMeanSmooth + uniformSEMSmooth, ...
                     lowMeanSmooth + lowSEMSmooth, ...
                     highMeanSmooth + highSEMSmooth];

        lowerVals = [uniformMeanSmooth - uniformSEMSmooth, ...
                     lowMeanSmooth - lowSEMSmooth, ...
                     highMeanSmooth - highSEMSmooth];

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
                'Color', [0.2 0.2 0.2], ...
                'LineWidth', 4)

            text(mean([sigStart sigEnd]), sigY + 0.08 * yRange, sigLabels{ss}, ...
                'Color', [0.1 0.1 0.1], ...
                'FontSize', 9, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom')

        end

        if ~isempty(sigSegments)
            ylim([curveLow - 0.15 * yRange, sigY + 0.55 * yRange])
        end

        xlim([-windowBeforeSec windowAfterSec])
        xlabel('time from choice onset (s)')
        ylabel('baseline z-scored firing rate')
        title('PSTH: uniform vs low vs high | mean ± sem');

        box off

        pdfName = sprintf('%s_%s_ch%d_unit%d_choice_distribution.pdf', ...
            ptID, areaNameClean, chanNum, unitNum);

        pdfPath = fullfile(OutputFolder, pdfName);

        exportgraphics(fig, pdfPath, 'ContentType', 'vector');
        close(fig);

        fprintf('saved: %s \n', pdfName);

    end
end

statsPath = fullfile(OutputFolder, 'choice_distribution_ANOVA_FDR_stats.csv');
writetable(allStats, statsPath);

clusterStatsPath = fullfile(OutputFolder, 'choice_distribution_cluster_pairwise_stats.csv');
writetable(allClusterStats, clusterStatsPath);

fprintf('\nSaved ANOVA/FDR stats to:\n%s\n', statsPath);
fprintf('\nSaved cluster/pairwise stats to:\n%s\n', clusterStatsPath);


function fdrP = bh_fdr(pVals)

    pVals = pVals(:);
    fdrP = nan(size(pVals));

    validIdx = ~isnan(pVals);
    p = pVals(validIdx);

    if isempty(p)
        return
    end

    [pSorted, sortIdx] = sort(p);
    m = length(pSorted);

    qSorted = pSorted .* m ./ (1:m)';

    for i = m-1:-1:1
        qSorted(i) = min(qSorted(i), qSorted(i+1));
    end

    qSorted(qSorted > 1) = 1;

    q = nan(size(p));
    q(sortIdx) = qSorted;

    fdrP(validIdx) = q;

end


function p = permutation_pair_p(x1, x2, nPerm)

    x1 = x1(~isnan(x1));
    x2 = x2(~isnan(x2));

    if isempty(x1) || isempty(x2)
        p = NaN;
        return
    end

    observedDiff = abs(mean(x1, 'omitnan') - mean(x2, 'omitnan'));

    allVals = [x1(:); x2(:)];
    n1 = length(x1);
    nTotal = length(allVals);

    permDiffs = nan(nPerm, 1);

    for pp = 1:nPerm

        idx = randperm(nTotal);

        permX1 = allVals(idx(1:n1));
        permX2 = allVals(idx(n1+1:end));

        permDiffs(pp) = abs(mean(permX1, 'omitnan') - mean(permX2, 'omitnan'));

    end

    p = (sum(permDiffs >= observedDiff) + 1) / (nPerm + 1);

end