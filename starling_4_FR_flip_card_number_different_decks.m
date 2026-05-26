clc;
clear;
close all;

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};

inputFolder = '\\155.100.91.44\d\Data\Nill\starling\spikes\spike_data\';
OutputFolder = '\\155.100.91.44\d\Code\Nill\Starling_units_analysis\starling_4_FR_flip_card_number_different_decks\';

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

testWindowStartSec = -windowBeforeSec;
testWindowEndSec   = 0;

permStartTimes = testWindowStartSec:permStrideSec:(testWindowEndSec - permWindowSec);
permEndTimes = permStartTimes + permWindowSec;
nWindows = length(permStartTimes);

distNames = {'uniform', 'low', 'high'};
distPlotNames = {'Uniform', 'Low', 'High'};

distBaseColors = [
    0.50 0.50 0.50
    0.85 0.35 0.00
    0.25 0.65 0.25
];

groupCardNums = cell(3, 3);

groupCardNums{1,1} = [1 2 3];
groupCardNums{1,2} = [4 5 6];
groupCardNums{1,3} = [7 8 9];

groupCardNums{2,1} = [1 2];
groupCardNums{2,2} = [3 4 5];
groupCardNums{2,3} = [6 7 8 9];

groupCardNums{3,1} = [1 2 3 4];
groupCardNums{3,2} = [5 6 7];
groupCardNums{3,3} = [8 9];

nDists = length(distNames);
nGroups = 3;

groupLabels = cell(nDists, nGroups);
for dd = 1:nDists
    for gg = 1:nGroups
        groupLabels{dd, gg} = sprintf('%d-%d', ...
            groupCardNums{dd, gg}(1), groupCardNums{dd, gg}(end));
    end
end

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

    myCard = double(bhvData.myCard);
    distribution = bhvData.distribution;

    myCard = myCard(:);

    distMask = false(length(myCard), nDists);

    for dd = 1:nDists
        distMask(:, dd) = get_distribution_mask(distribution, distNames{dd}, dd, length(myCard));
    end

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

        rasterX = cell(nDists, nGroups);
        rasterY = cell(nDists, nGroups);
        trialFR = cell(nDists, nGroups);
        trialSpikeCounts = cell(nDists, nGroups);
        validTrials = zeros(nDists, nGroups);

        for dd = 1:nDists

            for gg = 1:nGroups

                theseTrials = find(distMask(:, dd) & ismember(myCard, groupCardNums{dd, gg}));

                for i = 1:length(theseTrials)

                    tr = theseTrials(i);

                    if tr > length(choiceTime) || isnan(choiceTime(tr))
                        continue
                    end

                    validTrials(dd, gg) = validTrials(dd, gg) + 1;

                    windowStart = choiceTime(tr) - windowBeforeSec * SampleRes;
                    windowEnd   = choiceTime(tr) + windowAfterSec  * SampleRes;

                    spikesInWindow = unitSpikeTimes(unitSpikeTimes >= windowStart & unitSpikeTimes <= windowEnd);
                    relSpikesSec = (spikesInWindow - choiceTime(tr)) ./ SampleRes;

                    rasterX{dd, gg} = [rasterX{dd, gg}; relSpikesSec(:)];
                    rasterY{dd, gg} = [rasterY{dd, gg}; validTrials(dd, gg) .* ones(length(relSpikesSec), 1)];

                    counts = histcounts(relSpikesSec, timeEdges);

                    trialFR{dd, gg}(validTrials(dd, gg), :) = counts ./ binSizeSec;
                    trialSpikeCounts{dd, gg}(validTrials(dd, gg), 1) = sum(counts);

                end
            end
        end

        emptyGroup = false;

        for dd = 1:nDists
            for gg = 1:nGroups
                if isempty(trialFR{dd, gg})
                    emptyGroup = true;
                end
            end
        end

        if emptyGroup
            continue
        end

        allTrialSpikeCounts = vertcat(trialSpikeCounts{:});

        nValidTrialsTotal = length(allTrialSpikeCounts);
        totalSpikesInWindow = sum(allTrialSpikeCounts);

        meanFRHz = totalSpikesInWindow / (nValidTrialsTotal * analysisWindowSec);
        fracTrialsWithSpikes = sum(allTrialSpikeCounts > 0) / nValidTrialsTotal;

        if meanFRHz < minMeanFRHz || fracTrialsWithSpikes < minFracTrialsWithSpikes
            continue
        end

        allTrialFR = vertcat(trialFR{:});

        baselineVals = allTrialFR(:, baselineIdx);
        baselineMean = mean(baselineVals(:), 'omitnan');
        baselineStd  = std(baselineVals(:), 'omitnan');

        if baselineStd == 0 || isnan(baselineStd)
            baselineStd = 1;
        end

        trialZ = cell(nDists, nGroups);
        meanZ = cell(nDists, nGroups);
        semZ = cell(nDists, nGroups);
        meanSmooth = cell(nDists, nGroups);
        semSmooth = cell(nDists, nGroups);

        for dd = 1:nDists

            for gg = 1:nGroups

                trialZ{dd, gg} = (trialFR{dd, gg} - baselineMean) ./ baselineStd;

                meanZ{dd, gg} = mean(trialZ{dd, gg}, 1, 'omitnan');
                semZ{dd, gg} = std(trialZ{dd, gg}, 0, 1, 'omitnan') ./ sqrt(size(trialZ{dd, gg}, 1));

                meanSmooth{dd, gg} = smoothdata(meanZ{dd, gg}, 'gaussian', smoothBins);
                semSmooth{dd, gg} = smoothdata(semZ{dd, gg}, 'gaussian', smoothBins);

            end
        end

        sigSegments = cell(nDists, 1);
        sigLabels = cell(nDists, 1);

        for dd = 1:nDists

            thisTrialZ = trialZ(dd, :);
            thisGroupLabels = groupLabels(dd, :);

            [tmpStats, tmpClusterStats, sigSegments{dd}, sigLabels{dd}] = ...
                run_cluster_anova_3groups( ...
                thisTrialZ, ...
                thisGroupLabels, ...
                timeCenters, ...
                permStartTimes, ...
                permEndTimes, ...
                nPerm, ...
                alphaPerm);

            nRows = height(tmpStats);

            tmpStats.Patient = repmat(string(ptID), nRows, 1);
            tmpStats.Area = repmat(string(areaNameClean), nRows, 1);
            tmpStats.Channel = repmat(chanNum, nRows, 1);
            tmpStats.Unit = repmat(unitNum, nRows, 1);
            tmpStats.Distribution = repmat(string(distNames{dd}), nRows, 1);

            tmpStats = movevars(tmpStats, ...
                {'Patient', 'Area', 'Channel', 'Unit', 'Distribution'}, ...
                'Before', 1);

            allStats = [allStats; tmpStats];

            if ~isempty(tmpClusterStats)

                nClusterRows = height(tmpClusterStats);

                tmpClusterStats.Patient = repmat(string(ptID), nClusterRows, 1);
                tmpClusterStats.Area = repmat(string(areaNameClean), nClusterRows, 1);
                tmpClusterStats.Channel = repmat(chanNum, nClusterRows, 1);
                tmpClusterStats.Unit = repmat(unitNum, nClusterRows, 1);
                tmpClusterStats.Distribution = repmat(string(distNames{dd}), nClusterRows, 1);

                tmpClusterStats = movevars(tmpClusterStats, ...
                    {'Patient', 'Area', 'Channel', 'Unit', 'Distribution'}, ...
                    'Before', 1);

                allClusterStats = [allClusterStats; tmpClusterStats];

            end
        end

        fig = figure('Visible', 'off', 'Color', 'w');
        set(fig, 'Position', [100 100 1550 850]);

        tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

        for dd = 1:nDists

            groupColors = make_group_colors(distBaseColors(dd, :), nGroups);

            axRaster = nexttile(dd);
            hold(axRaster, 'on');

            trialOffset = 0;
            hLeg = gobjects(nGroups, 1);

            for gg = 1:nGroups

                hLeg(gg) = scatter(axRaster, nan, nan, 30, ...
                    groupColors(gg, :), ...
                    'filled');

                scatter(axRaster, rasterX{dd, gg}, rasterY{dd, gg} + trialOffset, 6, ...
                    'filled', ...
                    'MarkerFaceColor', groupColors(gg, :), ...
                    'MarkerEdgeColor', 'none', ...
                    'MarkerFaceAlpha', 0.4);

                if gg < nGroups
                    yline(axRaster, trialOffset + validTrials(dd, gg) + 0.5, '--k');
                end

                trialOffset = trialOffset + validTrials(dd, gg);

            end

            xline(axRaster, 0, '--k');

            xlim(axRaster, [-windowBeforeSec windowAfterSec]);
            xlabel(axRaster, 'time from choice onset (s)');
            ylabel(axRaster, 'trials');

            title(axRaster, sprintf('%s spikes', distPlotNames{dd}), ...
                'Interpreter', 'none');

            legend(axRaster, hLeg, groupLabels(dd, :), ...
                'Location', 'best', ...
                'Box', 'off');

            box(axRaster, 'off');

            axPSTH = nexttile(dd + 3);
            hold(axPSTH, 'on');

            for gg = 1:nGroups

                fill(axPSTH, ...
                    [timeCenters fliplr(timeCenters)], ...
                    [meanSmooth{dd, gg} + semSmooth{dd, gg}, ...
                    fliplr(meanSmooth{dd, gg} - semSmooth{dd, gg})], ...
                    groupColors(gg, :), ...
                    'FaceAlpha', 0.12, ...
                    'EdgeColor', 'none', ...
                    'HandleVisibility', 'off');

                plot(axPSTH, timeCenters, meanSmooth{dd, gg}, ...
                    'Color', groupColors(gg, :), ...
                    'LineWidth', 1.5, ...
                    'HandleVisibility', 'off');

            end

            xline(axPSTH, 0, '--k');
            yline(axPSTH, 0, ':k');

            upperVals = [];
            lowerVals = [];

            for gg = 1:nGroups
                upperVals = [upperVals, meanSmooth{dd, gg} + semSmooth{dd, gg}];
                lowerVals = [lowerVals, meanSmooth{dd, gg} - semSmooth{dd, gg}];
            end

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
            maxSigY = sigY;

            for ss = 1:size(sigSegments{dd}, 1)

                sigStart = sigSegments{dd}(ss, 1);
                sigEnd   = sigSegments{dd}(ss, 2);

                thisSigY = sigY + (ss - 1) * 0.18 * yRange;
                maxSigY = max(maxSigY, thisSigY);

                plot(axPSTH, [sigStart sigEnd], [thisSigY thisSigY], ...
                    '-', ...
                    'Color', [0.2 0.2 0.2], ...
                    'LineWidth', 4, ...
                    'HandleVisibility', 'off');

                text(axPSTH, mean([sigStart sigEnd]), thisSigY + 0.06 * yRange, ...
                    sigLabels{dd}{ss}, ...
                    'Color', [0.1 0.1 0.1], ...
                    'FontSize', 8, ...
                    'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'Interpreter', 'none');

            end

            if ~isempty(sigSegments{dd})
                ylim(axPSTH, [curveLow - 0.15 * yRange, maxSigY + 0.45 * yRange]);
            end

            xlim(axPSTH, [-windowBeforeSec windowAfterSec]);
            xlabel(axPSTH, 'time from choice onset (s)');
            ylabel(axPSTH, 'baseline z-scored firing rate');

            title(axPSTH, sprintf('%s PSTH', distPlotNames{dd}), ...
                'Interpreter', 'none');

            legend(axPSTH, 'off');
            box(axPSTH, 'off');

        end

        sgtitle(sprintf('%s | %s | Ch %d Unit %d | grouped myCard by distribution | FR %.2f Hz', ...
            ptID, areaNameClean, chanNum, unitNum, meanFRHz), ...
            'Interpreter', 'none');

        pdfName = sprintf('%s_%s_ch%d_unit%d_choice_distribution_cardGroups.pdf', ...
            ptID, areaNameClean, chanNum, unitNum);

        pdfPath = fullfile(OutputFolder, pdfName);

        exportgraphics(fig, pdfPath, 'ContentType', 'vector');
        close(fig);

        fprintf('saved: %s \n', pdfName);

    end
end

statsPath = fullfile(OutputFolder, 'choice_distribution_cardGroups_ANOVA_FDR_stats.csv');
writetable(allStats, statsPath);

clusterStatsPath = fullfile(OutputFolder, 'choice_distribution_cardGroups_cluster_pairwise_stats.csv');
writetable(allClusterStats, clusterStatsPath);

fprintf('\nSaved ANOVA/FDR stats to:\n%s\n', statsPath);
fprintf('\nSaved cluster/pairwise stats to:\n%s\n', clusterStatsPath);


function mask = get_distribution_mask(distribution, distName, distIndex, nTrials)

    mask = false(nTrials, 1);

    if isnumeric(distribution) || islogical(distribution)

        distribution = double(distribution(:));

        nUse = min(length(distribution), nTrials);
        mask(1:nUse) = distribution(1:nUse) == distIndex;

        return
    end

    if iscategorical(distribution)
        distStr = string(distribution);
    elseif isstring(distribution)
        distStr = distribution;
    elseif iscell(distribution)
        distStr = string(distribution);
    elseif ischar(distribution)
        distStr = string(cellstr(distribution));
    else
        error('distribution field type is not recognized.');
    end

    distStr = lower(strtrim(distStr(:)));

    nUse = min(length(distStr), nTrials);

    switch lower(distName)

        case 'uniform'
            thisMask = distStr(1:nUse) == "uniform" | ...
                       distStr(1:nUse) == "uni" | ...
                       distStr(1:nUse) == "unif" | ...
                       contains(distStr(1:nUse), "uniform");

        case 'low'
            thisMask = distStr(1:nUse) == "low" | ...
                       contains(distStr(1:nUse), "low");

        case 'high'
            thisMask = distStr(1:nUse) == "high" | ...
                       contains(distStr(1:nUse), "high");

        otherwise
            error('Unknown distribution name.');
    end

    mask(1:nUse) = thisMask;

end


function [windowStats, clusterStats, sigSegments, sigLabels] = run_cluster_anova_3groups( ...
    trialZ, groupLabels, timeCenters, permStartTimes, permEndTimes, nPerm, alphaPerm)

    nGroups = length(trialZ);
    nWindows = length(permStartTimes);

    windowMat = cell(nGroups, 1);

    for gg = 1:nGroups
        windowMat{gg} = nan(size(trialZ{gg}, 1), nWindows);
    end

    for ww = 1:nWindows

        thisStart = permStartTimes(ww);
        thisEnd   = permEndTimes(ww);

        thisIdx = timeCenters >= thisStart & timeCenters < thisEnd;

        if sum(thisIdx) < 1
            continue
        end

        for gg = 1:nGroups
            windowMat{gg}(:, ww) = mean(trialZ{gg}(:, thisIdx), 2, 'omitnan');
        end
    end

    observedF = nan(1, nWindows);
    observedP = nan(1, nWindows);

    for ww = 1:nWindows

        y = [];
        g = [];

        for gg = 1:nGroups
            y = [y; windowMat{gg}(:, ww)];
            g = [g; gg .* ones(size(windowMat{gg}, 1), 1)];
        end

        validIdx = ~isnan(y);
        y = y(validIdx);
        g = g(validIdx);

        if numel(unique(g)) < nGroups
            continue
        end

        [p, tbl] = anova1(y, g, 'off');

        observedP(ww) = p;
        observedF(ww) = tbl{2,5};

    end

    fdrP = bh_fdr(observedP);
    fdrSig = fdrP < alphaPerm;

    allWindowMat = vertcat(windowMat{:});

    originalGroupLabels = [];

    for gg = 1:nGroups
        originalGroupLabels = [originalGroupLabels; gg .* ones(size(windowMat{gg}, 1), 1)];
    end

    nTotalHere = length(originalGroupLabels);

    permF = nan(nPerm, nWindows);

    for pp = 1:nPerm

        shuffledLabels = originalGroupLabels(randperm(nTotalHere));

        for ww = 1:nWindows

            y = allWindowMat(:, ww);
            g = shuffledLabels;

            validIdx = ~isnan(y);
            y = y(validIdx);
            g = g(validIdx);

            if numel(unique(g)) < nGroups
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
    clusterStats = table();

    for cl = 1:length(observedClusterMasses)

        clusterP = (sum(maxPermClusterMass >= observedClusterMasses(cl)) + 1) / (nPerm + 1);

        if clusterP < alphaPerm

            sigStart = observedClusterStarts(cl);
            sigEnd   = observedClusterEnds(cl);

            idx1 = observedClusterStartIdx(cl);
            idx2 = observedClusterEndIdx(cl);

            clusterVals = cell(nGroups, 1);

            for gg = 1:nGroups
                clusterVals{gg} = mean(windowMat{gg}(:, idx1:idx2), 2, 'omitnan');
            end

            pairNames = {};
            pairP = [];

            for a = 1:nGroups-1
                for b = a+1:nGroups
                    pairNames{end+1,1} = sprintf('%s vs %s', groupLabels{a}, groupLabels{b});
                    pairP(end+1,1) = permutation_pair_p(clusterVals{a}, clusterVals{b}, nPerm);
                end
            end

            pairPFDR = bh_fdr(pairP);

            sigPairIdx = find(pairPFDR < alphaPerm);

            if isempty(sigPairIdx)
                labelTxt = 'ANOVA';
            else
                labelTxt = strjoin(pairNames(sigPairIdx), '; ');
            end

            sigSegments = [sigSegments; sigStart sigEnd];
            sigLabels{end+1,1} = labelTxt;

            tmpCluster = table();

            tmpCluster.ClusterStartSec = sigStart;
            tmpCluster.ClusterEndSec = sigEnd;
            tmpCluster.ClusterP = clusterP;
            tmpCluster.SignificantPairs = string(labelTxt);

            tmpCluster.Pair_1_vs_2 = string(pairNames{1});
            tmpCluster.Pair_1_vs_2_p = pairP(1);
            tmpCluster.Pair_1_vs_2_FDR = pairPFDR(1);

            tmpCluster.Pair_1_vs_3 = string(pairNames{2});
            tmpCluster.Pair_1_vs_3_p = pairP(2);
            tmpCluster.Pair_1_vs_3_FDR = pairPFDR(2);

            tmpCluster.Pair_2_vs_3 = string(pairNames{3});
            tmpCluster.Pair_2_vs_3_p = pairP(3);
            tmpCluster.Pair_2_vs_3_FDR = pairPFDR(3);

            clusterStats = [clusterStats; tmpCluster];

        end
    end

    windowStats = table();

    windowStats.WindowStartSec = permStartTimes(:);
    windowStats.WindowEndSec = permEndTimes(:);
    windowStats.ANOVA_F = observedF(:);
    windowStats.ANOVA_p = observedP(:);
    windowStats.FDR_p = fdrP(:);
    windowStats.FDR_sig = fdrSig(:);

end


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


function groupColors = make_group_colors(baseColor, nGroups)

    baseColor = baseColor(:)';

    mixVals = linspace(0.65, 0.05, nGroups);
    groupColors = nan(nGroups, 3);

    for gg = 1:nGroups

        mixWhite = mixVals(gg);
        groupColors(gg, :) = baseColor .* (1 - mixWhite) + [1 1 1] .* mixWhite;

    end

    groupColors(groupColors > 1) = 1;
    groupColors(groupColors < 0) = 0;

end