clc;
clear;
close all;

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};

inputFolder = '\\155.100.91.44\d\Data\Nill\starling\spikes\spike_data\';
OutputFolder = '\\155.100.91.44\d\Code\Nill\Starling_units_analysis\starling_4_FR_flip_card_number\';

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

cardNums = 1:9;
nCards = length(cardNums);
cardColors = lines(nCards);

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

    cardTrials = cell(nCards, 1);
    for cc = 1:nCards
        cardTrials{cc} = find(myCard == cardNums(cc));
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

        rasterX = cell(nCards, 1);
        rasterY = cell(nCards, 1);
        trialFR = cell(nCards, 1);
        trialSpikeCounts = cell(nCards, 1);
        validTrials = zeros(nCards, 1);

        for cc = 1:nCards

            theseTrials = cardTrials{cc};

            for i = 1:length(theseTrials)

                tr = theseTrials(i);

                if tr > length(choiceTime) || isnan(choiceTime(tr))
                    continue
                end

                validTrials(cc) = validTrials(cc) + 1;

                windowStart = choiceTime(tr) - windowBeforeSec * SampleRes;
                windowEnd   = choiceTime(tr) + windowAfterSec  * SampleRes;

                spikesInWindow = unitSpikeTimes(unitSpikeTimes >= windowStart & unitSpikeTimes <= windowEnd);
                relSpikesSec = (spikesInWindow - choiceTime(tr)) ./ SampleRes;

                rasterX{cc} = [rasterX{cc}; relSpikesSec(:)];
                rasterY{cc} = [rasterY{cc}; validTrials(cc) .* ones(length(relSpikesSec), 1)];

                counts = histcounts(relSpikesSec, timeEdges);

                trialFR{cc}(validTrials(cc), :) = counts ./ binSizeSec;
                trialSpikeCounts{cc}(validTrials(cc), 1) = sum(counts);

            end
        end

        emptyGroup = false;
        for cc = 1:nCards
            if isempty(trialFR{cc})
                emptyGroup = true;
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

        trialZ = cell(nCards, 1);
        meanZ = cell(nCards, 1);
        semZ = cell(nCards, 1);
        meanSmooth = cell(nCards, 1);
        semSmooth = cell(nCards, 1);

        for cc = 1:nCards
            trialZ{cc} = (trialFR{cc} - baselineMean) ./ baselineStd;

            meanZ{cc} = mean(trialZ{cc}, 1, 'omitnan');
            semZ{cc} = std(trialZ{cc}, 0, 1, 'omitnan') ./ sqrt(size(trialZ{cc}, 1));

            meanSmooth{cc} = smoothdata(meanZ{cc}, 'gaussian', smoothBins);
            semSmooth{cc} = smoothdata(semZ{cc}, 'gaussian', smoothBins);
        end

        permStartTimes = -windowBeforeSec:permStrideSec:(0 - permWindowSec);
        permEndTimes = permStartTimes + permWindowSec;
        nWindows = length(permStartTimes);

        windowMat = cell(nCards, 1);

        for cc = 1:nCards
            windowMat{cc} = nan(size(trialZ{cc}, 1), nWindows);
        end

        for ww = 1:nWindows

            thisStart = permStartTimes(ww);
            thisEnd   = permEndTimes(ww);

            thisIdx = timeCenters >= thisStart & timeCenters < thisEnd;

            if sum(thisIdx) < 1
                continue
            end

            for cc = 1:nCards
                windowMat{cc}(:, ww) = mean(trialZ{cc}(:, thisIdx), 2, 'omitnan');
            end
        end

        observedF = nan(1, nWindows);
        observedP = nan(1, nWindows);

        for ww = 1:nWindows

            y = [];
            g = [];

            for cc = 1:nCards
                y = [y; windowMat{cc}(:, ww)];
                g = [g; cc .* ones(size(windowMat{cc}, 1), 1)];
            end

            validIdx = ~isnan(y);
            y = y(validIdx);
            g = g(validIdx);

            if numel(unique(g)) < nCards
                continue
            end

            [p, tbl] = anova1(y, g, 'off');

            observedP(ww) = p;
            observedF(ww) = tbl{2,5};

        end

        fdrP = bh_fdr(observedP);
        fdrSig = fdrP < alphaPerm;

        allWindowMat = vertcat(windowMat{:});

        groupLabels = [];
        for cc = 1:nCards
            groupLabels = [groupLabels; cc .* ones(size(windowMat{cc}, 1), 1)];
        end

        nTotalHere = length(groupLabels);

        permF = nan(nPerm, nWindows);

        for pp = 1:nPerm

            shuffledLabels = groupLabels(randperm(nTotalHere));

            for ww = 1:nWindows

                y = allWindowMat(:, ww);
                g = shuffledLabels;

                validIdx = ~isnan(y);
                y = y(validIdx);
                g = g(validIdx);

                if numel(unique(g)) < nCards
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

        for cl = 1:length(observedClusterMasses)

            clusterP = (sum(maxPermClusterMass >= observedClusterMasses(cl)) + 1) / (nPerm + 1);

            if clusterP < alphaPerm

                sigStart = observedClusterStarts(cl);
                sigEnd   = observedClusterEnds(cl);

                sigStart = max(sigStart, -windowBeforeSec);
                sigEnd   = min(sigEnd, 0);

                if sigStart < 0

                    idx1 = observedClusterStartIdx(cl);
                    idx2 = observedClusterEndIdx(cl);

                    clusterVals = cell(nCards, 1);

                    for cc = 1:nCards
                        clusterVals{cc} = mean(windowMat{cc}(:, idx1:idx2), 2, 'omitnan');
                    end

                    pairNames = {};
                    pairP = [];

                    for a = 1:nCards-1
                        for b = a+1:nCards
                            pairNames{end+1,1} = sprintf('%d-%d', cardNums(a), cardNums(b));
                            pairP(end+1,1) = permutation_pair_p(clusterVals{a}, clusterVals{b}, nPerm);
                        end
                    end

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
                    tmpCluster.SignificantPairs = string(labelTxt);

                    for ppair = 1:length(pairNames)
                        cleanPairName = regexprep(pairNames{ppair}, '-', '_');
                        tmpCluster.(sprintf('Pair_%s_p', cleanPairName)) = pairP(ppair);
                        tmpCluster.(sprintf('Pair_%s_FDR', cleanPairName)) = pairPFDR(ppair);
                    end

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

        trialOffset = 0;

        for cc = 1:nCards

            scatter(rasterX{cc}, rasterY{cc} + trialOffset, 6, ...
                'filled', ...
                'MarkerFaceColor', cardColors(cc,:), ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', 0.4)

            if cc < nCards
                yline(trialOffset + validTrials(cc) + 0.5, '--k')
            end

            trialOffset = trialOffset + validTrials(cc);

        end

        xline(0, '--k')

        xlim([-windowBeforeSec windowAfterSec])
        xlabel('time from choice onset (s)')
        ylabel('trials')

        title(sprintf('%s | %s | Ch %d Unit %d | Raster by myCard | FR %.2f Hz', ...
            ptID, areaNameClean, chanNum, unitNum, meanFRHz), ...
            'Interpreter', 'none')

        box off

        subplot(2,1,2)
        hold on

        for cc = 1:nCards

            fill([timeCenters fliplr(timeCenters)], ...
                [meanSmooth{cc} + semSmooth{cc} fliplr(meanSmooth{cc} - semSmooth{cc})], ...
                cardColors(cc,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none')

            plot(timeCenters, meanSmooth{cc}, ...
                'Color', cardColors(cc,:), 'LineWidth', 1.5)

        end

        xline(0, '--k')
        yline(0, ':k')

        upperVals = [];
        lowerVals = [];

        for cc = 1:nCards
            upperVals = [upperVals, meanSmooth{cc} + semSmooth{cc}];
            lowerVals = [lowerVals, meanSmooth{cc} - semSmooth{cc}];
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

        for ss = 1:size(sigSegments, 1)

            sigStart = sigSegments(ss, 1);
            sigEnd   = sigSegments(ss, 2);

            plot([sigStart sigEnd], [sigY sigY], ...
                '-', ...
                'Color', [0.2 0.2 0.2], ...
                'LineWidth', 4)

            text(mean([sigStart sigEnd]), sigY + 0.08 * yRange, sigLabels{ss}, ...
                'Color', [0.1 0.1 0.1], ...
                'FontSize', 8, ...
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
        title('PSTH: myCard 1-9 | mean ± sem');

        box off

        pdfName = sprintf('%s_%s_ch%d_unit%d_choice_myCard.pdf', ...
            ptID, areaNameClean, chanNum, unitNum);

        pdfPath = fullfile(OutputFolder, pdfName);

        exportgraphics(fig, pdfPath, 'ContentType', 'vector');
        close(fig);

        fprintf('saved: %s \n', pdfName);

    end
end

statsPath = fullfile(OutputFolder, 'choice_myCard_ANOVA_FDR_stats.csv');
writetable(allStats, statsPath);

clusterStatsPath = fullfile(OutputFolder, 'choice_myCard_cluster_pairwise_stats.csv');
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