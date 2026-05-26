% pay attention: don't freak out about the windowBefore that is 3 sec. the
% start of the time is from choice and feedback time. 

clc;
clear;
close all;

rng('shuffle');

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};

inputFolder = '\\155.100.91.44\d\Data\Nill\starling\spikes\spike_data\';
OutputFolder = '\\155.100.91.44\d\Code\Nill\Starling_units_analysis\starling_7_FR_flip_card_number_fix_mix\';

if ~exist(OutputFolder, 'dir')
    mkdir(OutputFolder);
end

windowBeforeSec = 3;
windowAfterSec  = 0.01;
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

cardValues = 1:9;
nCards = length(cardValues);

cardColors = lines(nCards);

% Same color for fix and mix, different alphas
% Fix = kam rang tar / lighter
% Mix = por rang tar / stronger

fixRasterAlpha = 0.25;
mixRasterAlpha = 0.65;

fixLineAlpha = 0.45;
mixLineAlpha = 0.95;

fixSEMAlpha = 0.08;
mixSEMAlpha = 0.18;

allStats = table();
allClusterStats = table();

for pt = 1:length(microPts)

    ptID = microPts{pt};
    fprintf('\nprocessing patient %s\n', ptID);

    dataPath = fullfile(inputFolder, sprintf('%s_spikeData.mat', ptID));

    if ~exist(dataPath, 'file')
        fprintf('missing file: %s\n', dataPath);
        continue
    end

    load(dataPath);

    ChanUnitTimestamp = spikeData.ChanUnitTimestamp;
    eventTime = spikeData.eventTimes;
    bhvData = eventTime.bhvData;
    inclChans = spikeData.inclChans;
    microLabels = spikeData.microLabels;

    alignField = 'choiceAndFeedbackTime';
    alignTime = double(spikeData.eventTimes.(alignField)(:));

    SampleRes = double(spikeData.SampleRes);

    myCard = to_numeric_vector(bhvData.myCard);
    block = to_numeric_vector(bhvData.block);

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

        if iscell(areaName)
            areaName = areaName{1};
        end

        areaNameClean = regexprep(char(string(areaName)), '[^\w]', '_');

        unitSpikeTimes = ChanUnitTimestamp( ...
            double(ChanUnitTimestamp(:,1)) == chanNum & ...
            double(ChanUnitTimestamp(:,2)) == unitNum, 3);

        unitSpikeTimes = double(unitSpikeTimes);

        if isempty(unitSpikeTimes)
            continue
        end

        fixRasterX = cell(1, nCards);
        fixRasterY = cell(1, nCards);
        mixRasterX = cell(1, nCards);
        mixRasterY = cell(1, nCards);

        fixTrialFR = cell(1, nCards);
        mixTrialFR = cell(1, nCards);

        fixTrialSpikeCounts = cell(1, nCards);
        mixTrialSpikeCounts = cell(1, nCards);

        validFixTrials = zeros(1, nCards);
        validMixTrials = zeros(1, nCards);

        for cardIdx = 1:nCards

            thisCard = cardValues(cardIdx);

            fixTrials = find(myCard == thisCard & ismember(block, [1 2 3]));
            mixTrials = find(myCard == thisCard & block == 4);

            [fixRasterX{cardIdx}, fixRasterY{cardIdx}, fixTrialFR{cardIdx}, fixTrialSpikeCounts{cardIdx}, validFixTrials(cardIdx)] = ...
                collect_trials_for_condition( ...
                    fixTrials, ...
                    alignTime, ...
                    unitSpikeTimes, ...
                    SampleRes, ...
                    timeEdges, ...
                    windowBeforeSec, ...
                    windowAfterSec, ...
                    binSizeSec);

            [mixRasterX{cardIdx}, mixRasterY{cardIdx}, mixTrialFR{cardIdx}, mixTrialSpikeCounts{cardIdx}, validMixTrials(cardIdx)] = ...
                collect_trials_for_condition( ...
                    mixTrials, ...
                    alignTime, ...
                    unitSpikeTimes, ...
                    SampleRes, ...
                    timeEdges, ...
                    windowBeforeSec, ...
                    windowAfterSec, ...
                    binSizeSec);

        end

        allTrialFR = [];
        allTrialSpikeCounts = [];

        for cardIdx = 1:nCards

            if ~isempty(fixTrialFR{cardIdx})
                allTrialFR = [allTrialFR; fixTrialFR{cardIdx}];
                allTrialSpikeCounts = [allTrialSpikeCounts; fixTrialSpikeCounts{cardIdx}];
            end

            if ~isempty(mixTrialFR{cardIdx})
                allTrialFR = [allTrialFR; mixTrialFR{cardIdx}];
                allTrialSpikeCounts = [allTrialSpikeCounts; mixTrialSpikeCounts{cardIdx}];
            end

        end

        if isempty(allTrialFR) || isempty(allTrialSpikeCounts)
            continue
        end

        nValidTrialsTotal = length(allTrialSpikeCounts);
        totalSpikesInWindow = sum(allTrialSpikeCounts);

        meanFRHz = totalSpikesInWindow / (nValidTrialsTotal * analysisWindowSec);
        fracTrialsWithSpikes = sum(allTrialSpikeCounts > 0) / nValidTrialsTotal;

        if meanFRHz < minMeanFRHz || fracTrialsWithSpikes < minFracTrialsWithSpikes
            continue
        end

        baselineVals = allTrialFR(:, baselineIdx);
        baselineMean = mean(baselineVals(:), 'omitnan');
        baselineStd  = std(baselineVals(:), 'omitnan');

        if baselineStd == 0 || isnan(baselineStd)
            baselineStd = 1;
        end

        fixZ = cell(1, nCards);
        mixZ = cell(1, nCards);

        for cardIdx = 1:nCards

            if ~isempty(fixTrialFR{cardIdx})
                fixZ{cardIdx} = (fixTrialFR{cardIdx} - baselineMean) ./ baselineStd;
            else
                fixZ{cardIdx} = [];
            end

            if ~isempty(mixTrialFR{cardIdx})
                mixZ{cardIdx} = (mixTrialFR{cardIdx} - baselineMean) ./ baselineStd;
            else
                mixZ{cardIdx} = [];
            end

        end

        testOut = cell(1, nCards);

        for cardIdx = 1:nCards

            thisCard = cardValues(cardIdx);

            testOut{cardIdx} = run_fix_mix_cluster_test( ...
                fixZ{cardIdx}, ...
                mixZ{cardIdx}, ...
                timeCenters, ...
                windowBeforeSec, ...
                permWindowSec, ...
                permStrideSec, ...
                nPerm, ...
                alphaPerm);

            nWindows = length(testOut{cardIdx}.permStartTimes);

            tmpStats = table();

            tmpStats.Patient = repmat(string(ptID), nWindows, 1);
            tmpStats.Area = repmat(string(areaNameClean), nWindows, 1);
            tmpStats.Channel = repmat(chanNum, nWindows, 1);
            tmpStats.Unit = repmat(unitNum, nWindows, 1);
            tmpStats.MyCard = repmat(thisCard, nWindows, 1);
            tmpStats.FixTrials_n = repmat(validFixTrials(cardIdx), nWindows, 1);
            tmpStats.MixTrials_n = repmat(validMixTrials(cardIdx), nWindows, 1);
            tmpStats.WindowStartSec = testOut{cardIdx}.permStartTimes(:);
            tmpStats.WindowEndSec = testOut{cardIdx}.permEndTimes(:);
            tmpStats.FixMean = testOut{cardIdx}.fixWindowMean(:);
            tmpStats.MixMean = testOut{cardIdx}.mixWindowMean(:);
            tmpStats.MixMinusFix = testOut{cardIdx}.observedDiff(:);
            tmpStats.AbsT = testOut{cardIdx}.observedT(:);
            tmpStats.Perm_p = testOut{cardIdx}.permP(:);
            tmpStats.FDR_p = testOut{cardIdx}.fdrP(:);
            tmpStats.FDR_sig = testOut{cardIdx}.fdrSig(:);
            tmpStats.ClusterCorrected_sig = testOut{cardIdx}.clusterSigWindowIdx(:);

            allStats = [allStats; tmpStats];

            if ~isempty(testOut{cardIdx}.clusterP)

                for clIdx = 1:length(testOut{cardIdx}.clusterP)

                    tmpCluster = table();

                    tmpCluster.Patient = string(ptID);
                    tmpCluster.Area = string(areaNameClean);
                    tmpCluster.Channel = chanNum;
                    tmpCluster.Unit = unitNum;
                    tmpCluster.MyCard = thisCard;
                    tmpCluster.FixTrials_n = validFixTrials(cardIdx);
                    tmpCluster.MixTrials_n = validMixTrials(cardIdx);
                    tmpCluster.ClusterStartSec = testOut{cardIdx}.sigSegments(clIdx, 1);
                    tmpCluster.ClusterEndSec = testOut{cardIdx}.sigSegments(clIdx, 2);
                    tmpCluster.ClusterMass = testOut{cardIdx}.clusterMass(clIdx);
                    tmpCluster.ClusterP = testOut{cardIdx}.clusterP(clIdx);
                    tmpCluster.Direction = string(testOut{cardIdx}.sigLabels{clIdx});

                    allClusterStats = [allClusterStats; tmpCluster];

                end

            end

        end

        fig = figure('Visible', 'off', 'Color', 'w');
        set(fig, 'Position', [100 100 3200 900]);

        tiledlayout(2, nCards, 'TileSpacing', 'compact', 'Padding', 'compact');

        for cardIdx = 1:nCards

            thisCard = cardValues(cardIdx);
            thisColor = cardColors(cardIdx, :);

            nexttile(cardIdx)
            hold on

            hFixRaster = scatter(fixRasterX{cardIdx}, fixRasterY{cardIdx}, 6, ...
                'filled', ...
                'MarkerFaceColor', thisColor, ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', fixRasterAlpha, ...
                'MarkerEdgeAlpha', fixRasterAlpha);

            hMixRaster = scatter(mixRasterX{cardIdx}, mixRasterY{cardIdx} + validFixTrials(cardIdx), 6, ...
                'filled', ...
                'MarkerFaceColor', thisColor, ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', mixRasterAlpha, ...
                'MarkerEdgeAlpha', mixRasterAlpha);

            xline(0, '--k', 'HandleVisibility', 'off')

            if validFixTrials(cardIdx) > 0 && validMixTrials(cardIdx) > 0
                yline(validFixTrials(cardIdx) + 0.5, '--k', 'HandleVisibility', 'off')
            end

            totalRasterTrials = validFixTrials(cardIdx) + validMixTrials(cardIdx);

            if totalRasterTrials > 0
                ylim([0.5 totalRasterTrials + 0.5])
            end

            xlim([-windowBeforeSec windowAfterSec])

            title(sprintf('myCard %d raster | fix n=%d | mix n=%d', ...
                thisCard, validFixTrials(cardIdx), validMixTrials(cardIdx)), ...
                'Interpreter', 'none')

            xlabel(sprintf('time from %s (s)', alignField))

            if cardIdx == 1
                ylabel('trials')
            else
                ylabel('')
            end

            legend([hFixRaster hMixRaster], ...
                {'fix', 'mix'}, ...
                'Location', 'best')

            box off

        end

        for cardIdx = 1:nCards

            thisCard = cardValues(cardIdx);
            thisColor = cardColors(cardIdx, :);

            [fixMeanSmooth, fixSEMSmooth] = mean_sem_smooth(fixZ{cardIdx}, smoothBins);
            [mixMeanSmooth, mixSEMSmooth] = mean_sem_smooth(mixZ{cardIdx}, smoothBins);

            nexttile(cardIdx + nCards)
            hold on

            legendHandles = gobjects(0);
            legendLabels = {};

            if ~isempty(fixMeanSmooth)

                fill([timeCenters fliplr(timeCenters)], ...
                    [fixMeanSmooth + fixSEMSmooth fliplr(fixMeanSmooth - fixSEMSmooth)], ...
                    thisColor, ...
                    'FaceAlpha', fixSEMAlpha, ...
                    'EdgeColor', 'none', ...
                    'HandleVisibility', 'off')

                hFixLine = plot_alpha_line(timeCenters, fixMeanSmooth, thisColor, 1.8, fixLineAlpha);
                legendHandles(end+1) = hFixLine;
                legendLabels{end+1} = 'fix';

            end

            if ~isempty(mixMeanSmooth)

                fill([timeCenters fliplr(timeCenters)], ...
                    [mixMeanSmooth + mixSEMSmooth fliplr(mixMeanSmooth - mixSEMSmooth)], ...
                    thisColor, ...
                    'FaceAlpha', mixSEMAlpha, ...
                    'EdgeColor', 'none', ...
                    'HandleVisibility', 'off')

                hMixLine = plot_alpha_line(timeCenters, mixMeanSmooth, thisColor, 1.8, mixLineAlpha);
                legendHandles(end+1) = hMixLine;
                legendLabels{end+1} = 'mix';

            end

            xline(0, '--k', 'HandleVisibility', 'off')
            yline(0, ':k', 'HandleVisibility', 'off')

            curveVals = [];

            if ~isempty(fixMeanSmooth)
                curveVals = [curveVals, fixMeanSmooth + fixSEMSmooth, fixMeanSmooth - fixSEMSmooth];
            end

            if ~isempty(mixMeanSmooth)
                curveVals = [curveVals, mixMeanSmooth + mixSEMSmooth, mixMeanSmooth - mixSEMSmooth];
            end

            curveVals = curveVals(~isnan(curveVals));

            if isempty(curveVals)
                curveHigh = 1;
                curveLow = -1;
            else
                curveHigh = max(curveVals);
                curveLow = min(curveVals);
            end

            yRange = curveHigh - curveLow;

            if yRange == 0 || isnan(yRange)
                yRange = 1;
            end

            sigY = curveHigh + 0.15 * yRange;

            for ss = 1:size(testOut{cardIdx}.sigSegments, 1)

                sigStart = testOut{cardIdx}.sigSegments(ss, 1);
                sigEnd   = testOut{cardIdx}.sigSegments(ss, 2);

                plot([sigStart sigEnd], [sigY sigY], ...
                    '-', ...
                    'Color', [0.15 0.15 0.15], ...
                    'LineWidth', 4, ...
                    'HandleVisibility', 'off')

                text(mean([sigStart sigEnd]), sigY + 0.08 * yRange, testOut{cardIdx}.sigLabels{ss}, ...
                    'Color', [0.05 0.05 0.05], ...
                    'FontSize', 9, ...
                    'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom')

            end

            if ~isempty(testOut{cardIdx}.sigSegments)
                ylim([curveLow - 0.15 * yRange, sigY + 0.55 * yRange])
            else
                ylim([curveLow - 0.15 * yRange, curveHigh + 0.25 * yRange])
            end

            xlim([-windowBeforeSec windowAfterSec])

            title(sprintf('myCard %d PSTH: fix vs mix', thisCard), ...
                'Interpreter', 'none')

            xlabel(sprintf('time from %s (s)', alignField))

            if cardIdx == 1
                ylabel('baseline z-scored firing rate')
            else
                ylabel('')
            end

            if ~isempty(legendHandles)
                legend(legendHandles, legendLabels, 'Location', 'best')
            end

            box off

        end

        sgtitle(sprintf('%s | %s | Ch %d Unit %d | FR %.2f Hz | Fix blocks vs Mix block by myCard', ...
            ptID, areaNameClean, chanNum, unitNum, meanFRHz), ...
            'Interpreter', 'none', ...
            'FontWeight', 'bold')

        pdfName = sprintf('%s_%s_ch%d_unit%d_choice_fix_vs_mix_by_myCard.pdf', ...
            ptID, areaNameClean, chanNum, unitNum);

        pdfPath = fullfile(OutputFolder, pdfName);

        exportgraphics(fig, pdfPath, 'ContentType', 'vector');
        close(fig);

        fprintf('saved: %s \n', pdfName);

    end

end

statsPath = fullfile(OutputFolder, 'choice_fix_vs_mix_by_myCard_window_stats.csv');
writetable(allStats, statsPath);

clusterStatsPath = fullfile(OutputFolder, 'choice_fix_vs_mix_by_myCard_cluster_stats.csv');
writetable(allClusterStats, clusterStatsPath);

fprintf('\nSaved window stats to:\n%s\n', statsPath);
fprintf('\nSaved cluster stats to:\n%s\n', clusterStatsPath);


function x = to_numeric_vector(x)

    if iscell(x) || isstring(x) || ischar(x) || iscategorical(x)
        x = str2double(string(x));
    else
        x = double(x);
    end

    x = x(:);

end


function [rasterX, rasterY, trialFR, trialSpikeCounts, validTrials] = collect_trials_for_condition( ...
    trials, alignTime, unitSpikeTimes, SampleRes, timeEdges, windowBeforeSec, windowAfterSec, binSizeSec)

    rasterX = [];
    rasterY = [];
    trialFR = [];
    trialSpikeCounts = [];
    validTrials = 0;

    for ii = 1:length(trials)

        tr = trials(ii);

        if tr > length(alignTime) || isnan(alignTime(tr))
            continue
        end

        validTrials = validTrials + 1;

        windowStart = alignTime(tr) - windowBeforeSec * SampleRes;
        windowEnd   = alignTime(tr) + windowAfterSec  * SampleRes;

        spikesInWindow = unitSpikeTimes(unitSpikeTimes >= windowStart & unitSpikeTimes <= windowEnd);
        relSpikesSec = (spikesInWindow - alignTime(tr)) ./ SampleRes;

        rasterX = [rasterX; relSpikesSec(:)];
        rasterY = [rasterY; validTrials .* ones(length(relSpikesSec), 1)];

        counts = histcounts(relSpikesSec, timeEdges);

        trialFR(validTrials, :) = counts ./ binSizeSec;
        trialSpikeCounts(validTrials, 1) = sum(counts);

    end

end


function [meanSmooth, semSmooth] = mean_sem_smooth(zMat, smoothBins)

    if isempty(zMat)
        meanSmooth = [];
        semSmooth = [];
        return
    end

    meanVals = mean(zMat, 1, 'omitnan');

    nVals = sum(~isnan(zMat), 1);
    semVals = std(zMat, 0, 1, 'omitnan') ./ sqrt(nVals);
    semVals(nVals == 0) = NaN;

    meanSmooth = smoothdata(meanVals, 'gaussian', smoothBins);
    semSmooth  = smoothdata(semVals,  'gaussian', smoothBins);

end


function h = plot_alpha_line(x, y, colorVal, lineWidthVal, alphaVal)

    x = x(:)';
    y = y(:)';

    h = surface( ...
        [x; x], ...
        [y; y], ...
        zeros(2, numel(x)), ...
        'FaceColor', 'none', ...
        'EdgeColor', colorVal, ...
        'LineWidth', lineWidthVal, ...
        'EdgeAlpha', alphaVal);

end


function testOut = run_fix_mix_cluster_test( ...
    fixZ, mixZ, timeCenters, windowBeforeSec, permWindowSec, permStrideSec, nPerm, alphaPerm)

    permStartTimes = -windowBeforeSec:permStrideSec:(0 - permWindowSec);
    permEndTimes = permStartTimes + permWindowSec;
    nWindows = length(permStartTimes);

    testOut = struct();

    testOut.hasTest = false;
    testOut.permStartTimes = permStartTimes;
    testOut.permEndTimes = permEndTimes;

    testOut.fixWindowMean = nan(1, nWindows);
    testOut.mixWindowMean = nan(1, nWindows);
    testOut.observedDiff = nan(1, nWindows);
    testOut.observedT = nan(1, nWindows);
    testOut.permP = nan(1, nWindows);
    testOut.fdrP = nan(1, nWindows);
    testOut.fdrSig = false(1, nWindows);
    testOut.clusterSigWindowIdx = false(1, nWindows);

    testOut.sigSegments = [];
    testOut.sigLabels = {};
    testOut.clusterP = [];
    testOut.clusterMass = [];

    if isempty(fixZ) || isempty(mixZ)
        return
    end

    if size(fixZ, 1) < 2 || size(mixZ, 1) < 2
        return
    end

    fixWindowMat = nan(size(fixZ, 1), nWindows);
    mixWindowMat = nan(size(mixZ, 1), nWindows);

    for ww = 1:nWindows

        thisStart = permStartTimes(ww);
        thisEnd   = permEndTimes(ww);

        thisIdx = timeCenters >= thisStart & timeCenters < thisEnd;

        if sum(thisIdx) < 1
            continue
        end

        fixWindowMat(:, ww) = mean(fixZ(:, thisIdx), 2, 'omitnan');
        mixWindowMat(:, ww) = mean(mixZ(:, thisIdx), 2, 'omitnan');

    end

    observedT = nan(1, nWindows);
    observedDiff = nan(1, nWindows);
    fixWindowMean = nan(1, nWindows);
    mixWindowMean = nan(1, nWindows);

    for ww = 1:nWindows

        xFix = fixWindowMat(:, ww);
        xMix = mixWindowMat(:, ww);

        fixWindowMean(ww) = mean(xFix, 'omitnan');
        mixWindowMean(ww) = mean(xMix, 'omitnan');

        observedDiff(ww) = mixWindowMean(ww) - fixWindowMean(ww);
        observedT(ww) = calc_abs_t(xFix, xMix);

    end

    allWindowMat = [fixWindowMat; mixWindowMat];

    nFixHere = size(fixWindowMat, 1);
    nMixHere = size(mixWindowMat, 1);
    nTotalHere = nFixHere + nMixHere;

    groupLabels = [ones(nFixHere, 1); 2 .* ones(nMixHere, 1)];

    permT = nan(nPerm, nWindows);

    for pp = 1:nPerm

        shuffledLabels = groupLabels(randperm(nTotalHere));

        for ww = 1:nWindows

            y = allWindowMat(:, ww);

            xFixPerm = y(shuffledLabels == 1);
            xMixPerm = y(shuffledLabels == 2);

            permT(pp, ww) = calc_abs_t(xFixPerm, xMixPerm);

        end

    end

    permP = nan(1, nWindows);

    for ww = 1:nWindows

        validPerm = permT(:, ww);
        validPerm = validPerm(~isnan(validPerm));

        if isempty(validPerm) || isnan(observedT(ww))
            continue
        end

        permP(ww) = (sum(validPerm >= observedT(ww)) + 1) / (length(validPerm) + 1);

    end

    fdrP = bh_fdr(permP);
    fdrSig = fdrP < alphaPerm;

    permThresholds = nan(1, nWindows);

    for ww = 1:nWindows

        validPerm = permT(:, ww);
        validPerm = validPerm(~isnan(validPerm));

        if isempty(validPerm)
            continue
        end

        permThresholds(ww) = prctile(validPerm, 100 * (1 - alphaPerm));

    end

    sigWindowIdx = observedT > permThresholds;
    sigWindowIdx(isnan(sigWindowIdx)) = false;

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

            clusterMass = sum(observedT(clusterStartIdx:clusterEndIdx), 'omitnan');

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

        permSigWindowIdx = permT(pp, :) > permThresholds;
        permSigWindowIdx(isnan(permSigWindowIdx)) = false;

        permClusterMasses = [];

        ww = 1;

        while ww <= nWindows

            if permSigWindowIdx(ww)

                clusterStartIdx = ww;

                while ww <= nWindows && permSigWindowIdx(ww)
                    ww = ww + 1;
                end

                clusterEndIdx = ww - 1;

                clusterMass = sum(permT(pp, clusterStartIdx:clusterEndIdx), 'omitnan');
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
    clusterPvals = [];
    clusterMassVals = [];
    clusterSigWindowIdx = false(1, nWindows);

    for cc = 1:length(observedClusterMasses)

        clusterP = (sum(maxPermClusterMass >= observedClusterMasses(cc)) + 1) / (nPerm + 1);

        if clusterP < alphaPerm

            sigStart = observedClusterStarts(cc);
            sigEnd   = observedClusterEnds(cc);

            sigStart = max(sigStart, -windowBeforeSec);
            sigEnd   = min(sigEnd, 0);

            idx1 = observedClusterStartIdx(cc);
            idx2 = observedClusterEndIdx(cc);

            clusterFixVals = mean(fixWindowMat(:, idx1:idx2), 2, 'omitnan');
            clusterMixVals = mean(mixWindowMat(:, idx1:idx2), 2, 'omitnan');

            meanFix = mean(clusterFixVals, 'omitnan');
            meanMix = mean(clusterMixVals, 'omitnan');

            if meanMix > meanFix
                labelTxt = 'Mix > Fix';
            else
                labelTxt = 'Fix > Mix';
            end

            sigSegments = [sigSegments; sigStart sigEnd];
            sigLabels{end+1, 1} = labelTxt;

            clusterPvals = [clusterPvals; clusterP];
            clusterMassVals = [clusterMassVals; observedClusterMasses(cc)];

            clusterSigWindowIdx(idx1:idx2) = true;

        end

    end

    testOut.hasTest = true;
    testOut.fixWindowMean = fixWindowMean;
    testOut.mixWindowMean = mixWindowMean;
    testOut.observedDiff = observedDiff;
    testOut.observedT = observedT;
    testOut.permP = permP;
    testOut.fdrP = fdrP;
    testOut.fdrSig = fdrSig;
    testOut.clusterSigWindowIdx = clusterSigWindowIdx;

    testOut.sigSegments = sigSegments;
    testOut.sigLabels = sigLabels;
    testOut.clusterP = clusterPvals;
    testOut.clusterMass = clusterMassVals;

end


function tVal = calc_abs_t(x1, x2)

    x1 = x1(~isnan(x1));
    x2 = x2(~isnan(x2));

    if length(x1) < 2 || length(x2) < 2
        tVal = NaN;
        return
    end

    m1 = mean(x1, 'omitnan');
    m2 = mean(x2, 'omitnan');

    v1 = var(x1, 0, 'omitnan');
    v2 = var(x2, 0, 'omitnan');

    n1 = length(x1);
    n2 = length(x2);

    denom = sqrt((v1 / n1) + (v2 / n2));

    if denom == 0 || isnan(denom)
        tVal = 0;
    else
        tVal = abs((m1 - m2) / denom);
    end

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