clc;
clear;
close all;

rng('shuffle');

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};

inputFolder = '\\155.100.91.44\d\Data\Nill\starling\spikes\spike_data\';
OutputFolder = '\\155.100.91.44\d\Code\Nill\Starling_units_analysis\starling_7_FR_outcome_fix_mix\';

if ~exist(OutputFolder, 'dir')
    mkdir(OutputFolder);
end

windowBeforeSec = 0.05;
windowAfterSec  = 2;
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

baselineIdx = timeCenters >= -windowBeforeSec & timeCenters < 0;

outcomeNames  = ["win", "lose"];
outcomeTitles = ["Win", "Lose"];

winColor  = [0.10 0.45 0.80];
loseColor = [0.80 0.20 0.20];

outcomeColors = [
    winColor
    loseColor
];

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

    block = to_numeric_vector(bhvData.block);

    outcome = bhvData.outcome;
    outcome = outcome(:);

    if isnumeric(outcome) || islogical(outcome)

        outcomeNum = double(outcome);

        if any(outcomeNum < 0, 'all')
            isWinTrial = outcomeNum > 0;
            isLoseTrial = outcomeNum < 0;
        else
            isWinTrial = outcomeNum == 1;
            isLoseTrial = outcomeNum == 0;
        end

    else

        outcomeStr = lower(strtrim(string(outcome)));
        outcomeNum = str2double(outcomeStr);

        if any(~isnan(outcomeNum))

            if any(outcomeNum < 0, 'all')
                isWinTrial = outcomeNum > 0;
                isLoseTrial = outcomeNum < 0;
            else
                isWinTrial = outcomeNum == 1;
                isLoseTrial = outcomeNum == 0;
            end

        else

            isWinTrial = ...
                outcomeStr == "win" | ...
                outcomeStr == "won" | ...
                outcomeStr == "correct" | ...
                outcomeStr == "reward" | ...
                outcomeStr == "rewarded" | ...
                outcomeStr == "gain" | ...
                contains(outcomeStr, "win") | ...
                contains(outcomeStr, "correct") | ...
                contains(outcomeStr, "reward");

            isLoseTrial = ...
                outcomeStr == "lose" | ...
                outcomeStr == "loss" | ...
                outcomeStr == "lost" | ...
                outcomeStr == "incorrect" | ...
                outcomeStr == "wrong" | ...
                outcomeStr == "no_reward" | ...
                contains(outcomeStr, "lose") | ...
                contains(outcomeStr, "loss") | ...
                contains(outcomeStr, "incorrect") | ...
                contains(outcomeStr, "wrong");

        end

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

        fixRasterX = cell(1, 2);
        fixRasterY = cell(1, 2);
        mixRasterX = cell(1, 2);
        mixRasterY = cell(1, 2);

        fixTrialFR = cell(1, 2);
        mixTrialFR = cell(1, 2);

        fixTrialSpikeCounts = cell(1, 2);
        mixTrialSpikeCounts = cell(1, 2);

        validFixTrials = zeros(1, 2);
        validMixTrials = zeros(1, 2);

        for oo = 1:2

            if outcomeNames(oo) == "win"
                thisOutcomeMask = isWinTrial;
            else
                thisOutcomeMask = isLoseTrial;
            end

            fixTrials = find(thisOutcomeMask & ismember(block, [1 2 3]));
            mixTrials = find(thisOutcomeMask & block == 4);

            [fixRasterX{oo}, fixRasterY{oo}, fixTrialFR{oo}, fixTrialSpikeCounts{oo}, validFixTrials(oo)] = ...
                collect_trials_for_condition( ...
                    fixTrials, ...
                    alignTime, ...
                    unitSpikeTimes, ...
                    SampleRes, ...
                    timeEdges, ...
                    windowBeforeSec, ...
                    windowAfterSec, ...
                    binSizeSec);

            [mixRasterX{oo}, mixRasterY{oo}, mixTrialFR{oo}, mixTrialSpikeCounts{oo}, validMixTrials(oo)] = ...
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

        for oo = 1:2

            if ~isempty(fixTrialFR{oo})
                allTrialFR = [allTrialFR; fixTrialFR{oo}];
                allTrialSpikeCounts = [allTrialSpikeCounts; fixTrialSpikeCounts{oo}];
            end

            if ~isempty(mixTrialFR{oo})
                allTrialFR = [allTrialFR; mixTrialFR{oo}];
                allTrialSpikeCounts = [allTrialSpikeCounts; mixTrialSpikeCounts{oo}];
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

        fixZ = cell(1, 2);
        mixZ = cell(1, 2);

        for oo = 1:2

            if ~isempty(fixTrialFR{oo})
                fixZ{oo} = (fixTrialFR{oo} - baselineMean) ./ baselineStd;
            else
                fixZ{oo} = [];
            end

            if ~isempty(mixTrialFR{oo})
                mixZ{oo} = (mixTrialFR{oo} - baselineMean) ./ baselineStd;
            else
                mixZ{oo} = [];
            end

        end

        testOut = cell(1, 2);

        for oo = 1:2

            testOut{oo} = run_fix_mix_cluster_test( ...
                fixZ{oo}, ...
                mixZ{oo}, ...
                timeCenters, ...
                windowBeforeSec, ...
                windowAfterSec, ...
                permWindowSec, ...
                permStrideSec, ...
                nPerm, ...
                alphaPerm);

            nWindows = length(testOut{oo}.permStartTimes);

            tmpStats = table();

            tmpStats.Patient = repmat(string(ptID), nWindows, 1);
            tmpStats.Area = repmat(string(areaNameClean), nWindows, 1);
            tmpStats.Channel = repmat(chanNum, nWindows, 1);
            tmpStats.Unit = repmat(unitNum, nWindows, 1);
            tmpStats.Outcome = repmat(outcomeNames(oo), nWindows, 1);
            tmpStats.FixTrials_n = repmat(validFixTrials(oo), nWindows, 1);
            tmpStats.MixTrials_n = repmat(validMixTrials(oo), nWindows, 1);
            tmpStats.WindowStartSec = testOut{oo}.permStartTimes(:);
            tmpStats.WindowEndSec = testOut{oo}.permEndTimes(:);
            tmpStats.FixMean = testOut{oo}.fixWindowMean(:);
            tmpStats.MixMean = testOut{oo}.mixWindowMean(:);
            tmpStats.MixMinusFix = testOut{oo}.observedDiff(:);
            tmpStats.AbsT = testOut{oo}.observedT(:);
            tmpStats.Perm_p = testOut{oo}.permP(:);
            tmpStats.FDR_p = testOut{oo}.fdrP(:);
            tmpStats.FDR_sig = testOut{oo}.fdrSig(:);
            tmpStats.ClusterCorrected_sig = testOut{oo}.clusterSigWindowIdx(:);

            allStats = [allStats; tmpStats];

            if ~isempty(testOut{oo}.clusterP)

                for cc = 1:length(testOut{oo}.clusterP)

                    tmpCluster = table();

                    tmpCluster.Patient = string(ptID);
                    tmpCluster.Area = string(areaNameClean);
                    tmpCluster.Channel = chanNum;
                    tmpCluster.Unit = unitNum;
                    tmpCluster.Outcome = outcomeNames(oo);
                    tmpCluster.FixTrials_n = validFixTrials(oo);
                    tmpCluster.MixTrials_n = validMixTrials(oo);
                    tmpCluster.ClusterStartSec = testOut{oo}.sigSegments(cc, 1);
                    tmpCluster.ClusterEndSec = testOut{oo}.sigSegments(cc, 2);
                    tmpCluster.ClusterMass = testOut{oo}.clusterMass(cc);
                    tmpCluster.ClusterP = testOut{oo}.clusterP(cc);
                    tmpCluster.Direction = string(testOut{oo}.sigLabels{cc});

                    allClusterStats = [allClusterStats; tmpCluster];

                end

            end

        end

        fig = figure('Visible', 'off', 'Color', 'w');
        set(fig, 'Position', [100 100 1250 850]);

        tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        for oo = 1:2

            thisColor = outcomeColors(oo, :);

            nexttile(oo)
            hold on

            hFixRaster = scatter(fixRasterX{oo}, fixRasterY{oo}, 6, ...
                'filled', ...
                'MarkerFaceColor', thisColor, ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', fixRasterAlpha, ...
                'MarkerEdgeAlpha', fixRasterAlpha);

            hMixRaster = scatter(mixRasterX{oo}, mixRasterY{oo} + validFixTrials(oo), 6, ...
                'filled', ...
                'MarkerFaceColor', thisColor, ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', mixRasterAlpha, ...
                'MarkerEdgeAlpha', mixRasterAlpha);

            xline(0, '--k', 'HandleVisibility', 'off')

            if validFixTrials(oo) > 0 && validMixTrials(oo) > 0
                yline(validFixTrials(oo) + 0.5, '--k', 'HandleVisibility', 'off')
            end

            totalRasterTrials = validFixTrials(oo) + validMixTrials(oo);

            if totalRasterTrials > 0
                ylim([0.5 totalRasterTrials + 0.5])
            end

            xlim([-windowBeforeSec windowAfterSec])

            title(sprintf('%s raster | fix n=%d | mix n=%d', ...
                outcomeTitles(oo), validFixTrials(oo), validMixTrials(oo)), ...
                'Interpreter', 'none')

            xlabel(sprintf('time from %s (s)', alignField))

            if oo == 1
                ylabel('trials')
            else
                ylabel('')
            end

            legend([hFixRaster hMixRaster], ...
                {'fix', 'mix'}, ...
                'Location', 'best')

            box off

        end

        for oo = 1:2

            thisColor = outcomeColors(oo, :);

            [fixMeanSmooth, fixSEMSmooth] = mean_sem_smooth(fixZ{oo}, smoothBins);
            [mixMeanSmooth, mixSEMSmooth] = mean_sem_smooth(mixZ{oo}, smoothBins);

            nexttile(oo + 2)
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

            for ss = 1:size(testOut{oo}.sigSegments, 1)

                sigStart = testOut{oo}.sigSegments(ss, 1);
                sigEnd   = testOut{oo}.sigSegments(ss, 2);

                plot([sigStart sigEnd], [sigY sigY], ...
                    '-', ...
                    'Color', [0.15 0.15 0.15], ...
                    'LineWidth', 4, ...
                    'HandleVisibility', 'off')

                text(mean([sigStart sigEnd]), sigY + 0.08 * yRange, testOut{oo}.sigLabels{ss}, ...
                    'Color', [0.05 0.05 0.05], ...
                    'FontSize', 9, ...
                    'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom')

            end

            if ~isempty(testOut{oo}.sigSegments)
                ylim([curveLow - 0.15 * yRange, sigY + 0.55 * yRange])
            else
                ylim([curveLow - 0.15 * yRange, curveHigh + 0.25 * yRange])
            end

            xlim([-windowBeforeSec windowAfterSec])

            title(sprintf('%s PSTH: fix vs mix', outcomeTitles(oo)), ...
                'Interpreter', 'none')

            xlabel(sprintf('time from %s (s)', alignField))

            if oo == 1
                ylabel('baseline z-scored firing rate')
            else
                ylabel('')
            end

            if ~isempty(legendHandles)
                legend(legendHandles, legendLabels, 'Location', 'best')
            end

            box off

        end

        sgtitle(sprintf('%s | %s | Ch %d Unit %d | FR %.2f Hz | outcome fix vs mix ', ...
            ptID, areaNameClean, chanNum, unitNum, meanFRHz), ...
            'Interpreter', 'none', ...
            'FontWeight', 'bold')

        pdfName = sprintf('%s_%s_ch%d_unit%d_choice_feedback_win_lose_fix_vs_mix.pdf', ...
            ptID, areaNameClean, chanNum, unitNum);

        pdfPath = fullfile(OutputFolder, pdfName);

        exportgraphics(fig, pdfPath, 'ContentType', 'vector');
        close(fig);

        fprintf('saved: %s \n', pdfName);

    end

end

statsPath = fullfile(OutputFolder, 'choice_feedback_win_lose_fix_vs_mix_window_stats.csv');
writetable(allStats, statsPath);

clusterStatsPath = fullfile(OutputFolder, 'choice_feedback_win_lose_fix_vs_mix_cluster_stats.csv');
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
    fixZ, mixZ, timeCenters, windowBeforeSec, windowAfterSec, permWindowSec, permStrideSec, nPerm, alphaPerm)

    permStartTimes = -windowBeforeSec:permStrideSec:(windowAfterSec - permWindowSec);
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
            sigEnd   = min(sigEnd, windowAfterSec);

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