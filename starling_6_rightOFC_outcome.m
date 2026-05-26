clc;
clear;
close all;

% ========================================================================
% Starling task: RIGHT OFC outcome encoding analysis
%
% This script focuses only on outcome encoding in right OFC.
%
% It creates:
%   1) Per-unit win vs lose outcome PSTH/raster with cluster permutation test
%   2) Right OFC population outcome PSTH using unit as the observation
%   3) Outcome modulation index per unit: mean(win) - mean(lose), 0-1 s
%   4) AUROC per unit for win vs lose from 0-1 s firing rate
%   5) Binomial prevalence test: is # significant units > chance at 5%?
%   6) Participant-level robustness summary
%   7) Leave-one-participant-out robustness summary
%   8) Trial-level regression per unit:
%          FR_0_1s ~ outcome + choice + distribution + myCard + trialIndex
%
% IMPORTANT:
%   - eventTime.choiceAndFeedbackTime is used as outcome onset.
%   - FR is baseline z-scored using -1 to -0.25 s before outcome onset.
%   - Per-unit statistics use cluster-based permutation correction over time.
%   - Population statistics use sign-flip cluster permutation across units.
%   - Units are pooled descriptively across participants, but participant-level
%     and leave-one-participant-out summaries are included.
% ========================================================================

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};

inputFolder = '\\155.100.91.44\d\Data\Nill\starling\spikes\spike_data\';
OutputFolder = '\\155.100.91.44\d\Code\Nill\Starling_units_analysis\starling_6_rightOFC_outcome\';

if ~exist(OutputFolder, 'dir')
    mkdir(OutputFolder);
end

figFolder = fullfile(OutputFolder, 'Figures');
unitFigFolder = fullfile(OutputFolder, 'Unit_Outcome_Figures');
tableFolder = fullfile(OutputFolder, 'Tables');

if ~exist(figFolder, 'dir'); mkdir(figFolder); end
if ~exist(unitFigFolder, 'dir'); mkdir(unitFigFolder); end
if ~exist(tableFolder, 'dir'); mkdir(tableFolder); end

% ============================= settings =================================
windowBeforeSec = 1;
windowAfterSec  = 2;
binSizeSec = 0.05;
smoothBins = 7;

baselineWindowSec = [-1 -0.25];
outcomeSummaryWindowSec = [0 1];
populationTestWindowSec = [0 2];

minMeanFRHz = 0.5;
minFracTrialsWithSpikes = 0.10;
minTrialsPerOutcome = 5;

permWindowMs = 200;
permStrideMs = 50;
nPerm = 1000;
alphaPerm = 0.05;
clusterFormingAlpha = 0.05;

saveAllUnitFigures = true;
% If this makes too many PDFs, set saveAllUnitFigures = false to save only
% units with corrected significant outcome clusters.

rng(1);
% ========================================================================

permWindowSec = permWindowMs / 1000;
permStrideSec = permStrideMs / 1000;
analysisWindowSec = windowBeforeSec + windowAfterSec;

timeEdges = -windowBeforeSec:binSizeSec:windowAfterSec;
timeCenters = timeEdges(1:end-1) + binSizeSec/2;
baselineIdx = timeCenters >= baselineWindowSec(1) & timeCenters < baselineWindowSec(2);
outcomeSummaryIdx = timeCenters >= outcomeSummaryWindowSec(1) & timeCenters < outcomeSummaryWindowSec(2);

winColor  = [0 0.60 0];
loseColor = [0.85 0 0];

allUnitSummary = table();
allClusterStats = table();
allRegressionSummary = table();

% For right OFC population plots/statistics. Each row = one unit.
popWinMean = [];
popLoseMean = [];
popDiffMean = [];
popUnitIDs = strings(0,1);
popPatientIDs = strings(0,1);
popAreas = strings(0,1);
popIsSig = [];
popOutcomeDelta = [];
popAUC = [];

% ========================================================================
% Main loop
% ========================================================================

for pt = 1:length(microPts)

    ptID = microPts{pt};
    fprintf('\nProcessing patient %s\n', ptID);

    dataPath = fullfile(inputFolder, sprintf('%s_spikeData.mat', ptID));

    if ~exist(dataPath, 'file')
        warning('Missing file: %s. Skipping patient.', dataPath);
        continue
    end

    load(dataPath, 'spikeData');

    ChanUnitTimestamp = spikeData.ChanUnitTimestamp;
    eventTime = spikeData.eventTimes;
    bhvData = eventTime.bhvData;
    inclChans = spikeData.inclChans;
    microLabels = spikeData.microLabels;
    SampleRes = double(spikeData.SampleRes);

    if ~isfield(eventTime, 'choiceAndFeedbackTime')
        warning('Patient %s has no eventTime.choiceAndFeedbackTime. Skipping.', ptID);
        continue
    end

    outcomeTime = double(eventTime.choiceAndFeedbackTime(:));

    outcome = get_string_field(bhvData, 'outcome');
    winTrialsAll  = find(outcome == "win" | outcome == "correct" | outcome == "reward" | outcome == "1");
    loseTrialsAll = find(outcome == "lose" | outcome == "loss" | outcome == "incorrect" | outcome == "0");

    if isempty(winTrialsAll) || isempty(loseTrialsAll)
        warning('Patient %s does not have both win and lose trials. Skipping.', ptID);
        continue
    end

    allChanUnits = unique(ChanUnitTimestamp(:,1:2), 'rows');
    allChanUnits = double(allChanUnits);

    % Remove unit 255 because it means no waveform was saved.
    allChanUnits(allChanUnits(:,2) == 255, :) = [];

    % Remove channels greater than included channels.
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

        areaNameRaw = char(areaName);
        areaNameClean = regexprep(areaNameRaw, '[^\w]', '_');

        % Focus only on right OFC.
        if ~is_right_ofc(areaNameRaw)
            continue
        end

        unitSpikeTimes = ChanUnitTimestamp( ...
            double(ChanUnitTimestamp(:,1)) == chanNum & ...
            double(ChanUnitTimestamp(:,2)) == unitNum, 3);

        unitSpikeTimes = double(unitSpikeTimes);

        if isempty(unitSpikeTimes)
            continue
        end

        unitID = string(sprintf('%s_%s_ch%d_u%d', ptID, areaNameClean, chanNum, unitNum));

        % ---------------------- Get trial firing rates ------------------
        [winRasterX, winRasterY, winTrialFR, winTrialSpikeCounts, validWinTrials, usedWinTrialIDs] = ...
            get_trial_fr(winTrialsAll, outcomeTime, unitSpikeTimes, SampleRes, ...
            windowBeforeSec, windowAfterSec, timeEdges, binSizeSec);

        [loseRasterX, loseRasterY, loseTrialFR, loseTrialSpikeCounts, validLoseTrials, usedLoseTrialIDs] = ...
            get_trial_fr(loseTrialsAll, outcomeTime, unitSpikeTimes, SampleRes, ...
            windowBeforeSec, windowAfterSec, timeEdges, binSizeSec);

        if validWinTrials < minTrialsPerOutcome || validLoseTrials < minTrialsPerOutcome
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

        % ---------------------- baseline z-score ------------------------
        allTrialFR = [winTrialFR; loseTrialFR];
        baselineVals = allTrialFR(:, baselineIdx);
        baselineMean = mean(baselineVals(:), 'omitnan');
        baselineStd = std(baselineVals(:), 'omitnan');

        if baselineStd == 0 || isnan(baselineStd)
            baselineStd = 1;
        end

        winZ  = (winTrialFR  - baselineMean) ./ baselineStd;
        loseZ = (loseTrialFR - baselineMean) ./ baselineStd;

        winMean = mean(winZ, 1, 'omitnan');
        loseMean = mean(loseZ, 1, 'omitnan');
        winSEM = std(winZ, 0, 1, 'omitnan') ./ sqrt(size(winZ, 1));
        loseSEM = std(loseZ, 0, 1, 'omitnan') ./ sqrt(size(loseZ, 1));

        winMeanSmooth = smoothdata(winMean, 'gaussian', smoothBins);
        loseMeanSmooth = smoothdata(loseMean, 'gaussian', smoothBins);
        winSEMSmooth = smoothdata(winSEM, 'gaussian', smoothBins);
        loseSEMSmooth = smoothdata(loseSEM, 'gaussian', smoothBins);

        % ---------------------- per-unit cluster test -------------------
        [sigSegments, clusterTable, observedDiffs, windowStartTimes, windowEndTimes] = ...
            run_two_group_cluster_permutation( ...
            winZ, loseZ, timeCenters, [0 windowAfterSec], ...
            permWindowSec, permStrideSec, nPerm, alphaPerm);

        isSignificant = ~isempty(sigSegments);

        if isempty(clusterTable)
            minClusterP = NaN;
            firstSigStartSec = NaN;
            firstSigEndSec = NaN;
        else
            minClusterP = min(clusterTable.ClusterP, [], 'omitnan');
            sigRows = find(clusterTable.ClusterSignificant);
            if isempty(sigRows)
                firstSigStartSec = NaN;
                firstSigEndSec = NaN;
            else
                firstSigStartSec = clusterTable.ClusterStartSec(sigRows(1));
                firstSigEndSec = clusterTable.ClusterEndSec(sigRows(1));
            end
        end

        % ---------------------- effect size summaries -------------------
        winFR_0_1 = mean(winZ(:, outcomeSummaryIdx), 2, 'omitnan');
        loseFR_0_1 = mean(loseZ(:, outcomeSummaryIdx), 2, 'omitnan');
        outcomeDelta = mean(winFR_0_1, 'omitnan') - mean(loseFR_0_1, 'omitnan');

        if outcomeDelta > 0
            preferredOutcome = "win_preferring";
        elseif outcomeDelta < 0
            preferredOutcome = "lose_preferring";
        else
            preferredOutcome = "no_difference";
        end

        trialFR_0_1 = [winFR_0_1; loseFR_0_1];
        outcomeLabels = [ones(length(winFR_0_1), 1); zeros(length(loseFR_0_1), 1)];
        aucWinLose = compute_auc(trialFR_0_1, outcomeLabels);

        % Peak time of absolute win-lose difference after outcome.
        diffTrace = winMean - loseMean;
        afterIdx = timeCenters >= 0 & timeCenters < windowAfterSec;
        [peakAbsDiff, peakLocalIdx] = max(abs(diffTrace(afterIdx)), [], 'omitnan');
        afterTimes = timeCenters(afterIdx);
        if isempty(peakLocalIdx) || isnan(peakAbsDiff)
            peakTimeSec = NaN;
        else
            peakTimeSec = afterTimes(peakLocalIdx);
        end

        % ---------------------- regression ------------------------------
        usedTrialIDs = [usedWinTrialIDs; usedLoseTrialIDs];
        [regResult, designInfo] = run_outcome_regression_for_unit( ...
            trialFR_0_1, outcomeLabels, usedTrialIDs, bhvData);

        % ---------------------- save unit summary -----------------------
        tmpSummary = table();
        tmpSummary.Patient = string(ptID);
        tmpSummary.Area = string(areaNameClean);
        tmpSummary.Channel = chanNum;
        tmpSummary.Unit = unitNum;
        tmpSummary.UnitID = unitID;
        tmpSummary.NWinTrials = validWinTrials;
        tmpSummary.NLoseTrials = validLoseTrials;
        tmpSummary.NTrialsTotal = nValidTrialsTotal;
        tmpSummary.MeanFRHz = meanFRHz;
        tmpSummary.FracTrialsWithSpikes = fracTrialsWithSpikes;
        tmpSummary.IsOutcomeSignificant = isSignificant;
        tmpSummary.NSignificantClusters = size(sigSegments, 1);
        tmpSummary.MinClusterP = minClusterP;
        tmpSummary.FirstSigStartSec = firstSigStartSec;
        tmpSummary.FirstSigEndSec = firstSigEndSec;
        tmpSummary.OutcomeDelta_WinMinusLose_0to1s = outcomeDelta;
        tmpSummary.PreferredOutcome = preferredOutcome;
        tmpSummary.AUROC_WinVsLose_0to1s = aucWinLose;
        tmpSummary.PeakAbsDiff_0to2s = peakAbsDiff;
        tmpSummary.PeakAbsDiffTimeSec = peakTimeSec;

        allUnitSummary = [allUnitSummary; tmpSummary]; %#ok<AGROW>

        if ~isempty(clusterTable)
            nRows = height(clusterTable);
            clusterTable.Patient = repmat(string(ptID), nRows, 1);
            clusterTable.Area = repmat(string(areaNameClean), nRows, 1);
            clusterTable.Channel = repmat(chanNum, nRows, 1);
            clusterTable.Unit = repmat(unitNum, nRows, 1);
            clusterTable.UnitID = repmat(unitID, nRows, 1);
            clusterTable = movevars(clusterTable, ...
                {'Patient', 'Area', 'Channel', 'Unit', 'UnitID'}, 'Before', 1);
            allClusterStats = [allClusterStats; clusterTable]; %#ok<AGROW>
        end

        tmpReg = table();
        tmpReg.Patient = string(ptID);
        tmpReg.Area = string(areaNameClean);
        tmpReg.Channel = chanNum;
        tmpReg.Unit = unitNum;
        tmpReg.UnitID = unitID;
        tmpReg.NTrialsUsed = regResult.NTrialsUsed;
        tmpReg.NPredictors = regResult.NPredictors;
        tmpReg.OutcomeCoef = regResult.OutcomeCoef;
        tmpReg.OutcomeT = regResult.OutcomeT;
        tmpReg.OutcomeP = regResult.OutcomeP;
        tmpReg.ModelR2 = regResult.R2;
        tmpReg.DesignColumns = string(strjoin(designInfo.columnNames, ', '));
        allRegressionSummary = [allRegressionSummary; tmpReg]; %#ok<AGROW>

        % ---------------------- store for population --------------------
        popWinMean = [popWinMean; winMean]; %#ok<AGROW>
        popLoseMean = [popLoseMean; loseMean]; %#ok<AGROW>
        popDiffMean = [popDiffMean; diffTrace]; %#ok<AGROW>
        popUnitIDs(end+1,1) = unitID; %#ok<AGROW>
        popPatientIDs(end+1,1) = string(ptID); %#ok<AGROW>
        popAreas(end+1,1) = string(areaNameClean); %#ok<AGROW>
        popIsSig(end+1,1) = isSignificant; %#ok<AGROW>
        popOutcomeDelta(end+1,1) = outcomeDelta; %#ok<AGROW>
        popAUC(end+1,1) = aucWinLose; %#ok<AGROW>

        % ---------------------- unit figure -----------------------------
        if saveAllUnitFigures || isSignificant
            save_unit_outcome_figure(unitFigFolder, ptID, areaNameClean, chanNum, unitNum, ...
                timeCenters, winRasterX, winRasterY, loseRasterX, loseRasterY, ...
                validWinTrials, winMeanSmooth, loseMeanSmooth, winSEMSmooth, loseSEMSmooth, ...
                sigSegments, winColor, loseColor, windowBeforeSec, windowAfterSec, ...
                meanFRHz, outcomeDelta, aucWinLose, preferredOutcome);
        end

        fprintf('Right OFC unit saved: %s | sig=%d | delta=%.3f | AUC=%.3f\n', ...
            unitID, isSignificant, outcomeDelta, aucWinLose);
    end
end

% ========================================================================
% Save raw unit-level tables
% ========================================================================

if isempty(allUnitSummary)
    warning('No valid right OFC units found. Check area labels and filters.');
    return
end

% FDR correction for regression outcome p-values across right OFC units.
allRegressionSummary.OutcomeQ_BHFDR = bh_fdr(allRegressionSummary.OutcomeP);
allRegressionSummary.OutcomeSignificant_FDR = allRegressionSummary.OutcomeQ_BHFDR < 0.05;

writetable(allUnitSummary, fullfile(tableFolder, 'rightOFC_outcome_unit_summary.csv'));
writetable(allClusterStats, fullfile(tableFolder, 'rightOFC_outcome_cluster_stats.csv'));
writetable(allRegressionSummary, fullfile(tableFolder, 'rightOFC_outcome_regression_summary.csv'));

% ========================================================================
% Summary statistics
% ========================================================================

nRightOFCUnits = height(allUnitSummary);
nSigOutcomeUnits = sum(allUnitSummary.IsOutcomeSignificant);
percentSigOutcome = 100 * nSigOutcomeUnits / nRightOFCUnits;
binomialP_5percent = binomial_upper_tail(nSigOutcomeUnits, nRightOFCUnits, 0.05);

nWinPrefSig = sum(allUnitSummary.IsOutcomeSignificant & allUnitSummary.OutcomeDelta_WinMinusLose_0to1s > 0);
nLosePrefSig = sum(allUnitSummary.IsOutcomeSignificant & allUnitSummary.OutcomeDelta_WinMinusLose_0to1s < 0);

fprintf('\n============================================================\n');
fprintf('RIGHT OFC OUTCOME ENCODING SUMMARY\n');
fprintf('Valid right OFC units: %d\n', nRightOFCUnits);
fprintf('Outcome significant units: %d/%d = %.1f%%\n', nSigOutcomeUnits, nRightOFCUnits, percentSigOutcome);
fprintf('Binomial p against 5%% chance: %.6g\n', binomialP_5percent);
fprintf('Significant win-preferring units: %d\n', nWinPrefSig);
fprintf('Significant lose-preferring units: %d\n', nLosePrefSig);
fprintf('============================================================\n');

summaryMain = table();
summaryMain.Area = "right_orbitofrontal_cortex";
summaryMain.NValidUnits = nRightOFCUnits;
summaryMain.NOutcomeSignificantUnits = nSigOutcomeUnits;
summaryMain.PercentOutcomeSignificant = percentSigOutcome;
summaryMain.BinomialP_Against_5Percent = binomialP_5percent;
summaryMain.NSignificantWinPreferring = nWinPrefSig;
summaryMain.NSignificantLosePreferring = nLosePrefSig;
writetable(summaryMain, fullfile(tableFolder, 'rightOFC_main_outcome_summary.csv'));

% ========================================================================
% Participant-level summary
% ========================================================================

patientSummary = make_patient_summary(allUnitSummary);
writetable(patientSummary, fullfile(tableFolder, 'rightOFC_outcome_patient_summary.csv'));
save_patient_summary_figure(patientSummary, figFolder);

% ========================================================================
% Leave-one-participant-out summary
% ========================================================================

leaveOneSummary = make_leave_one_patient_out_summary(allUnitSummary);
writetable(leaveOneSummary, fullfile(tableFolder, 'rightOFC_outcome_leave_one_patient_out.csv'));
save_leave_one_out_figure(leaveOneSummary, figFolder);

% ========================================================================
% Population PSTH and sign-flip cluster permutation across units
% ========================================================================

[popSigSegments, popClusterStats] = run_population_signflip_cluster( ...
    popDiffMean, timeCenters, populationTestWindowSec, ...
    permWindowSec, permStrideSec, nPerm, alphaPerm, clusterFormingAlpha);

writetable(popClusterStats, fullfile(tableFolder, 'rightOFC_population_cluster_stats.csv'));

save_population_psth_figure(figFolder, timeCenters, popWinMean, popLoseMean, ...
    popSigSegments, winColor, loseColor, windowBeforeSec, windowAfterSec, nRightOFCUnits);

% ========================================================================
% Modulation/AUROC/regression figures
% ========================================================================

save_outcome_delta_figure(figFolder, allUnitSummary);
save_auc_figure(figFolder, allUnitSummary);
save_regression_figure(figFolder, allRegressionSummary);

fprintf('\nDONE. Results saved in:\n%s\n', OutputFolder);

% ========================================================================
% Local functions
% ========================================================================

function tf = is_right_ofc(areaNameRaw)
    areaLower = lower(regexprep(char(areaNameRaw), '[^a-zA-Z0-9]', ' '));
    tf = contains(areaLower, 'right') && ...
         (contains(areaLower, 'orbitofrontal') || contains(areaLower, 'orbito frontal') || contains(areaLower, 'ofc')) && ...
         ~contains(areaLower, 'left');
end

function [rasterX, rasterY, trialFR, trialSpikeCounts, validTrials, usedTrialIDs] = get_trial_fr( ...
    trials, alignTime, unitSpikeTimes, SampleRes, ...
    windowBeforeSec, windowAfterSec, timeEdges, binSizeSec)

    rasterX = [];
    rasterY = [];
    trialFR = [];
    trialSpikeCounts = [];
    usedTrialIDs = [];
    validTrials = 0;

    for i = 1:length(trials)
        tr = trials(i);

        if tr > length(alignTime) || isnan(alignTime(tr))
            continue
        end

        validTrials = validTrials + 1;
        usedTrialIDs(validTrials, 1) = tr; %#ok<AGROW>

        windowStart = alignTime(tr) - windowBeforeSec * SampleRes;
        windowEnd = alignTime(tr) + windowAfterSec * SampleRes;

        spikesInWindow = unitSpikeTimes(unitSpikeTimes >= windowStart & unitSpikeTimes <= windowEnd);
        relSpikesSec = (spikesInWindow - alignTime(tr)) ./ SampleRes;

        rasterX = [rasterX; relSpikesSec(:)]; %#ok<AGROW>
        rasterY = [rasterY; validTrials .* ones(length(relSpikesSec), 1)]; %#ok<AGROW>

        counts = histcounts(relSpikesSec, timeEdges);
        trialFR(validTrials, :) = counts ./ binSizeSec; %#ok<AGROW>
        trialSpikeCounts(validTrials, 1) = sum(counts); %#ok<AGROW>
    end
end

function [sigSegments, clusterTable, observedDiffs, permStartTimes, permEndTimes] = ...
    run_two_group_cluster_permutation(winZ, loseZ, timeCenters, testWindowSec, ...
    permWindowSec, permStrideSec, nPerm, alphaPerm)

    sigSegments = [];
    clusterTable = table();

    permStartTimes = testWindowSec(1):permStrideSec:(testWindowSec(2) - permWindowSec);
    permEndTimes = permStartTimes + permWindowSec;
    nWindows = length(permStartTimes);

    observedDiffs = nan(1, nWindows);

    if isempty(winZ) || isempty(loseZ) || nWindows < 1
        return
    end

    winWindowMat = nan(size(winZ, 1), nWindows);
    loseWindowMat = nan(size(loseZ, 1), nWindows);

    for ww = 1:nWindows
        thisIdx = timeCenters >= permStartTimes(ww) & timeCenters < permEndTimes(ww);
        if sum(thisIdx) < 1
            continue
        end
        winWindowMat(:, ww) = mean(winZ(:, thisIdx), 2, 'omitnan');
        loseWindowMat(:, ww) = mean(loseZ(:, thisIdx), 2, 'omitnan');
    end

    for ww = 1:nWindows
        observedDiffs(ww) = mean(winWindowMat(:, ww), 'omitnan') - mean(loseWindowMat(:, ww), 'omitnan');
    end

    allWindowMat = [winWindowMat; loseWindowMat];
    nWinHere = size(winWindowMat, 1);
    nLoseHere = size(loseWindowMat, 1);
    nTotalHere = nWinHere + nLoseHere;

    permDiffs = nan(nPerm, nWindows);

    for pp = 1:nPerm
        shuffledIdx = randperm(nTotalHere);
        permWinIdx = shuffledIdx(1:nWinHere);
        permLoseIdx = shuffledIdx(nWinHere+1:end);

        permWinMat = allWindowMat(permWinIdx, :);
        permLoseMat = allWindowMat(permLoseIdx, :);

        for ww = 1:nWindows
            permDiffs(pp, ww) = mean(permWinMat(:, ww), 'omitnan') - mean(permLoseMat(:, ww), 'omitnan');
        end
    end

    % Pointwise cluster-forming threshold from the permutation distribution.
    permThresholds = prctile(abs(permDiffs), 100 * (1 - alphaPerm), 1);
    sigWindowIdx = abs(observedDiffs) > permThresholds;

    [obsStartIdx, obsEndIdx, obsMass] = find_clusters(sigWindowIdx, abs(observedDiffs));

    maxPermClusterMass = zeros(nPerm, 1);

    for pp = 1:nPerm
        permSigWindowIdx = abs(permDiffs(pp, :)) > permThresholds;
        [~, ~, permMass] = find_clusters(permSigWindowIdx, abs(permDiffs(pp, :)));
        if isempty(permMass)
            maxPermClusterMass(pp) = 0;
        else
            maxPermClusterMass(pp) = max(permMass, [], 'omitnan');
        end
    end

    nObsClusters = length(obsMass);
    clusterP = nan(nObsClusters, 1);
    clusterSignificant = false(nObsClusters, 1);
    clusterStartSec = nan(nObsClusters, 1);
    clusterEndSec = nan(nObsClusters, 1);

    for cc = 1:nObsClusters
        clusterP(cc) = (sum(maxPermClusterMass >= obsMass(cc)) + 1) / (nPerm + 1);
        clusterSignificant(cc) = clusterP(cc) < alphaPerm;
        clusterStartSec(cc) = permStartTimes(obsStartIdx(cc));
        clusterEndSec(cc) = permEndTimes(obsEndIdx(cc));

        if clusterSignificant(cc)
            sigSegments = [sigSegments; clusterStartSec(cc) clusterEndSec(cc) clusterP(cc) obsMass(cc)]; %#ok<AGROW>
        end
    end

    if nObsClusters > 0
        clusterTable = table((1:nObsClusters)', clusterStartSec, clusterEndSec, obsMass(:), clusterP, clusterSignificant, ...
            'VariableNames', {'ClusterIndex', 'ClusterStartSec', 'ClusterEndSec', 'ClusterMass', 'ClusterP', 'ClusterSignificant'});
    end
end

function [popSigSegments, popClusterStats] = run_population_signflip_cluster( ...
    diffByUnit, timeCenters, testWindowSec, permWindowSec, permStrideSec, ...
    nPerm, alphaPerm, clusterFormingAlpha)

    popSigSegments = [];
    popClusterStats = table();

    nUnits = size(diffByUnit, 1);
    permStartTimes = testWindowSec(1):permStrideSec:(testWindowSec(2) - permWindowSec);
    permEndTimes = permStartTimes + permWindowSec;
    nWindows = length(permStartTimes);

    if nUnits < 2 || nWindows < 1
        return
    end

    windowDiff = nan(nUnits, nWindows);

    for ww = 1:nWindows
        thisIdx = timeCenters >= permStartTimes(ww) & timeCenters < permEndTimes(ww);
        windowDiff(:, ww) = mean(diffByUnit(:, thisIdx), 2, 'omitnan');
    end

    observedT = nan(1, nWindows);
    observedP = nan(1, nWindows);

    for ww = 1:nWindows
        x = windowDiff(:, ww);
        [observedT(ww), observedP(ww)] = one_sample_t_and_p(x);
    end

    observedSig = observedP < clusterFormingAlpha;
    [obsStartIdx, obsEndIdx, obsMass] = find_clusters(observedSig, abs(observedT));

    maxPermClusterMass = zeros(nPerm, 1);

    for pp = 1:nPerm
        signs = randi([0 1], nUnits, 1);
        signs(signs == 0) = -1;
        permWindowDiff = windowDiff .* signs;

        permT = nan(1, nWindows);
        permP = nan(1, nWindows);

        for ww = 1:nWindows
            [permT(ww), permP(ww)] = one_sample_t_and_p(permWindowDiff(:, ww));
        end

        permSig = permP < clusterFormingAlpha;
        [~, ~, permMass] = find_clusters(permSig, abs(permT));

        if isempty(permMass)
            maxPermClusterMass(pp) = 0;
        else
            maxPermClusterMass(pp) = max(permMass, [], 'omitnan');
        end
    end

    nObsClusters = length(obsMass);
    clusterP = nan(nObsClusters, 1);
    clusterSignificant = false(nObsClusters, 1);
    clusterStartSec = nan(nObsClusters, 1);
    clusterEndSec = nan(nObsClusters, 1);

    for cc = 1:nObsClusters
        clusterP(cc) = (sum(maxPermClusterMass >= obsMass(cc)) + 1) / (nPerm + 1);
        clusterSignificant(cc) = clusterP(cc) < alphaPerm;
        clusterStartSec(cc) = permStartTimes(obsStartIdx(cc));
        clusterEndSec(cc) = permEndTimes(obsEndIdx(cc));

        if clusterSignificant(cc)
            popSigSegments = [popSigSegments; clusterStartSec(cc) clusterEndSec(cc) clusterP(cc) obsMass(cc)]; %#ok<AGROW>
        end
    end

    if nObsClusters > 0
        popClusterStats = table((1:nObsClusters)', clusterStartSec, clusterEndSec, obsMass(:), clusterP, clusterSignificant, ...
            'VariableNames', {'ClusterIndex', 'ClusterStartSec', 'ClusterEndSec', 'ClusterMass', 'ClusterP', 'ClusterSignificant'});
    end
end

function [clusterStartIdx, clusterEndIdx, clusterMasses] = find_clusters(sigIdx, statVals)
    clusterStartIdx = [];
    clusterEndIdx = [];
    clusterMasses = [];

    sigIdx(isnan(sigIdx)) = false;
    statVals(isnan(statVals)) = 0;

    ii = 1;
    n = length(sigIdx);

    while ii <= n
        if sigIdx(ii)
            startIdx = ii;
            while ii <= n && sigIdx(ii)
                ii = ii + 1;
            end
            endIdx = ii - 1;
            clusterStartIdx = [clusterStartIdx; startIdx]; %#ok<AGROW>
            clusterEndIdx = [clusterEndIdx; endIdx]; %#ok<AGROW>
            clusterMasses = [clusterMasses; sum(statVals(startIdx:endIdx), 'omitnan')]; %#ok<AGROW>
        else
            ii = ii + 1;
        end
    end
end

function [tVal, pVal] = one_sample_t_and_p(x)
    x = x(:);
    x = x(~isnan(x));
    n = length(x);

    if n < 2
        tVal = NaN;
        pVal = NaN;
        return
    end

    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');

    if sd == 0 || isnan(sd)
        tVal = NaN;
        pVal = NaN;
        return
    end

    tVal = mu / (sd / sqrt(n));
    pVal = 2 * (1 - tcdf(abs(tVal), n - 1));
end

function [regResult, designInfo] = run_outcome_regression_for_unit(yFR, outcomeLabels, usedTrialIDs, bhvData)
    % Builds a trial-level OLS model:
    %   FR_0_1s ~ outcome + choice + distribution + myCard + trialIndex
    % Missing covariates are automatically skipped.

    y = yFR(:);
    outcomeLabels = outcomeLabels(:);
    usedTrialIDs = usedTrialIDs(:);

    X = ones(length(y), 1);
    columnNames = {'Intercept'};

    X = [X outcomeLabels];
    columnNames{end+1} = 'Outcome_win1_lose0';
    outcomeColumnIndex = 2;

    % Choice: arrowup vs arrowdown if available.
    if has_bhv_field(bhvData, 'choice')
        choiceRaw = get_string_field(bhvData, 'choice');
        choiceThis = choiceRaw(usedTrialIDs);
        choiceUp = double(choiceThis == "arrowup" | choiceThis == "up" | choiceThis == "higher" | choiceThis == "high");
        choiceDown = double(choiceThis == "arrowdown" | choiceThis == "down" | choiceThis == "lower" | choiceThis == "low");
        choiceUp(~(choiceUp == 1 | choiceDown == 1)) = NaN;
        X = [X choiceUp];
        columnNames{end+1} = 'Choice_arrowup';
    end

    % Distribution: use two dummy variables, low and high, with uniform as reference.
    if has_bhv_field(bhvData, 'distribution')
        distRaw = get_bhv_field(bhvData, 'distribution');
        distRaw = distRaw(:);
        distThis = distRaw(usedTrialIDs);
        [lowDummy, highDummy] = distribution_to_dummies(distThis);
        X = [X lowDummy highDummy];
        columnNames{end+1} = 'Distribution_low';
        columnNames{end+1} = 'Distribution_high';
    end

    if has_bhv_field(bhvData, 'myCard')
        myCard = get_numeric_field(bhvData, 'myCard');
        cardThis = myCard(usedTrialIDs);
        cardThis = zscore_safe(cardThis);
        X = [X cardThis];
        columnNames{end+1} = 'MyCard_z';
    end

    if has_bhv_field(bhvData, 'trialIndex')
        trialIndex = get_numeric_field(bhvData, 'trialIndex');
        trialThis = trialIndex(usedTrialIDs);
        trialThis = zscore_safe(trialThis);
        X = [X trialThis];
        columnNames{end+1} = 'TrialIndex_z';
    else
        trialThis = zscore_safe(usedTrialIDs);
        X = [X trialThis];
        columnNames{end+1} = 'TrialNumber_z_from_order';
    end

    validRows = ~isnan(y) & all(~isnan(X), 2);
    y = y(validRows);
    X = X(validRows, :);

    % Remove columns with zero variance except intercept.
    keepCols = true(1, size(X, 2));
    for cc = 2:size(X, 2)
        if std(X(:, cc), 0, 'omitnan') == 0
            keepCols(cc) = false;
        end
    end

    X = X(:, keepCols);
    keptNames = columnNames(keepCols);

    outcomeColumnAfterKeep = find(strcmp(keptNames, 'Outcome_win1_lose0'), 1);

    regResult.NTrialsUsed = length(y);
    regResult.NPredictors = size(X, 2) - 1;
    regResult.OutcomeCoef = NaN;
    regResult.OutcomeT = NaN;
    regResult.OutcomeP = NaN;
    regResult.R2 = NaN;

    designInfo.columnNames = keptNames;

    if length(y) <= size(X, 2) + 1 || isempty(outcomeColumnAfterKeep)
        return
    end

    [beta, se, tVals, pVals, R2] = local_ols(y, X);

    regResult.OutcomeCoef = beta(outcomeColumnAfterKeep);
    regResult.OutcomeT = tVals(outcomeColumnAfterKeep);
    regResult.OutcomeP = pVals(outcomeColumnAfterKeep);
    regResult.R2 = R2;
end

function [beta, se, tVals, pVals, R2] = local_ols(y, X)
    y = y(:);
    beta = X \ y;
    yHat = X * beta;
    resid = y - yHat;

    n = length(y);
    p = size(X, 2);
    df = n - p;

    SSE = sum(resid.^2, 'omitnan');
    SST = sum((y - mean(y, 'omitnan')).^2, 'omitnan');

    if SST == 0 || df <= 0
        R2 = NaN;
        se = nan(size(beta));
        tVals = nan(size(beta));
        pVals = nan(size(beta));
        return
    end

    sigma2 = SSE / df;
    covBeta = sigma2 * pinv(X' * X);
    se = sqrt(diag(covBeta));
    tVals = beta ./ se;
    pVals = 2 * (1 - tcdf(abs(tVals), df));
    R2 = 1 - SSE / SST;
end

function [lowDummy, highDummy] = distribution_to_dummies(distThis)
    n = length(distThis);
    lowDummy = nan(n, 1);
    highDummy = nan(n, 1);

    if isnumeric(distThis) || islogical(distThis)
        d = double(distThis(:));
        % Assumption from earlier code: 1=uniform, 2=low, 3=high.
        lowDummy = double(d == 2);
        highDummy = double(d == 3);
        unknown = ~(d == 1 | d == 2 | d == 3);
        lowDummy(unknown) = NaN;
        highDummy(unknown) = NaN;
        return
    end

    d = lower(string(distThis(:)));
    isUniform = d == "uniform" | d == "uni" | d == "u" | d == "gray" | d == "grey";
    isLow = d == "low" | d == "l" | d == "orange";
    isHigh = d == "high" | d == "h" | d == "green";

    lowDummy = double(isLow);
    highDummy = double(isHigh);

    unknown = ~(isUniform | isLow | isHigh);
    lowDummy(unknown) = NaN;
    highDummy(unknown) = NaN;
end

function patientSummary = make_patient_summary(allUnitSummary)
    patients = unique(allUnitSummary.Patient, 'stable');

    patientSummary = table();

    for i = 1:length(patients)
        p = patients(i);
        idx = allUnitSummary.Patient == p;
        nValid = sum(idx);
        nSig = sum(allUnitSummary.IsOutcomeSignificant(idx));

        tmp = table();
        tmp.Patient = p;
        tmp.NValidRightOFCUnits = nValid;
        tmp.NOutcomeSignificantUnits = nSig;
        tmp.PercentOutcomeSignificant = 100 * nSig / nValid;
        tmp.NWinPreferringSignificant = sum(allUnitSummary.IsOutcomeSignificant(idx) & allUnitSummary.OutcomeDelta_WinMinusLose_0to1s(idx) > 0);
        tmp.NLosePreferringSignificant = sum(allUnitSummary.IsOutcomeSignificant(idx) & allUnitSummary.OutcomeDelta_WinMinusLose_0to1s(idx) < 0);
        patientSummary = [patientSummary; tmp]; %#ok<AGROW>
    end
end

function leaveOneSummary = make_leave_one_patient_out_summary(allUnitSummary)
    patients = unique(allUnitSummary.Patient, 'stable');
    leaveOneSummary = table();

    for i = 1:length(patients)
        excludedPatient = patients(i);
        idx = allUnitSummary.Patient ~= excludedPatient;

        nValid = sum(idx);
        nSig = sum(allUnitSummary.IsOutcomeSignificant(idx));

        tmp = table();
        tmp.ExcludedPatient = excludedPatient;
        tmp.NValidUnitsRemaining = nValid;
        tmp.NOutcomeSignificantRemaining = nSig;
        tmp.PercentOutcomeSignificantRemaining = 100 * nSig / nValid;
        tmp.BinomialP_Against_5Percent = binomial_upper_tail(nSig, nValid, 0.05);

        leaveOneSummary = [leaveOneSummary; tmp]; %#ok<AGROW>
    end
end

function save_unit_outcome_figure(unitFigFolder, ptID, areaNameClean, chanNum, unitNum, ...
    timeCenters, winRasterX, winRasterY, loseRasterX, loseRasterY, ...
    validWinTrials, winMeanSmooth, loseMeanSmooth, winSEMSmooth, loseSEMSmooth, ...
    sigSegments, winColor, loseColor, windowBeforeSec, windowAfterSec, ...
    meanFRHz, outcomeDelta, aucWinLose, preferredOutcome)

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100 100 950 750]);

    subplot(2,1,1)
    hold on

    scatter(winRasterX, winRasterY, 6, 'filled', ...
        'MarkerFaceColor', winColor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.40)

    scatter(loseRasterX, loseRasterY + validWinTrials, 6, 'filled', ...
        'MarkerFaceColor', loseColor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.40)

    xline(0, '--k')
    yline(validWinTrials + 0.5, '--k')
    xlim([-windowBeforeSec windowAfterSec])
    xlabel('time from outcome onset (s)')
    ylabel('trials')
    title(sprintf('%s | %s | Ch %d Unit %d | Raster | FR %.2f Hz', ...
        ptID, areaNameClean, chanNum, unitNum, meanFRHz), 'Interpreter', 'none')
    legend({'win', 'lose'}, 'Location', 'best')
    box off

    subplot(2,1,2)
    hold on

    fill([timeCenters fliplr(timeCenters)], ...
        [winMeanSmooth + winSEMSmooth fliplr(winMeanSmooth - winSEMSmooth)], ...
        winColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

    fill([timeCenters fliplr(timeCenters)], ...
        [loseMeanSmooth + loseSEMSmooth fliplr(loseMeanSmooth - loseSEMSmooth)], ...
        loseColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

    plot(timeCenters, winMeanSmooth, 'Color', winColor, 'LineWidth', 1.8)
    plot(timeCenters, loseMeanSmooth, 'Color', loseColor, 'LineWidth', 1.8)

    xline(0, '--k')
    yline(0, ':k')

    add_sig_segments_to_axis(gca, sigSegments, {winMeanSmooth, loseMeanSmooth}, {winSEMSmooth, loseSEMSmooth});

    xlim([-windowBeforeSec windowAfterSec])
    xlabel('time from outcome onset (s)')
    ylabel('baseline z-scored firing rate')
    title(sprintf('PSTH: mean ± SEM | Delta 0-1s = %.3f | AUROC = %.3f | %s', ...
        outcomeDelta, aucWinLose, preferredOutcome), 'Interpreter', 'none')
    box off

    pdfName = sprintf('%s_%s_ch%d_unit%d_rightOFC_outcome.pdf', ptID, areaNameClean, chanNum, unitNum);
    exportgraphics(fig, fullfile(unitFigFolder, pdfName), 'ContentType', 'vector');
    close(fig);
end

function save_population_psth_figure(figFolder, timeCenters, popWinMean, popLoseMean, ...
    popSigSegments, winColor, loseColor, windowBeforeSec, windowAfterSec, nUnits)

    winPopMean = mean(popWinMean, 1, 'omitnan');
    losePopMean = mean(popLoseMean, 1, 'omitnan');
    winPopSEM = std(popWinMean, 0, 1, 'omitnan') ./ sqrt(size(popWinMean, 1));
    losePopSEM = std(popLoseMean, 0, 1, 'omitnan') ./ sqrt(size(popLoseMean, 1));

    winPopMeanSmooth = smoothdata(winPopMean, 'gaussian', 7);
    losePopMeanSmooth = smoothdata(losePopMean, 'gaussian', 7);
    winPopSEMSmooth = smoothdata(winPopSEM, 'gaussian', 7);
    losePopSEMSmooth = smoothdata(losePopSEM, 'gaussian', 7);

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100 100 900 550]);
    hold on

    fill([timeCenters fliplr(timeCenters)], ...
        [winPopMeanSmooth + winPopSEMSmooth fliplr(winPopMeanSmooth - winPopSEMSmooth)], ...
        winColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

    fill([timeCenters fliplr(timeCenters)], ...
        [losePopMeanSmooth + losePopSEMSmooth fliplr(losePopMeanSmooth - losePopSEMSmooth)], ...
        loseColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

    plot(timeCenters, winPopMeanSmooth, 'Color', winColor, 'LineWidth', 2.2)
    plot(timeCenters, losePopMeanSmooth, 'Color', loseColor, 'LineWidth', 2.2)

    xline(0, '--k')
    yline(0, ':k')
    add_sig_segments_to_axis(gca, popSigSegments, {winPopMeanSmooth, losePopMeanSmooth}, {winPopSEMSmooth, losePopSEMSmooth});

    xlim([-windowBeforeSec windowAfterSec])
    xlabel('time from outcome onset (s)')
    ylabel('baseline z-scored firing rate')
    title(sprintf('Right OFC population outcome PSTH | unit as observation | n = %d units', nUnits))
    legend({'win SEM', 'lose SEM', 'win', 'lose'}, 'Location', 'best', 'Box', 'off')
    box off

    exportgraphics(fig, fullfile(figFolder, 'rightOFC_population_outcome_PSTH.pdf'), 'ContentType', 'vector');
    close(fig);
end

function save_outcome_delta_figure(figFolder, allUnitSummary)
    delta = allUnitSummary.OutcomeDelta_WinMinusLose_0to1s;
    isSig = allUnitSummary.IsOutcomeSignificant;

    [deltaSorted, sortIdx] = sort(delta, 'ascend');
    isSigSorted = isSig(sortIdx);

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100 100 850 520]);
    hold on

    x = 1:length(deltaSorted);
    scatter(x(~isSigSorted), deltaSorted(~isSigSorted), 45, 'o', 'LineWidth', 1.2)
    scatter(x(isSigSorted), deltaSorted(isSigSorted), 65, 'filled')
    yline(0, '--k')

    xlabel('right OFC units sorted by outcome modulation')
    ylabel('mean zFR(win) - mean zFR(lose), 0-1 s')
    title('Right OFC outcome modulation per unit')
    legend({'not significant', 'cluster-corrected significant'}, 'Location', 'best', 'Box', 'off')
    box off

    exportgraphics(fig, fullfile(figFolder, 'rightOFC_outcome_modulation_delta.pdf'), 'ContentType', 'vector');
    close(fig);
end

function save_auc_figure(figFolder, allUnitSummary)
    auc = allUnitSummary.AUROC_WinVsLose_0to1s;
    isSig = allUnitSummary.IsOutcomeSignificant;

    [aucSorted, sortIdx] = sort(auc, 'ascend');
    isSigSorted = isSig(sortIdx);

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100 100 850 520]);
    hold on

    x = 1:length(aucSorted);
    scatter(x(~isSigSorted), aucSorted(~isSigSorted), 45, 'o', 'LineWidth', 1.2)
    scatter(x(isSigSorted), aucSorted(isSigSorted), 65, 'filled')
    yline(0.5, '--k')

    xlabel('right OFC units sorted by AUROC')
    ylabel('AUROC for win vs lose, 0-1 s')
    title('Single-unit win/lose discriminability in right OFC')
    legend({'not significant', 'cluster-corrected significant'}, 'Location', 'best', 'Box', 'off')
    box off

    exportgraphics(fig, fullfile(figFolder, 'rightOFC_AUROC_win_vs_lose.pdf'), 'ContentType', 'vector');
    close(fig);
end

function save_patient_summary_figure(patientSummary, figFolder)
    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100 100 850 520]);

    bar(patientSummary.PercentOutcomeSignificant)
    hold on
    ylabel('% significant outcome units')
    xlabel('participant')
    title('Right OFC outcome encoding by participant')
    xticks(1:height(patientSummary))
    xticklabels(patientSummary.Patient)
    xtickangle(35)
    ylim([0 max(100, max(patientSummary.PercentOutcomeSignificant, [], 'omitnan') + 10)])

    for i = 1:height(patientSummary)
        text(i, patientSummary.PercentOutcomeSignificant(i) + 2, ...
            sprintf('%d/%d', patientSummary.NOutcomeSignificantUnits(i), patientSummary.NValidRightOFCUnits(i)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold')
    end
    box off

    exportgraphics(fig, fullfile(figFolder, 'rightOFC_outcome_patient_summary.pdf'), 'ContentType', 'vector');
    close(fig);
end

function save_leave_one_out_figure(leaveOneSummary, figFolder)
    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100 100 900 520]);

    bar(leaveOneSummary.PercentOutcomeSignificantRemaining)
    hold on
    ylabel('% significant outcome units remaining')
    xlabel('excluded participant')
    title('Leave-one-participant-out robustness: right OFC outcome encoding')
    xticks(1:height(leaveOneSummary))
    xticklabels(leaveOneSummary.ExcludedPatient)
    xtickangle(35)
    ylim([0 max(100, max(leaveOneSummary.PercentOutcomeSignificantRemaining, [], 'omitnan') + 10)])

    for i = 1:height(leaveOneSummary)
        text(i, leaveOneSummary.PercentOutcomeSignificantRemaining(i) + 2, ...
            sprintf('%d/%d', leaveOneSummary.NOutcomeSignificantRemaining(i), leaveOneSummary.NValidUnitsRemaining(i)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold')
    end
    box off

    exportgraphics(fig, fullfile(figFolder, 'rightOFC_outcome_leave_one_patient_out.pdf'), 'ContentType', 'vector');
    close(fig);
end

function save_regression_figure(figFolder, allRegressionSummary)
    coef = allRegressionSummary.OutcomeCoef;
    q = allRegressionSummary.OutcomeQ_BHFDR;
    isSig = q < 0.05;

    [coefSorted, sortIdx] = sort(coef, 'ascend');
    isSigSorted = isSig(sortIdx);

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100 100 850 520]);
    hold on

    x = 1:length(coefSorted);
    scatter(x(~isSigSorted), coefSorted(~isSigSorted), 45, 'o', 'LineWidth', 1.2)
    scatter(x(isSigSorted), coefSorted(isSigSorted), 65, 'filled')
    yline(0, '--k')

    xlabel('right OFC units sorted by outcome regression coefficient')
    ylabel('outcome coefficient: win vs lose')
    title('Trial-level regression controlling for choice, distribution, card value, and trial index')
    legend({'not FDR significant', 'FDR significant'}, 'Location', 'best', 'Box', 'off')
    box off

    exportgraphics(fig, fullfile(figFolder, 'rightOFC_outcome_regression_coefficients.pdf'), 'ContentType', 'vector');
    close(fig);
end

function add_sig_segments_to_axis(ax, sigSegments, meanCurves, semCurves)
    upperVals = [];
    lowerVals = [];

    for gg = 1:length(meanCurves)
        upperVals = [upperVals, meanCurves{gg} + semCurves{gg}]; %#ok<AGROW>
        lowerVals = [lowerVals, meanCurves{gg} - semCurves{gg}]; %#ok<AGROW>
    end

    curveHigh = max(upperVals, [], 'omitnan');
    curveLow = min(lowerVals, [], 'omitnan');

    if isempty(curveHigh) || isnan(curveHigh); curveHigh = 1; end
    if isempty(curveLow) || isnan(curveLow); curveLow = -1; end

    yRange = curveHigh - curveLow;
    if yRange == 0 || isnan(yRange); yRange = 1; end

    if isempty(sigSegments)
        ylim(ax, [curveLow - 0.15*yRange, curveHigh + 0.25*yRange]);
        return
    end

    for ss = 1:size(sigSegments, 1)
        sigStart = sigSegments(ss, 1);
        sigEnd = sigSegments(ss, 2);
        sigP = sigSegments(ss, 3);
        sigY = curveHigh + (0.15 + 0.12*(ss-1)) * yRange;

        plot(ax, [sigStart sigEnd], [sigY sigY], '-', ...
            'Color', [0.5 0.5 0.5], 'LineWidth', 4, 'HandleVisibility', 'off')

        text(ax, mean([sigStart sigEnd]), sigY + 0.04*yRange, ...
            sprintf('p=%.3f', sigP), ...
            'Color', [0.35 0.35 0.35], 'FontSize', 8, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    end

    ylim(ax, [curveLow - 0.15*yRange, curveHigh + (0.30 + 0.12*size(sigSegments,1))*yRange]);
end

function auc = compute_auc(scores, labels)
    scores = scores(:);
    labels = labels(:);
    validIdx = ~isnan(scores) & ~isnan(labels);
    scores = scores(validIdx);
    labels = labels(validIdx);

    pos = labels == 1;
    neg = labels == 0;
    nPos = sum(pos);
    nNeg = sum(neg);

    if nPos < 1 || nNeg < 1
        auc = NaN;
        return
    end

    ranks = local_tiedrank(scores);
    auc = (sum(ranks(pos)) - nPos * (nPos + 1) / 2) / (nPos * nNeg);
end

function ranks = local_tiedrank(x)
    x = x(:);
    [xs, sortIdx] = sort(x);
    ranksSorted = nan(size(xs));

    i = 1;
    n = length(xs);
    while i <= n
        j = i;
        while j < n && xs(j+1) == xs(i)
            j = j + 1;
        end
        ranksSorted(i:j) = mean(i:j);
        i = j + 1;
    end

    ranks = nan(size(x));
    ranks(sortIdx) = ranksSorted;
end

function p = binomial_upper_tail(kObs, n, p0)
    if n <= 0 || kObs < 0
        p = NaN;
        return
    end

    ks = kObs:n;
    logTerms = gammaln(n + 1) - gammaln(ks + 1) - gammaln(n - ks + 1) + ...
        ks .* log(p0) + (n - ks) .* log(1 - p0);
    p = sum(exp(logTerms));
    p = min(max(p, 0), 1);
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

function z = zscore_safe(x)
    x = double(x(:));
    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');
    if sd == 0 || isnan(sd)
        z = zeros(size(x));
    else
        z = (x - mu) ./ sd;
    end
end

function tf = has_bhv_field(bhvData, fieldName)
    if istable(bhvData)
        tf = ismember(fieldName, bhvData.Properties.VariableNames);
    elseif isstruct(bhvData)
        tf = isfield(bhvData, fieldName);
    else
        tf = false;
    end
end

function x = get_bhv_field(bhvData, fieldName)
    if istable(bhvData)
        if ~ismember(fieldName, bhvData.Properties.VariableNames)
            error('Could not find bhvData.%s. Available fields are: %s', ...
                fieldName, strjoin(bhvData.Properties.VariableNames, ', '));
        end
        x = bhvData.(fieldName);
    elseif isstruct(bhvData)
        if ~isfield(bhvData, fieldName)
            error('Could not find bhvData.%s. Available fields are: %s', ...
                fieldName, strjoin(fieldnames(bhvData), ', '));
        end
        x = bhvData.(fieldName);
    else
        error('bhvData must be table or struct. Current class: %s', class(bhvData));
    end

    if isrow(x) && ~ischar(x) && ~isstring(x)
        x = x(:);
    end
end

function x = get_numeric_field(bhvData, fieldName)
    xRaw = get_bhv_field(bhvData, fieldName);
    if iscell(xRaw)
        xRaw = string(xRaw);
    end
    x = double(xRaw);
    x = x(:);
end

function x = get_string_field(bhvData, fieldName)
    xRaw = get_bhv_field(bhvData, fieldName);
    x = lower(string(xRaw));
    x = strtrim(x(:));
end
