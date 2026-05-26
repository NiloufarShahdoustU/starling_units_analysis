clc;
clear;
close all;

% ========================================================================
% Starling single-neuron 4-stage encoding analysis
%
% Stages:
%   1) Card backs appear  -> distribution / card color: uniform, low, high
%   2) Card is revealed   -> card number group: 1-3, 4-6, 7-9
%   3) Decision stage     -> choice: arrowup vs arrowdown, tested pre-choice
%   4) Outcome stage      -> outcome: win vs lose, tested post-outcome
%
% What this script saves:
%   - one compact 4-stage PSTH figure per significant unit
%   - stage_encoding_summary.csv
%   - stage_cluster_stats.csv
%   - stage_window_stats.csv
%   - summary_stage_counts.csv / summary_stage_counts.pdf
%   - summary_area_stage_percent.csv / summary_area_stage_percent.pdf
%   - summary_patient_stage_percent.csv / summary_patient_stage_percent.pdf
%
% Removed on purpose:
%   - summary_stage_heatmap.pdf
%   - summary_area_stage_counts.pdf
%   - summary_stage_overlap.pdf
%
% Notes:
%   - FR is baseline z-scored using -1 to -0.25 sec before each event.
%   - Time-resolved statistics use cluster-based permutation correction.
%   - Cluster-forming threshold is p < 0.05; corrected cluster p < 0.05.
%   - For choice, the test window is -0.75 to 0 sec, so it focuses on
%     pre-choice/decision activity instead of post-choice feedback.
% ========================================================================

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};

inputFolder = '\\155.100.91.44\d\Data\Nill\starling\spikes\spike_data\';
OutputFolder = '\\155.100.91.44\d\Code\Nill\Starling_units_analysis\starling_5_four_stage_encoding\';

if ~exist(OutputFolder, 'dir')
    mkdir(OutputFolder);
end

unitFigureFolder = fullfile(OutputFolder, 'Unit_4Stage_Figures');
if ~exist(unitFigureFolder, 'dir')
    mkdir(unitFigureFolder);
end

% Delete old summary visualizations that are no longer used.
oldSummaryFigures = { ...
    'summary_stage_heatmap.pdf', ...
    'summary_area_stage_counts.pdf', ...
    'summary_stage_overlap.pdf'};
for ff = 1:numel(oldSummaryFigures)
    oldPath = fullfile(OutputFolder, oldSummaryFigures{ff});
    if exist(oldPath, 'file')
        delete(oldPath);
    end
end

% Delete old overlap CSVs too, because the overlap summary is no longer saved.
oldSummaryCSVs = { ...
    'summary_stage_overlap.csv', ...
    'summary_unit_stage_overlap.csv'};
for ff = 1:numel(oldSummaryCSVs)
    oldPath = fullfile(OutputFolder, oldSummaryCSVs{ff});
    if exist(oldPath, 'file')
        delete(oldPath);
    end
end

% ============================= settings =================================
windowBeforeSec = 1;
windowAfterSec  = 1;
binSizeSec = 0.05;
smoothBins = 7;

baselineWindowSec = [-1 -0.25];

minMeanFRHz = 0.5;
minFracTrialsWithSpikes = 0.10;
minTrialsPerGroup = 5;

permWindowMs = 200;
permStrideMs = 50;
nPerm = 1000;
alphaPerm = 0.05;
clusterFormingAlpha = 0.05;

saveOnlySignificantUnitFigures = true;
% If this makes too many files, keep true. If you want every unit figure,
% set saveOnlySignificantUnitFigures = false.

rng(1); % reproducible permutations
% ========================================================================

permWindowSec = permWindowMs / 1000;
permStrideSec = permStrideMs / 1000;
analysisWindowSec = windowBeforeSec + windowAfterSec;

timeEdges = -windowBeforeSec:binSizeSec:windowAfterSec;
timeCenters = timeEdges(1:end-1) + binSizeSec/2;
baselineIdx = timeCenters >= baselineWindowSec(1) & timeCenters < baselineWindowSec(2);

% Colors
uniformColor = [0.50 0.50 0.50];
lowColor     = [0.85 0.35 0.00];
highColor    = [0.25 0.65 0.25];
numLowColor  = [0.20 0.45 0.85];
numMidColor  = [0.50 0.50 0.50];
numHighColor = [0.85 0.35 0.00];
upColor      = [0.00 0.45 0.85];
downColor    = [0.85 0.35 0.00];
winColor     = [0.00 0.60 0.00];
loseColor    = [0.85 0.00 0.00];

allUnitSummary = table();
allClusterStats = table();
allWindowStats = table();

for pt = 1:length(microPts)

    ptID = microPts{pt};

    fprintf('Processing patient %s\n', ptID);


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

    stages = build_stage_definitions(eventTime, bhvData, SampleRes, ...
        uniformColor, lowColor, highColor, ...
        numLowColor, numMidColor, numHighColor, ...
        upColor, downColor, winColor, loseColor);

    allChanUnits = unique(ChanUnitTimestamp(:,1:2), 'rows');
    allChanUnits = double(allChanUnits);

    % Unit 255 means no waveform was saved.
    allChanUnits(allChanUnits(:,2) == 255, :) = [];

    % Remove channels outside included channels.
    maxInclChan = max(double(inclChans(:)));
    allChanUnits(allChanUnits(:,1) > maxInclChan, :) = [];

    fprintf('Found %d candidate channel-units.\n', size(allChanUnits, 1));

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

        unitResults = cell(length(stages), 1);
        unitHasAnySignificantStage = false;

        for ss = 1:length(stages)

            stage = stages(ss);

            result = analyze_one_stage_for_unit( ...
                stage, unitSpikeTimes, SampleRes, ...
                windowBeforeSec, windowAfterSec, binSizeSec, ...
                timeEdges, timeCenters, baselineIdx, analysisWindowSec, ...
                minMeanFRHz, minFracTrialsWithSpikes, minTrialsPerGroup, ...
                permWindowSec, permStrideSec, nPerm, alphaPerm, clusterFormingAlpha, ...
                smoothBins);

            unitResults{ss} = result;

            if ~result.valid
                continue
            end

            unitHasAnySignificantStage = unitHasAnySignificantStage || result.isSignificant;

            % ---------------------- unit/stage summary -------------------
            tmpSummary = table();
            tmpSummary.Patient = string(ptID);
            tmpSummary.Area = string(areaNameClean);
            tmpSummary.Channel = chanNum;
            tmpSummary.Unit = unitNum;
            tmpSummary.UnitID = string(sprintf('%s_%s_ch%d_u%d', ptID, areaNameClean, chanNum, unitNum));
            tmpSummary.Stage = string(stage.name);
            tmpSummary.StageLabel = string(stage.longName);
            tmpSummary.EventFieldUsed = string(stage.eventFieldUsed);
            tmpSummary.TestWindowStartSec = stage.testWindowSec(1);
            tmpSummary.TestWindowEndSec = stage.testWindowSec(2);
            tmpSummary.MeanFRHz = result.meanFRHz;
            tmpSummary.FracTrialsWithSpikes = result.fracTrialsWithSpikes;
            tmpSummary.NGroups = length(stage.groupNames);
            tmpSummary.NTrialsTotal = result.nTrialsTotal;
            tmpSummary.GroupTrialCounts = string(strjoin(arrayfun(@num2str, result.validTrials(:)', 'UniformOutput', false), ','));
            tmpSummary.IsSignificant = result.isSignificant;
            tmpSummary.NSignificantClusters = result.nSignificantClusters;
            tmpSummary.MinClusterP = result.minClusterP;
            tmpSummary.MaxObservedF = max(result.observedF, [], 'omitnan');
            tmpSummary.EffectSizeMaxGroupDiff = result.effectSizeMaxGroupDiff;

            allUnitSummary = [allUnitSummary; tmpSummary]; %#ok<AGROW>

            % ---------------------- cluster stats ------------------------
            if ~isempty(result.clusterStats)
                tmpClusterStats = result.clusterStats;
                nRows = height(tmpClusterStats);
                tmpClusterStats.Patient = repmat(string(ptID), nRows, 1);
                tmpClusterStats.Area = repmat(string(areaNameClean), nRows, 1);
                tmpClusterStats.Channel = repmat(chanNum, nRows, 1);
                tmpClusterStats.Unit = repmat(unitNum, nRows, 1);
                tmpClusterStats.UnitID = repmat(string(sprintf('%s_%s_ch%d_u%d', ptID, areaNameClean, chanNum, unitNum)), nRows, 1);
                tmpClusterStats.Stage = repmat(string(stage.name), nRows, 1);
                tmpClusterStats.EventFieldUsed = repmat(string(stage.eventFieldUsed), nRows, 1);
                tmpClusterStats = movevars(tmpClusterStats, ...
                    {'Patient', 'Area', 'Channel', 'Unit', 'UnitID', 'Stage', 'EventFieldUsed'}, ...
                    'Before', 1);
                allClusterStats = [allClusterStats; tmpClusterStats]; %#ok<AGROW>
            end

            % ---------------------- window stats -------------------------
            if ~isempty(result.windowStats)
                tmpWindowStats = result.windowStats;
                nRows = height(tmpWindowStats);
                tmpWindowStats.Patient = repmat(string(ptID), nRows, 1);
                tmpWindowStats.Area = repmat(string(areaNameClean), nRows, 1);
                tmpWindowStats.Channel = repmat(chanNum, nRows, 1);
                tmpWindowStats.Unit = repmat(unitNum, nRows, 1);
                tmpWindowStats.UnitID = repmat(string(sprintf('%s_%s_ch%d_u%d', ptID, areaNameClean, chanNum, unitNum)), nRows, 1);
                tmpWindowStats.Stage = repmat(string(stage.name), nRows, 1);
                tmpWindowStats.EventFieldUsed = repmat(string(stage.eventFieldUsed), nRows, 1);
                tmpWindowStats = movevars(tmpWindowStats, ...
                    {'Patient', 'Area', 'Channel', 'Unit', 'UnitID', 'Stage', 'EventFieldUsed'}, ...
                    'Before', 1);
                allWindowStats = [allWindowStats; tmpWindowStats]; %#ok<AGROW>
            end
        end

        if ~saveOnlySignificantUnitFigures || unitHasAnySignificantStage
            save_4stage_unit_figure(unitResults, stages, timeCenters, ...
                ptID, areaNameClean, chanNum, unitNum, unitFigureFolder);
        end

        if unitHasAnySignificantStage
            fprintf('Significant unit: %s | %s | ch %d unit %d\n', ...
                ptID, areaNameClean, chanNum, unitNum);
        end
    end
end

% =============================== save tables =============================
% Add an optional second-level FDR correction across all valid unit-stage tests.
% The main IsSignificant column is already corrected over time within each
% unit/stage by cluster-based permutation. This global FDR column is more
% conservative and is useful for summary reporting across many units.
if ~isempty(allUnitSummary)
    allUnitSummary.GlobalFDR_Q = bh_fdr(allUnitSummary.MinClusterP);
    allUnitSummary.IsSignificant_GlobalFDR = allUnitSummary.GlobalFDR_Q < alphaPerm;
end

summaryPath = fullfile(OutputFolder, 'stage_encoding_summary.csv');
writetable(allUnitSummary, summaryPath);
fprintf('\nSaved unit/stage summary:\n%s\n', summaryPath);

clusterPath = fullfile(OutputFolder, 'stage_cluster_stats.csv');
writetable(allClusterStats, clusterPath);
fprintf('Saved cluster stats:\n%s\n', clusterPath);

windowPath = fullfile(OutputFolder, 'stage_window_stats.csv');
writetable(allWindowStats, windowPath);
fprintf('Saved window stats:\n%s\n', windowPath);

% ============================ summary figures ============================
if ~isempty(allUnitSummary)
    save_summary_figures(allUnitSummary, OutputFolder);
end

fprintf('\nDONE. Results saved in:\n%s\n', OutputFolder);

% ========================================================================
% Local functions
% ========================================================================

function stages = build_stage_definitions(eventTime, bhvData, SampleRes, ...
    uniformColor, lowColor, highColor, ...
    numLowColor, numMidColor, numHighColor, ...
    upColor, downColor, winColor, loseColor)

    stages = struct([]);

    % --------------------------------------------------------------------
    % Stage 1: card backs / card color / deck identity
    % --------------------------------------------------------------------
    [cardShowTime, cardShowField] = get_event_time(eventTime, {'cardShowTime'}, 'CardColor');
    distribution = get_distribution_vector(bhvData);
    nTrials = length(distribution);

    uniformTrials = find(get_distribution_mask(distribution, 'uniform', 1, nTrials));
    lowTrials     = find(get_distribution_mask(distribution, 'low',     2, nTrials));
    highTrials    = find(get_distribution_mask(distribution, 'high',    3, nTrials));

    stages(1).name = 'CardColor';
    stages(1).longName = 'Backs of two cards / color-deck encoding';
    stages(1).alignTime = cardShowTime;
    stages(1).eventFieldUsed = cardShowField;
    stages(1).groupNames = {'uniform', 'low', 'high'};
    stages(1).groupTrials = {uniformTrials, lowTrials, highTrials};
    stages(1).groupColors = [uniformColor; lowColor; highColor];
    stages(1).testWindowSec = [0 1];
    stages(1).xLabel = 'time from two card backs appearing / color cue onset (s)';

    % --------------------------------------------------------------------
    % Stage 2: number processing after participant's card is revealed
    % There is no separate flip-time field in bhvData.
    % Your bhvData has spaceRT, so flip/reveal time is computed as:
    %     flipTime = cardShowTime + spaceRT
    % spaceRT is auto-detected as seconds or milliseconds.
    % --------------------------------------------------------------------
    spaceRT = get_numeric_field(bhvData, 'spaceRT');
    flipTime = add_rt_to_event_time(cardShowTime, spaceRT, SampleRes);
    flipField = 'cardShowTime_plus_bhvData.spaceRT';

    myCard = get_numeric_field(bhvData, 'myCard');
    myCard = myCard(:);

    numberLowTrials  = find(myCard >= 1 & myCard <= 3);
    numberMidTrials  = find(myCard >= 4 & myCard <= 6);
    numberHighTrials = find(myCard >= 7 & myCard <= 9);

    stages(2).name = 'CardNumber';
    stages(2).longName = 'Card reveal / number encoding';
    stages(2).alignTime = flipTime;
    stages(2).eventFieldUsed = flipField;
    stages(2).groupNames = {'card 1-3', 'card 4-6', 'card 7-9'};
    stages(2).groupTrials = {numberLowTrials, numberMidTrials, numberHighTrials};
    stages(2).groupColors = [numLowColor; numMidColor; numHighColor];
    stages(2).testWindowSec = [0 1];
    stages(2).xLabel = 'time from card reveal / flip event (s)';

    % --------------------------------------------------------------------
    % Stage 3: decision stage, tested before choice time
    % Use the exact event timestamp if choiceAndFeedbackTime exists.
    % Otherwise compute from behavioral RTs:
    %     choiceTime = cardShowTime + spaceRT + arrowRT
    % arrowRT is auto-detected as seconds or milliseconds.
    % --------------------------------------------------------------------
    arrowRT = get_numeric_field(bhvData, 'arrowRT');

    if isfield(eventTime, 'choiceAndFeedbackTime')
        choiceTime = double(eventTime.choiceAndFeedbackTime);
        choiceField = 'choiceAndFeedbackTime';
    else
        choiceTime = add_rt_to_event_time(flipTime, arrowRT, SampleRes);
        choiceField = 'cardShowTime_plus_spaceRT_plus_arrowRT';
    end

    choice = get_string_field(bhvData, 'choice');
    arrowUpTrials = find(choice == "arrowup" | choice == "up" | choice == "higher" | choice == "high");
    arrowDownTrials = find(choice == "arrowdown" | choice == "down" | choice == "lower" | choice == "low");

    stages(3).name = 'Choice';
    stages(3).longName = 'Decision / choice-direction encoding';
    stages(3).alignTime = choiceTime;
    stages(3).eventFieldUsed = choiceField;
    stages(3).groupNames = {'arrow up', 'arrow down'};
    stages(3).groupTrials = {arrowUpTrials, arrowDownTrials};
    stages(3).groupColors = [upColor; downColor];
    stages(3).testWindowSec = [-0.75 0];
    stages(3).xLabel = 'time from choice / response (s)';

    % --------------------------------------------------------------------
    % Stage 4: outcome feedback stage
    % In this task, feedback/outcome appears after the arrow choice.
    % Use choiceAndFeedbackTime if present; otherwise use the same derived
    % time as the decision response: cardShowTime + spaceRT + arrowRT.
    % --------------------------------------------------------------------
    if isfield(eventTime, 'choiceAndFeedbackTime')
        outcomeTime = double(eventTime.choiceAndFeedbackTime);
        outcomeField = 'choiceAndFeedbackTime';
    else
        outcomeTime = choiceTime;
        outcomeField = choiceField;
    end

    outcome = get_string_field(bhvData, 'outcome');
    winTrials = find(outcome == "win" | outcome == "correct" | outcome == "reward" | outcome == "1");
    loseTrials = find(outcome == "lose" | outcome == "loss" | outcome == "incorrect" | outcome == "0");

    stages(4).name = 'Outcome';
    stages(4).longName = 'Outcome / feedback encoding';
    stages(4).alignTime = outcomeTime;
    stages(4).eventFieldUsed = outcomeField;
    stages(4).groupNames = {'win', 'lose'};
    stages(4).groupTrials = {winTrials, loseTrials};
    stages(4).groupColors = [winColor; loseColor];
    stages(4).testWindowSec = [0 1];
    stages(4).xLabel = 'time from outcome / feedback onset (s)';

end

function result = analyze_one_stage_for_unit(stage, unitSpikeTimes, SampleRes, ...
    windowBeforeSec, windowAfterSec, binSizeSec, ...
    timeEdges, timeCenters, baselineIdx, analysisWindowSec, ...
    minMeanFRHz, minFracTrialsWithSpikes, minTrialsPerGroup, ...
    permWindowSec, permStrideSec, nPerm, alphaPerm, clusterFormingAlpha, ...
    smoothBins)

    result = initialize_empty_result();

    if isempty(stage.alignTime) || all(isnan(stage.alignTime))
        return
    end

    nGroups = length(stage.groupNames);

    rasterX = cell(nGroups, 1);
    rasterY = cell(nGroups, 1);
    trialFR = cell(nGroups, 1);
    trialSpikeCounts = cell(nGroups, 1);
    validTrials = zeros(nGroups, 1);

    for gg = 1:nGroups
        [rasterX{gg}, rasterY{gg}, trialFR{gg}, trialSpikeCounts{gg}, validTrials(gg)] = ...
            get_trial_fr(stage.groupTrials{gg}, stage.alignTime, unitSpikeTimes, SampleRes, ...
            windowBeforeSec, windowAfterSec, timeEdges, binSizeSec);
    end

    if any(validTrials < minTrialsPerGroup)
        return
    end

    allTrialSpikeCounts = vertcat(trialSpikeCounts{:});

    if isempty(allTrialSpikeCounts)
        return
    end

    nTrialsTotal = length(allTrialSpikeCounts);
    totalSpikesInWindow = sum(allTrialSpikeCounts);
    meanFRHz = totalSpikesInWindow / (nTrialsTotal * analysisWindowSec);
    fracTrialsWithSpikes = sum(allTrialSpikeCounts > 0) / nTrialsTotal;

    if meanFRHz < minMeanFRHz || fracTrialsWithSpikes < minFracTrialsWithSpikes
        return
    end

    allTrialFR = vertcat(trialFR{:});

    baselineVals = allTrialFR(:, baselineIdx);
    baselineMean = mean(baselineVals(:), 'omitnan');
    baselineStd  = std(baselineVals(:), 'omitnan');

    if baselineStd == 0 || isnan(baselineStd)
        baselineStd = 1;
    end

    trialZ = cell(nGroups, 1);
    meanZ = cell(nGroups, 1);
    semZ = cell(nGroups, 1);
    meanSmooth = cell(nGroups, 1);
    semSmooth = cell(nGroups, 1);

    for gg = 1:nGroups
        trialZ{gg} = (trialFR{gg} - baselineMean) ./ baselineStd;
        meanZ{gg} = mean(trialZ{gg}, 1, 'omitnan');
        semZ{gg} = std(trialZ{gg}, 0, 1, 'omitnan') ./ sqrt(size(trialZ{gg}, 1));
        meanSmooth{gg} = smoothdata(meanZ{gg}, 'gaussian', smoothBins);
        semSmooth{gg} = smoothdata(semZ{gg}, 'gaussian', smoothBins);
    end

    [windowStats, clusterStats, sigSegments, sigLabels, observedF] = ...
        run_cluster_permutation_anova( ...
        trialZ, stage.groupNames, timeCenters, stage.testWindowSec, ...
        permWindowSec, permStrideSec, nPerm, alphaPerm, clusterFormingAlpha);

    testIdx = timeCenters >= stage.testWindowSec(1) & timeCenters <= stage.testWindowSec(2);
    groupMeansInTestWindow = nan(nGroups, 1);
    for gg = 1:nGroups
        groupMeansInTestWindow(gg) = mean(trialZ{gg}(:, testIdx), 'all', 'omitnan');
    end

    pairwiseDiffs = [];
    for g1 = 1:nGroups-1
        for g2 = g1+1:nGroups
            pairwiseDiffs = [pairwiseDiffs; abs(groupMeansInTestWindow(g1) - groupMeansInTestWindow(g2))]; %#ok<AGROW>
        end
    end

    if isempty(pairwiseDiffs)
        effectSizeMaxGroupDiff = NaN;
    else
        effectSizeMaxGroupDiff = max(pairwiseDiffs, [], 'omitnan');
    end

    isSignificant = ~isempty(sigSegments);

    if isempty(clusterStats)
        minClusterP = NaN;
        nSignificantClusters = 0;
    else
        minClusterP = min(clusterStats.ClusterP, [], 'omitnan');
        nSignificantClusters = sum(clusterStats.ClusterSignificant);
    end

    result.valid = true;
    result.rasterX = rasterX;
    result.rasterY = rasterY;
    result.trialZ = trialZ;
    result.meanZ = meanZ;
    result.semZ = semZ;
    result.meanSmooth = meanSmooth;
    result.semSmooth = semSmooth;
    result.validTrials = validTrials;
    result.nTrialsTotal = nTrialsTotal;
    result.meanFRHz = meanFRHz;
    result.fracTrialsWithSpikes = fracTrialsWithSpikes;
    result.windowStats = windowStats;
    result.clusterStats = clusterStats;
    result.sigSegments = sigSegments;
    result.sigLabels = sigLabels;
    result.observedF = observedF;
    result.isSignificant = isSignificant;
    result.nSignificantClusters = nSignificantClusters;
    result.minClusterP = minClusterP;
    result.effectSizeMaxGroupDiff = effectSizeMaxGroupDiff;
end

function result = initialize_empty_result()
    result.valid = false;
    result.rasterX = {};
    result.rasterY = {};
    result.trialZ = {};
    result.meanZ = {};
    result.semZ = {};
    result.meanSmooth = {};
    result.semSmooth = {};
    result.validTrials = [];
    result.nTrialsTotal = 0;
    result.meanFRHz = NaN;
    result.fracTrialsWithSpikes = NaN;
    result.windowStats = table();
    result.clusterStats = table();
    result.sigSegments = [];
    result.sigLabels = {};
    result.observedF = NaN;
    result.isSignificant = false;
    result.nSignificantClusters = 0;
    result.minClusterP = NaN;
    result.effectSizeMaxGroupDiff = NaN;
end

function [rasterX, rasterY, trialFR, trialSpikeCounts, validTrials] = get_trial_fr( ...
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

        windowStart = alignTime(tr) - windowBeforeSec * SampleRes;
        windowEnd   = alignTime(tr) + windowAfterSec  * SampleRes;

        spikesInWindow = unitSpikeTimes( ...
            unitSpikeTimes >= windowStart & ...
            unitSpikeTimes <= windowEnd);

        relSpikesSec = (spikesInWindow - alignTime(tr)) ./ SampleRes;

        rasterX = [rasterX; relSpikesSec(:)]; %#ok<AGROW>
        rasterY = [rasterY; validTrials .* ones(length(relSpikesSec), 1)]; %#ok<AGROW>

        counts = histcounts(relSpikesSec, timeEdges);
        trialFR(validTrials, :) = counts ./ binSizeSec; %#ok<AGROW>
        trialSpikeCounts(validTrials, 1) = sum(counts); %#ok<AGROW>
    end
end

function [windowStats, clusterStats, sigSegments, sigLabels, observedF] = ...
    run_cluster_permutation_anova(trialZ, groupNames, timeCenters, testWindowSec, ...
    permWindowSec, permStrideSec, nPerm, alphaPerm, clusterFormingAlpha)

    nGroups = length(trialZ);
    sigSegments = [];
    sigLabels = {};
    clusterStats = table();
    windowStats = table();

    if nGroups < 2
        observedF = NaN;
        return
    end

    for gg = 1:nGroups
        if isempty(trialZ{gg}) || size(trialZ{gg}, 1) < 2
            observedF = NaN;
            return
        end
    end

    permStartTimes = testWindowSec(1):permStrideSec:(testWindowSec(2) - permWindowSec);
    permEndTimes = permStartTimes + permWindowSec;
    nWindows = length(permStartTimes);

    if nWindows < 1
        observedF = NaN;
        return
    end

    windowMat = cell(nGroups, 1);
    for gg = 1:nGroups
        windowMat{gg} = nan(size(trialZ{gg}, 1), nWindows);
    end

    for ww = 1:nWindows
        thisIdx = timeCenters >= permStartTimes(ww) & timeCenters < permEndTimes(ww);
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
        [y, g] = make_y_g_from_window_mat(windowMat, ww);
        [observedP(ww), observedF(ww)] = local_anova_p_f(y, g);
    end

    observedSig = observedP < clusterFormingAlpha;
    [obsStartIdx, obsEndIdx, obsMass] = find_clusters(observedSig, observedF);

    allWindowMat = vertcat(windowMat{:});
    groupLabels = [];
    for gg = 1:nGroups
        groupLabels = [groupLabels; gg .* ones(size(windowMat{gg}, 1), 1)]; %#ok<AGROW>
    end

    nTotal = length(groupLabels);
    maxPermClusterMass = zeros(nPerm, 1);

    for pp = 1:nPerm
        shuffledLabels = groupLabels(randperm(nTotal));
        permF = nan(1, nWindows);
        permP = nan(1, nWindows);

        for ww = 1:nWindows
            y = allWindowMat(:, ww);
            g = shuffledLabels;
            validIdx = ~isnan(y) & ~isnan(g);
            y = y(validIdx);
            g = g(validIdx);
            [permP(ww), permF(ww)] = local_anova_p_f(y, g);
        end

        permSig = permP < clusterFormingAlpha;
        [~, ~, permMass] = find_clusters(permSig, permF);

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
    pairwiseLabels = strings(nObsClusters, 1);

    for cc = 1:nObsClusters
        clusterP(cc) = (1 + sum(maxPermClusterMass >= obsMass(cc))) / (nPerm + 1);
        clusterSignificant(cc) = clusterP(cc) < alphaPerm;
        clusterStartSec(cc) = permStartTimes(obsStartIdx(cc));
        clusterEndSec(cc) = permEndTimes(obsEndIdx(cc));

        if clusterSignificant(cc)
            [pairwiseLabels(cc), ~] = pairwise_cluster_labels( ...
                trialZ, groupNames, timeCenters, [clusterStartSec(cc) clusterEndSec(cc)], nPerm, alphaPerm);

            sigSegments = [sigSegments; clusterStartSec(cc) clusterEndSec(cc) clusterP(cc) obsMass(cc)]; %#ok<AGROW>
            sigLabels{end+1} = char(pairwiseLabels(cc)); %#ok<AGROW>
        end
    end

    if nObsClusters > 0
        clusterStats = table((1:nObsClusters)', clusterStartSec, clusterEndSec, obsMass(:), clusterP, clusterSignificant, pairwiseLabels, ...
            'VariableNames', {'ClusterIndex', 'ClusterStartSec', 'ClusterEndSec', 'ClusterMass', 'ClusterP', 'ClusterSignificant', 'PairwiseLabel'});
    end

    windowStats = table((1:nWindows)', permStartTimes(:), permEndTimes(:), observedF(:), observedP(:), observedSig(:), ...
        'VariableNames', {'WindowIndex', 'WindowStartSec', 'WindowEndSec', 'ObservedF', 'ObservedP', 'ClusterFormingSignificant'});
end

function [y, g] = make_y_g_from_window_mat(windowMat, ww)
    y = [];
    g = [];
    for gg = 1:length(windowMat)
        y = [y; windowMat{gg}(:, ww)]; %#ok<AGROW>
        g = [g; gg .* ones(size(windowMat{gg}, 1), 1)]; %#ok<AGROW>
    end
    validIdx = ~isnan(y) & ~isnan(g);
    y = y(validIdx);
    g = g(validIdx);
end

function [p, F] = local_anova_p_f(y, g)
    y = y(:);
    g = g(:);
    validIdx = ~isnan(y) & ~isnan(g);
    y = y(validIdx);
    g = g(validIdx);

    uniqueGroups = unique(g);
    nGroups = length(uniqueGroups);
    nTotal = length(y);

    if nGroups < 2 || nTotal <= nGroups
        p = NaN;
        F = NaN;
        return
    end

    grandMean = mean(y, 'omitnan');
    ssBetween = 0;
    ssWithin = 0;

    for ii = 1:nGroups
        thisY = y(g == uniqueGroups(ii));
        thisY = thisY(~isnan(thisY));
        if isempty(thisY)
            continue
        end
        nThis = length(thisY);
        thisMean = mean(thisY, 'omitnan');
        ssBetween = ssBetween + nThis * (thisMean - grandMean)^2;
        ssWithin = ssWithin + sum((thisY - thisMean).^2, 'omitnan');
    end

    dfBetween = nGroups - 1;
    dfWithin = nTotal - nGroups;

    if dfWithin <= 0 || ssWithin <= 0
        p = NaN;
        F = NaN;
        return
    end

    msBetween = ssBetween / dfBetween;
    msWithin = ssWithin / dfWithin;
    F = msBetween / msWithin;
    p = 1 - fcdf(F, dfBetween, dfWithin);
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

function [labelText, pairTable] = pairwise_cluster_labels(trialZ, groupNames, timeCenters, clusterWindowSec, nPerm, alphaPerm)

    nGroups = length(trialZ);
    clusterIdx = timeCenters >= clusterWindowSec(1) & timeCenters < clusterWindowSec(2);

    pairNames = strings(0, 1);
    rawP = [];
    obsDiff = [];

    for g1 = 1:nGroups-1
        for g2 = g1+1:nGroups
            x1 = mean(trialZ{g1}(:, clusterIdx), 2, 'omitnan');
            x2 = mean(trialZ{g2}(:, clusterIdx), 2, 'omitnan');
            p = permutation_pair_p(x1, x2, nPerm);
            pairNames(end+1, 1) = string(sprintf('%s vs %s', groupNames{g1}, groupNames{g2})); %#ok<AGROW>
            rawP(end+1, 1) = p; %#ok<AGROW>
            obsDiff(end+1, 1) = mean(x1, 'omitnan') - mean(x2, 'omitnan'); %#ok<AGROW>
        end
    end

    qVals = bh_fdr(rawP);
    isSig = qVals < alphaPerm;

    pairTable = table(pairNames, rawP, qVals, obsDiff, isSig, ...
        'VariableNames', {'Pair', 'PairP', 'PairQ', 'ObservedDiff', 'PairSignificant'});

    if any(isSig)
        sigPairs = pairNames(isSig);
        labelText = strjoin(sigPairs, '; ');
    else
        labelText = "cluster";
    end
end

function p = permutation_pair_p(x1, x2, nPerm)
    x1 = x1(~isnan(x1));
    x2 = x2(~isnan(x2));

    if length(x1) < 2 || length(x2) < 2
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

    p = (1 + sum(permDiffs >= observedDiff)) / (nPerm + 1);
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

function save_4stage_unit_figure(unitResults, stages, timeCenters, ptID, areaNameClean, chanNum, unitNum, unitFigureFolder)

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100 100 1200 1000]);

    tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    for ss = 1:length(stages)
        ax = nexttile;
        hold(ax, 'on');

        stage = stages(ss);
        result = unitResults{ss};

        if ~result.valid
            text(ax, 0, 0.5, sprintf('%s: skipped / not enough trials or low FR', stage.name), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 11);
            axis(ax, 'off');
            continue
        end

        nGroups = length(stage.groupNames);

        for gg = 1:nGroups
            thisColor = stage.groupColors(gg, :);

            fill(ax, [timeCenters fliplr(timeCenters)], ...
                [result.meanSmooth{gg} + result.semSmooth{gg}, ...
                 fliplr(result.meanSmooth{gg} - result.semSmooth{gg})], ...
                thisColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

            plot(ax, timeCenters, result.meanSmooth{gg}, ...
                'Color', thisColor, 'LineWidth', 1.8);
        end

        xline(ax, 0, '--k', 'HandleVisibility', 'off');
        yline(ax, 0, ':k', 'HandleVisibility', 'off');

        add_sig_segments_to_axis(ax, result.sigSegments, result.sigLabels, result.meanSmooth, result.semSmooth);

        title(ax, sprintf('%s | %s | event: %s | n = [%s] | FR %.2f Hz', ...
            stage.name, stage.longName, stage.eventFieldUsed, ...
            strjoin(arrayfun(@num2str, result.validTrials(:)', 'UniformOutput', false), ','), ...
            result.meanFRHz), 'Interpreter', 'none');

        ylabel(ax, 'z FR');
        xlabel(ax, stage.xLabel);
        xlim(ax, [-1 1]);
        legend(ax, stage.groupNames, 'Location', 'eastoutside', 'Box', 'off');
        box(ax, 'off');
    end

    sgtitle(sprintf('%s | %s | Ch %d Unit %d | 4-stage encoding', ...
        ptID, areaNameClean, chanNum, unitNum), 'Interpreter', 'none');

    pdfName = sprintf('%s_%s_ch%d_unit%d_4stage_encoding.pdf', ptID, areaNameClean, chanNum, unitNum);
    pdfPath = fullfile(unitFigureFolder, pdfName);
    exportgraphics(fig, pdfPath, 'ContentType', 'vector');
    close(fig);
end

function add_sig_segments_to_axis(ax, sigSegments, sigLabels, meanSmooth, semSmooth)

    upperVals = [];
    lowerVals = [];

    for gg = 1:length(meanSmooth)
        upperVals = [upperVals, meanSmooth{gg} + semSmooth{gg}]; %#ok<AGROW>
        lowerVals = [lowerVals, meanSmooth{gg} - semSmooth{gg}]; %#ok<AGROW>
    end

    curveHigh = max(upperVals, [], 'omitnan');
    curveLow = min(lowerVals, [], 'omitnan');

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
            'Color', [0.2 0.2 0.2], 'LineWidth', 4, 'HandleVisibility', 'off');

        if ss <= length(sigLabels) && ~isempty(sigLabels{ss})
            labelText = sprintf('%s | p=%.3f', sigLabels{ss}, sigP);
        else
            labelText = sprintf('cluster p=%.3f', sigP);
        end

        text(ax, mean([sigStart sigEnd]), sigY + 0.04*yRange, labelText, ...
            'Color', [0.1 0.1 0.1], 'FontSize', 8, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'Interpreter', 'none');
    end

    ylim(ax, [curveLow - 0.15*yRange, curveHigh + (0.30 + 0.12*size(sigSegments,1))*yRange]);
end

function save_summary_figures(allUnitSummary, OutputFolder)

    stageOrder = {'CardColor', 'CardNumber', 'Choice', 'Outcome'};

    % ====================================================================
    % 1) Stage summary: clean percentage plot, no dual y-axis, no global FDR subplot
    % ====================================================================
    stageSummary = table();

    for ss = 1:length(stageOrder)
        idx = allUnitSummary.Stage == string(stageOrder{ss});
        nValid = sum(idx);
        nSig = sum(allUnitSummary.IsSignificant(idx));

        if nValid > 0
            pctSig = 100 * nSig / nValid;
        else
            pctSig = NaN;
        end

        tmp = table(string(stageOrder{ss}), nValid, nSig, pctSig, ...
            'VariableNames', {'Stage', 'NValidUnits', 'NSignificantUnits', 'PercentSignificant'});
        stageSummary = [stageSummary; tmp]; %#ok<AGROW>
    end

    writetable(stageSummary, fullfile(OutputFolder, 'summary_stage_counts.csv'));

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100 100 750 520]);

    ax = axes(fig);
    bar(ax, stageSummary.PercentSignificant);
    ylabel(ax, '% significant units');
    xticks(ax, 1:length(stageOrder));
    xticklabels(ax, stageOrder);
    xtickangle(ax, 25);
    title(ax, 'Single-neuron stage encoding summary');
    box(ax, 'off');

    maxPct = max(stageSummary.PercentSignificant, [], 'omitnan');
    if isempty(maxPct) || isnan(maxPct)
        maxPct = 10;
    end
    yMax = min(100, max(10, ceil((maxPct + 10) / 10) * 10));
    ylim(ax, [0 yMax]);

    for ss = 1:height(stageSummary)
        labelText = sprintf('%d/%d\n%.1f%%', ...
            stageSummary.NSignificantUnits(ss), ...
            stageSummary.NValidUnits(ss), ...
            stageSummary.PercentSignificant(ss));
        text(ax, ss, stageSummary.PercentSignificant(ss) + 0.03*yMax, labelText, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontWeight', 'bold', 'FontSize', 10);
    end

    exportgraphics(fig, fullfile(OutputFolder, 'summary_stage_counts.pdf'), 'ContentType', 'vector');
    close(fig);

    % ====================================================================
    % 2) Area summary: percentages, not raw counts
    %    This replaces the old summary_area_stage_counts heatmap.
    % ====================================================================
    areaNames = unique(allUnitSummary.Area, 'stable');
    areaStageSummary = table();

    for aa = 1:length(areaNames)
        for ss = 1:length(stageOrder)
            idx = allUnitSummary.Area == areaNames(aa) & allUnitSummary.Stage == string(stageOrder{ss});
            nValid = sum(idx);
            nSig = sum(allUnitSummary.IsSignificant(idx));

            if nValid > 0
                pctSig = 100 * nSig / nValid;
            else
                pctSig = NaN;
            end

            tmp = table(areaNames(aa), string(stageOrder{ss}), nValid, nSig, pctSig, ...
                'VariableNames', {'Area', 'Stage', 'NValidUnits', 'NSignificantUnits', 'PercentSignificant'});
            areaStageSummary = [areaStageSummary; tmp]; %#ok<AGROW>
        end
    end

    writetable(areaStageSummary, fullfile(OutputFolder, 'summary_area_stage_percent.csv'));

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100 100 1150 800]);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for ss = 1:length(stageOrder)
        ax = nexttile;
        idx = areaStageSummary.Stage == string(stageOrder{ss}) & areaStageSummary.NValidUnits > 0;

        if ~any(idx)
            text(ax, 0.5, 0.5, sprintf('%s: no valid units', stageOrder{ss}), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax, 'off');
            continue
        end

        vals = areaStageSummary.PercentSignificant(idx);
        names = areaStageSummary.Area(idx);
        nSig = areaStageSummary.NSignificantUnits(idx);
        nValid = areaStageSummary.NValidUnits(idx);

        [vals, sortIdx] = sort(vals, 'descend');
        names = names(sortIdx);
        nSig = nSig(sortIdx);
        nValid = nValid(sortIdx);

        barh(ax, vals);
        set(ax, 'YDir', 'reverse');
        yticks(ax, 1:length(names));
        yticklabels(ax, names);
        set(ax, 'TickLabelInterpreter', 'none');
        xlabel(ax, '% significant units');
        title(ax, stageOrder{ss}, 'Interpreter', 'none');
        box(ax, 'off');

        maxVal = max(vals, [], 'omitnan');
        if isempty(maxVal) || isnan(maxVal)
            maxVal = 10;
        end
        xMax = min(100, max(10, ceil((maxVal + 15) / 10) * 10));
        xlim(ax, [0 xMax]);

        for ii = 1:length(vals)
            text(ax, vals(ii) + 0.02*xMax, ii, sprintf('%d/%d', nSig(ii), nValid(ii)), ...
                'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'FontSize', 9);
        end
    end

    sgtitle('Area summary: percentage of significant units, not raw counts');
    exportgraphics(fig, fullfile(OutputFolder, 'summary_area_stage_percent.pdf'), 'ContentType', 'vector');
    close(fig);

    % ====================================================================
    % 3) Participant summary: useful because units come from participants
    % ====================================================================
    patientNames = unique(allUnitSummary.Patient, 'stable');
    patientStageSummary = table();

    for pp = 1:length(patientNames)
        for ss = 1:length(stageOrder)
            idx = allUnitSummary.Patient == patientNames(pp) & allUnitSummary.Stage == string(stageOrder{ss});
            nValid = sum(idx);
            nSig = sum(allUnitSummary.IsSignificant(idx));

            if nValid > 0
                pctSig = 100 * nSig / nValid;
            else
                pctSig = NaN;
            end

            tmp = table(patientNames(pp), string(stageOrder{ss}), nValid, nSig, pctSig, ...
                'VariableNames', {'Patient', 'Stage', 'NValidUnits', 'NSignificantUnits', 'PercentSignificant'});
            patientStageSummary = [patientStageSummary; tmp]; %#ok<AGROW>
        end
    end

    writetable(patientStageSummary, fullfile(OutputFolder, 'summary_patient_stage_percent.csv'));

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100 100 1150 800]);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for ss = 1:length(stageOrder)
        ax = nexttile;
        idx = patientStageSummary.Stage == string(stageOrder{ss}) & patientStageSummary.NValidUnits > 0;

        if ~any(idx)
            text(ax, 0.5, 0.5, sprintf('%s: no valid units', stageOrder{ss}), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax, 'off');
            continue
        end

        vals = patientStageSummary.PercentSignificant(idx);
        names = patientStageSummary.Patient(idx);
        nSig = patientStageSummary.NSignificantUnits(idx);
        nValid = patientStageSummary.NValidUnits(idx);

        bar(ax, vals);
        xticks(ax, 1:length(names));
        xticklabels(ax, names);
        xtickangle(ax, 35);
        ylabel(ax, '% significant units');
        title(ax, stageOrder{ss}, 'Interpreter', 'none');
        box(ax, 'off');

        maxVal = max(vals, [], 'omitnan');
        if isempty(maxVal) || isnan(maxVal)
            maxVal = 10;
        end
        yMax = min(100, max(10, ceil((maxVal + 15) / 10) * 10));
        ylim(ax, [0 yMax]);

        for ii = 1:length(vals)
            text(ax, ii, vals(ii) + 0.03*yMax, sprintf('%d/%d', nSig(ii), nValid(ii)), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9);
        end
    end

    sgtitle('Participant summary: percentage of significant units per patient');
    exportgraphics(fig, fullfile(OutputFolder, 'summary_patient_stage_percent.pdf'), 'ContentType', 'vector');
    close(fig);
end

function [alignTime, fieldUsed] = get_event_time(eventTime, candidates, stageName)
    alignTime = [];
    fieldUsed = '';

    for ii = 1:length(candidates)
        thisField = candidates{ii};
        if isfield(eventTime, thisField)
            alignTime = double(eventTime.(thisField));
            fieldUsed = thisField;
            return
        end
    end

    warning('No event time field found for stage %s. Tried: %s', stageName, strjoin(candidates, ', '));
end


function alignTimeOut = add_rt_to_event_time(baseEventTime, rtVals, SampleRes)
    % Adds a behavioral reaction time vector to an event-time vector.
    % baseEventTime is in neural sample units.
    % RT values are auto-detected:
    %   median RT > 20  -> milliseconds, converted to seconds
    %   median RT <= 20 -> seconds

    baseEventTime = double(baseEventTime(:));
    rtVals = double(rtVals(:));

    nUse = min(length(baseEventTime), length(rtVals));
    alignTimeOut = nan(size(baseEventTime));

    rtUse = rtVals(1:nUse);
    rtMed = median(rtUse(~isnan(rtUse)), 'omitnan');

    if isempty(rtMed) || isnan(rtMed)
        warning('RT vector is empty or all NaN. Returning NaN event times.');
        return
    end

    if rtMed > 20
        % Usually jsPsych RTs are saved in milliseconds.
        rtSec = rtUse ./ 1000;
    else
        % Already in seconds.
        rtSec = rtUse;
    end

    alignTimeOut(1:nUse) = baseEventTime(1:nUse) + rtSec .* SampleRes;
end

function tf = has_bhv_field(bhvData, fieldName)
    if istable(bhvData)
        tf = ismember(fieldName, bhvData.Properties.VariableNames);
    else
        tf = isfield(bhvData, fieldName);
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
        error('bhvData must be a table or struct. Current class: %s', class(bhvData));
    end

    if isrow(x) && ~ischar(x) && ~isstring(x)
        x = x(:);
    end
end

function distribution = get_distribution_vector(bhvData)
    % Actual field in your bhvData: distribution
    if has_bhv_field(bhvData, 'distribution')
        distribution = get_bhv_field(bhvData, 'distribution');
    else
        error('Could not find bhvData.distribution. Available fields are: %s', get_bhv_field_list(bhvData));
    end
end

function mask = get_distribution_mask(distribution, distName, distIndex, nTrials)
    mask = false(nTrials, 1);

    if isnumeric(distribution) || islogical(distribution)
        distribution = double(distribution(:));
        nUse = min(length(distribution), nTrials);
        mask(1:nUse) = distribution(1:nUse) == distIndex;
        return
    end

    distribution = lower(string(distribution(:)));
    nUse = min(length(distribution), nTrials);
    d = distribution(1:nUse);

    mask(1:nUse) = d == string(distName);

    if strcmp(distName, 'uniform')
        mask(1:nUse) = mask(1:nUse) | d == "uni" | d == "u" | d == "gray" | d == "grey";
    elseif strcmp(distName, 'low')
        mask(1:nUse) = mask(1:nUse) | d == "l" | d == "orange";
    elseif strcmp(distName, 'high')
        mask(1:nUse) = mask(1:nUse) | d == "h" | d == "green";
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

function fieldList = get_bhv_field_list(bhvData)
    if istable(bhvData)
        fieldList = strjoin(bhvData.Properties.VariableNames, ', ');
    else
        fieldList = strjoin(fieldnames(bhvData), ', ');
    end
end
