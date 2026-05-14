clc;
clear;
close all;

totalRewardShowTime = 1000;
cardsShowTime = 1000;
choiceTime = 3500;
timeoutMessageTime = 2000;

input_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\raw');
trigs_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\spikes\trigs');
output_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\spikes\eventTimes');

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

d = dir(input_folder);

isub = [d(:).isdir];
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = [];
ptIDs = string(subFolders);

for p = 1:numel(ptIDs)

    ptID = ptIDs{p};
    fprintf('\n--- Processing ptID: %s ---\n', ptID);

    input_folder_pt = fullfile(input_folder, ptID);

    bhvFiles = dir(fullfile(input_folder_pt, 'task_data*.csv'));

    if isempty(bhvFiles)
        warning('No task_data csv found for ptID %s. Skipping...', ptID);
        continue;
    end

    bhvFile = fullfile(bhvFiles(1).folder, bhvFiles(1).name);
    bhvData = readtable(bhvFile);

    bhv_ITI = bhvData.interTrialInterval;
    timeout_idx = find(strcmp(bhvData.trialType, 'timeout'));

    trigFile = fullfile(trigs_folder, [ptID '_trigs.mat']);

    if ~exist(trigFile, 'file')
        trigFile = fullfile(trigs_folder, 'trigs.mat');
    end

    if ~exist(trigFile, 'file')
        warning('No trigs file found for ptID %s. Skipping...', ptID);
        continue;
    end

    load(trigFile, 'trigs');

    trialStartTime         = trigs(2, trigs(1,:) == 1);
    cardShowTime           = trigs(2, trigs(1,:) == 2);
    instructionMessageTime = trigs(2, trigs(1,:) == 3);
    flipSpaceTime          = trigs(2, trigs(1,:) == 4);
    choiceAndFeedbackTime  = trigs(2, trigs(1,:) == 5);
    totalRewardTime        = trigs(2, trigs(1,:) == 6);

    cfTime = choiceAndFeedbackTime;
    trTime = totalRewardTime;

    for k = 1:numel(timeout_idx)
        idx = timeout_idx(k);

        cfTime = [cfTime(1:idx-1), NaN, cfTime(idx:end)];
        trTime = [trTime(1:idx-1), NaN, trTime(idx:end)];
    end

    choiceAndFeedbackTime = cfTime;
    totalRewardTime = trTime;

    nTrials = length(trialStartTime);
    photodiodeITI = nan(1, nTrials);

    for i = 1:nTrials-1

        if ismember(i, timeout_idx)
            trialEndTime = flipSpaceTime(i) + choiceTime + timeoutMessageTime + cardsShowTime;
        else
            trialEndTime = totalRewardTime(i) + totalRewardShowTime;
        end

        photodiodeITI(i+1) = trialStartTime(i+1) - trialEndTime;
    end

    save(fullfile(output_folder, [ptID '_eventTimes.mat']), ...
        'trialStartTime', ...
        'cardShowTime', ...
        'instructionMessageTime', ...
        'flipSpaceTime', ...
        'choiceAndFeedbackTime', ...
        'totalRewardTime', ...
        'photodiodeITI', ...
        'bhv_ITI', ...
        'bhvData');

end