% this code is for the purpose of getting familiar with the single neuron
% data.

% THIS CODE IS NOT CORRECT, YOU NEED TO CHANGE IT WHEN YOU TALKED TO
% ELLIOT. AND NOW THE SECOND VERION OF n2_create_spike_data is correct. 
% AUTHOR: Nill

clc;
clear;
close all;
%% loading neural and beahvioral data.

ptID = '202314a';


nevList = dir(sprintf('//155.100.91.44/d/Data/preProcessed/BART_preprocessed/%s/Data/*.nev',ptID));
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end


% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;

% loading behavioral matFile
matFile = sprintf('//155.100.91.44/d/Data/preProcessed/BART_preprocessed/%s/Data/%s.bartBHV.mat',ptID, ptID);
load(matFile)
pointsEarned = [data.points];

% standard 3-D [chan, unit, timestamp (seconds)] matrix.
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];

% channel deets.
inclChans = unique(ChanUnitTimestamp(:,1));
microLabels = microLabelsBART(ptID);
inclChans(inclChans-96>length(microLabels)*8) = []; % magic numbers for recording on bank D and number of BF micros.
nChans = length(inclChans);

%% I wanna find trial durations
% I did this out of curiosity!!! :D


TrialStartTimes  =  trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
TrialEndTimes  = trigTimes(trigs==120);

% removing first and last trials for safety margin
TrialEndTimes = TrialEndTimes(2:end-1);
TrialStartTimes = TrialStartTimes(2:end-1);

TrialDuration = TrialEndTimes - TrialStartTimes;
disp("--------------------------------------------------------------------------------------------------");
disp("min trial duration (this should be at least 1 sec and not less than 1): ");
disp(min(TrialDuration)); % this value should be at least 1 and should not be less than 1. 
disp("--------------------------------------------------------------------------------------------------");





%% task parameters in chronological order..
% There aren't any trigs that == 4
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
balloonIDs = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
nTrials = min([length(outcomeType) length(balloonIDs)]);
balloonTimes = balloonTimes(1:nTrials);
balloonIDs = balloonIDs(1:nTrials);
pre = 2;
post = 4;
%% Cleaning ChanUnitTimestamp
% there are a couple of unit numbers in ChanUnitTimestamp, we want to
% delete the lines that the amount in 255 because it means no unit
% Logical indexing approach
ChanUnitTimestamp_units = ChanUnitTimestamp(:,2); 
logicalIndex = (ChanUnitTimestamp_units ~= 255) & (ChanUnitTimestamp_units ~= 0);
ChanUnitTimestamp_new = ChanUnitTimestamp(logicalIndex,:);


%%
SignalLength = 3001; % I want it to have the same size as lfp signal (-2 to 4 sec of balloon onset)
% looping over number of units in the AP data
nUnits = length(unique(ChanUnitTimestamp_new(:,2)));
SpikeData = zeros(nTrials,nUnits,SignalLength);
% Old range (min and max)
oldMin = 0;
oldMax = 6;
% New range (min and max)
newMin = 1;
newMax = 3001;

% Scaling the array

for trial=1:nTrials
    for un = 1:nUnits
        StimOnsetTime = balloonTimes(trial);
        unitTimes = ChanUnitTimestamp_new(ChanUnitTimestamp_new(:,2)==un,3);
        unitTimesStimChunk = unitTimes(unitTimes>StimOnsetTime-pre & unitTimes<StimOnsetTime+post);
        CurrentTrialSpikeTime =  unitTimesStimChunk - repmat(StimOnsetTime-pre,length(unitTimesStimChunk),1); 
        CurrentTrialSpikeTimeScaled = floor(((CurrentTrialSpikeTime - oldMin) / (oldMax - oldMin)) * (newMax - newMin) + newMin);
        SpikesHappened = zeros(1,SignalLength);
        SpikesHappened(CurrentTrialSpikeTimeScaled) = 1;
        SpikeData(trial,un, :) = SpikesHappened;
    end % looping over units
end % looping over trials


%%
% Example matrix with 0s and 1s
data = squeeze(SpikeData(:,1,:));

% Number of trials and time points
[num_trials, num_timepoints] = size(data);

% Create figure
figure;

% Subplot 1: Raster Plot
subplot(2, 1, 1);
hold on;
imagesc(data);
xlim([1,3000]);
ylim([1 num_trials+1]);
xlabel('Time (ms)');
ylabel('Trial');
title('Raster Plot');
hold off;

% Subplot 2: PSTH (using a curve)
subplot(2, 1, 2);
% Sum across all trials to get spike counts at each time point
spike_counts = sum(data, 1);
% Smooth the spike counts
smoothed_counts = smooth(spike_counts, 200);

% Normalize the smoothed counts between 0 and 1
max_count = max(smoothed_counts); % Find the maximum count
normalized_counts = smoothed_counts / max_count; % Normalize by dividing by the max count

% Plotting the normalized PSTH with a smooth line
plot(1:num_timepoints, normalized_counts, 'k', 'LineWidth', 2);
xlim([1,3000]);
xlabel('Time (ms)');
ylabel('Normalized Spike Rate');
title('Normalized PSTH');

