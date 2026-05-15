% this code is for the purpose of creating spikes data based on the
% dimensions: trials, units, spikes.

% when you read the NEV file and you go to NEV.Data.Spikes.Unit,
% you'll see number 255 too between units, that means that there was
% no waveform saved at that time for any unit.

% attention: when new patient is added you also need to update microLabelsSTARLING
% function

% PAY ATTENTION: ChanUnitTimestamp MATRIX has many rows and 3 columns. Each
% row mean 1 spike that happened. column 1: the channel that spike has
% happened in! column 2: the unit of the correspondence channel that the
% spike has happened in. column 3: the time of the spike in samples.
% AUTHOR: Nill

% in this analysis 1 second is 30000 datapoints!!!

clc;
clear;
close all;

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};
difChanNumberPts = {'202421', '202511', '202512', '202518'};

outputFolderName = '\\155.100.91.44\d\Data\Nill\starling\spikes\';
eventTimesFolder = '\\155.100.91.44\d\Data\Nill\starling\spikes\eventTimes';
spikeDataFolder = '\\155.100.91.44\d\Data\Nill\starling\spikes\spike_data';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

if ~exist(spikeDataFolder, 'dir')
    mkdir(spikeDataFolder);
end

for pt = 1:length(microPts)
% for pt = 1:1

    ptID = microPts{pt};
    disp(ptID);

    eventTimesFile = fullfile(eventTimesFolder, sprintf('%s_eventTimes.mat', ptID));

    if ~exist(eventTimesFile, 'file')
        error('EventTimes file not found for ptID %s: %s', ptID, eventTimesFile);
    end

    eventTimesData = load(eventTimesFile);

    baseFolder = sprintf('\\\\155.100.91.44\\d\\Data\\UIC%s', ptID);

    folderList = dir(baseFolder);
    folderNames = {folderList([folderList.isdir]).name};

    starlingIdx = find(strcmpi(folderNames, 'starling'));

    if isempty(starlingIdx)
        error('No Starling folder found for patient %s', ptID)
    end

    nevFolder = fullfile(baseFolder, folderNames{starlingIdx});

    nevList = dir(fullfile(nevFolder, '*sortedNS*.nev'));

    if length(nevList) > 1
        error('Many sortedNS nev files available for patient %s. Please specify...', ptID)
    elseif isempty(nevList)
        error('No sortedNS nev file found for patient %s...', ptID)
    else
        nevFile = fullfile(nevList.folder, nevList.name);
    end

    disp(nevFile)

    NEV = openNEV(nevFile, 'overwrite');

    TimeRes = NEV.MetaTags.TimeRes;

    ChanUnitTimestamp = [ ...
        double(NEV.Data.Spikes.Electrode)' ...
        double(NEV.Data.Spikes.Unit)' ...
        double(NEV.Data.Spikes.TimeStamp)' ...
        ];

    waveForms = NEV.Data.Spikes.Waveform;

    inclChans = unique(ChanUnitTimestamp(:,1));

    microLabelsOriginal = microLabelsSTARLING(ptID);

    microLabels = cell(length(microLabelsOriginal) * 8, 1);
    
    idx = 1;
    
    for m = 1:length(microLabelsOriginal)
    
        for k = 1:8
    
            microLabels{idx} = microLabelsOriginal{m};
            idx = idx + 1;
    
        end
    
    end

    if any(strcmp(ptID, difChanNumberPts))
        inclChans(inclChans - 96 > length(microLabelsOriginal) * 8) = [];
        nChans = length(inclChans);
    else
        inclChans(inclChans - 192 > length(microLabelsOriginal) * 8) = [];
        nChans = length(inclChans);
    end

    NumberOfUnits = zeros(size(inclChans));

    for ch = 1:length(inclChans)

        thisChan = inclChans(ch);

        unitsThisChan = unique(ChanUnitTimestamp(ChanUnitTimestamp(:,1) == thisChan, 2));

        unitsThisChan(unitsThisChan == 255) = [];

        NumberOfUnits(ch) = length(unitsThisChan);

    end

    spikeData = struct();

    spikeData.ptID = ptID;
    spikeData.eventTimes = eventTimesData;
    spikeData.ChanUnitTimestamp = ChanUnitTimestamp;
    spikeData.inclChans = inclChans;
    spikeData.microLabels = microLabels;
    spikeData.nChans = nChans;
    spikeData.NumberOfUnits = NumberOfUnits;
    spikeData.TimeRes = TimeRes;
    spikeData.nevFile = nevFile;
    spikeData.WaveFroms = waveForms;

    saveFile = fullfile(spikeDataFolder, sprintf('%s_spikeData.mat', ptID));

    save(saveFile, 'spikeData', '-v7.3');

    fprintf('Saved spikeData for ptID %s to:\n%s\n', ptID, saveFile);

end