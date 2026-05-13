% this code is for the purpose of careating spikes data based on the
% dimensions: trials, units, spikes.


% PAY ATTENTION: ChanUnitTimestamp MATRIX has many rows and 3 columns. Each
% row mean 1 spike that happened. column 1: the channel that spike has
% happened in! column 2: the unit of the correspondence channel that the
% spike has happened in. column 3: the time of the spike in seconds. 
% AUTHOR: Nill


clc;
clear;
close all;

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};
outputFolderName = '\\155.100.91.44\d\Data\Nill\starling\spikes\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

% for pt = 1:length(microPts)

for pt = 1:1
    ptID = microPts{pt};
    disp(ptID);

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

    % TimeRes = NEV.MetaTags.TimeRes;
    ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp))'];
    inclChans = unique(ChanUnitTimestamp(:,1));
    % microLabels = microLabelsSTARLING(ptID);
    % inclChans(inclChans-96>length(microLabels)*8) = []; % magic numbers for recording on bank D and number of BF micros.
    % nChans = length(inclChans);
end



%% DEEEEBUUUG

