% this code is for the purpose of careating spikes data based on the
% dimensions: trials, units, spikes.
% THIS IS A CRAPPY VERSION FULL OF MISTAKES! AFTER YOU TALKED TO ELLIOT,
% YOU HAVE VERSION 2 OF THIS THAT WORKS WELL. NICE JOB NILOO KHARE!
% AUTHOR: Nill

clc;
clear;
close all;
%% loading neural and beahvioral data.


% microPts = {'202001','202002','202006','202007','202009','202011','202014',...
%     '202015','202016','202105','202107','202110','202114','202118',... % ,'202117' : issues with this patient's nev file... matbe try to re-sort???
%     '202201','202202','202205','202207','202212','202214','202215','202216',... % 202208: micros but no units.
%     '202217','202302','202306','202307','202308','202311','202314a','202314b','202401', '202405', '202406'...
%     '202407'}; %'202309' % was 202309 the pt that seized during BART?


microPts = { '202007'};

outputFolderName = '\\155.100.91.44\d\Data\Nill\BART\Spikes\'; 

for pt = 1:length(microPts)

    ptID = microPts{pt};
    disp(ptID);
    
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
    % we want to consider labels based on the length of the microlabels
    logicalIndex = (ChanUnitTimestamp_units ~= 0) & (ChanUnitTimestamp_units ~= 255);
    ChanUnitTimestamp_new = ChanUnitTimestamp(logicalIndex,:);
    
    
    %%
    SignalLength = 6001; % I want it to have the same size as lfp signal (-2 to 4 sec of balloon onset)
    % looping over number of units in the AP data
    [UnitsNumbers, ~, indices] = unique(ChanUnitTimestamp_new(:,2));

    %here I want to remove the units that occured less than 30 times.
    
    UnitsNumbers_counts = accumarray(indices, 1);
    threshold = 100;
    valid_units_filter = UnitsNumbers_counts >= threshold;
    UnitsNumbers = UnitsNumbers(valid_units_filter);


    nUnits = length(UnitsNumbers);
    disp("number of units: " + string(nUnits))
    SpikeData = zeros(nTrials,nUnits,SignalLength);
    ChannelData = zeros(nTrials,nUnits,SignalLength);
    Fnew = 1000;
    % Scaling the array
    
    for trial=1:nTrials
        for un = 1:nUnits
            StimOnsetTime = balloonTimes(trial);
            unitTimes = ChanUnitTimestamp_new(ChanUnitTimestamp_new(:,2)==UnitsNumbers(un),3);
            unitChanels = ChanUnitTimestamp_new(ChanUnitTimestamp_new(:,2)==UnitsNumbers(un),1);

            unitTimesStimChunk = unitTimes(unitTimes>StimOnsetTime-pre & unitTimes<StimOnsetTime+post);
            unitChannelsTimesStimChunk = unitChanels(unitTimes>StimOnsetTime-pre & unitTimes<StimOnsetTime+post);

            repmating = repmat(StimOnsetTime-pre,length(unitTimesStimChunk),1);


            CurrentTrialSpikeTime =  unitTimesStimChunk - repmating; 
            CurrentTrialSpikeTimeScaled = round(CurrentTrialSpikeTime*Fnew);



            % dealing with duplicates, because we don't want two spikes in
            % 1 milisecond!
            [~, unique_idx] = unique(CurrentTrialSpikeTimeScaled, 'stable');
            CurrentTrialSpikeTimeScaled = CurrentTrialSpikeTimeScaled(unique_idx);


            CurrentTrialSpikeTimeScaled(CurrentTrialSpikeTimeScaled>SignalLength | CurrentTrialSpikeTimeScaled<1) = [];
            SpikesHappened = zeros(1,SignalLength);
            SpikesHappened(CurrentTrialSpikeTimeScaled) = 1;

            if(length(CurrentTrialSpikeTimeScaled)~= sum(SpikesHappened))
                disp("oooopps! there are duplicate values in CurrentTrialSpikeTimeScaled")
                disp(length(CurrentTrialSpikeTimeScaled));
                disp(sum(SpikesHappened));
            end

            SpikeData(trial,un, :) = SpikesHappened;
        end % looping over units
    end % looping over trials
    % clear trial un StimOnsetTime unitTimes unitTimesStimChunk CurrentTrialSpikeTime
    % clear CurrentTrialSpikeTimeScaled SpikesHappened unitChanels unitChanelsStimChunk
    % % 
    
    SpikeStruct.Spikes = SpikeData;
    SpikeStruct.microLabels = microLabels;
    SpikeStruct.balloonTimes = balloonTimes;
    SpikeStruct.UnitsNumbers = UnitsNumbers;

    
    
    save([outputFolderName ptID '.spikes.mat'],'SpikeStruct')

    % clear SpikeData microLabels balloonTimes balloonType nTrials SpikeStruct UnitsNumbers balloonIDs UnitsNumbers_counts valid_units_filter

end
