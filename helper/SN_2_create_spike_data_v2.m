% this code is for the purpose of careating spikes data based on the
% dimensions: trials, units, spikes.
% THIS IS A CRAPPY VERSION FULL OF MISTAKES! AFTER YOU TALKED TO ELLIOT,
% YOU HAVE VERSION 2 OF THIS THAT WORKS WELL. NICE JOB NILOO KHARE!

% PAY ATTENTION: ChanUnitTimestamp MATRIX has many rows and 3 columns. Each
% row mean 1 spike that happened. column 1: the channel that spike has
% happened in! column 2: the unit of the correspondence channel that the
% spike has happened in. column 3: the time of the spike in seconds. 
% AUTHOR: Nill



clc;
clear;
close all;

%% loading neural and beahvioral data.


microPts = {'202001','202002','202006','202007','202009','202011','202014',...
    '202015','202016','202105','202107','202110','202114','202118',... % ,'202117' : issues with this patient's nev file... matbe try to re-sort???
    '202201','202202','202205','202207','202212','202214','202215','202216',... % 202208: micros but no units.
    '202217','202302','202306','202307','202308','202311','202314a','202314b','202401', '202405', '202406'...
    '202407'}; %'202309' % was 202309 the pt that seized during BART?

% microPts = { '202212'};

outputFolderName = '\\155.100.91.44\d\Data\Nill\BART\Spikes\'; 

%%
for pt = 1:length(microPts)

    ptID = microPts{pt};
    disp(ptID);
    
    nevList = dir(sprintf('//155.100.91.44/d/Data/preProcessed/BART_units/%s/Data/*.nev',ptID));
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
    matFile = sprintf('//155.100.91.44/d/Data/preProcessed/BART_units/%s/Data/%s.bartBHV.mat',ptID, ptID);
    load(matFile)

    pointsEarned = [data.points];
    
    % standard 3-D [chan, unit, timestamp (seconds)] matrix.
    ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];
    
    % channel deets.
    inclChans = unique(ChanUnitTimestamp(:,1));

    microLabels = microLabelsBART(ptID);
    inclChans(inclChans-96>length(microLabels)*8) = []; % magic numbers for recording on bank D and number of BF micros.
    nChans = length(inclChans);
    
    
    
   % task parameters in chronological order..
    % There aren't any trigs that == 4
    balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
    inflateTimes = trigTimes(trigs==23);
    balloonIDs = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
    outcomeTimes = trigTimes(trigs==25 | trigs==26);
    outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
    nTrials = min([length(outcomeType) length(balloonIDs)]);
    balloonTimes = balloonTimes(1:nTrials);
    ReactionTimes  = [data.rt];
    ReactionTimes = ReactionTimes(1:nTrials);
    pre = 2;
    post = 4;

    binWidth = 1;
    Fspikes = 1000; 
    SpikeLength = 6001;
    MaxnUnits = 5; % hypothetically there are 30 units

    SpikeData = nan(nChans,MaxnUnits, nTrials, SpikeLength);
    

    % looping over Channels
    for ch = 1:nChans
        % looping over number of units in the AP data
        

        Units = unique(ChanUnitTimestamp(inclChans(ch).*ones(size(ChanUnitTimestamp,1),1)==ChanUnitTimestamp(:,1),2));
        if ~any(Units == 0 & Units == 255)
         
            nUnits = length(Units);
            for un = 1:nUnits
                
                % getting unit times for the current channel and unit.
                unitTimes = ChanUnitTimestamp(ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un,3); % in seconds
                if ~isempty(unitTimes)
                    
                    % cue aligned spikes.
                    % loooping over trials
                    for trial = 1:nTrials
                        % disp(ch)
                            % putting the data in a structure
                            SpikesInTrialTimes = unitTimes(unitTimes>balloonTimes(trial)-pre & unitTimes<balloonTimes(trial)+post) - repmat(balloonTimes(trial)-pre,length(unitTimes(unitTimes>balloonTimes(trial)-pre & unitTimes<balloonTimes(trial)+post)),1);
                            SpikeStruct.Spikes.channel(ch).unit(un).trial(trial).spikesTimes = SpikesInTrialTimes;
                            SpikesInTrialTimesScaled = ceil(SpikesInTrialTimes*Fspikes);
                            SpikesInTrial = zeros(1,SpikeLength);
                            SpikesInTrial(SpikesInTrialTimesScaled) = 1;
                            SpikeData(ch,un,trial,:) = 0; % first turn nans to zeros
                            SpikeData(ch,un,trial,:) = SpikesInTrial;
                            if sum(SpikesInTrial>0)
                                SpikeStruct.Spikes.channel(ch).unit(un).trial(trial).spikes = SpikesInTrial;
                            else
                                SpikeStruct.Spikes.channel(ch).unit(un).trial(trial).spikes = [];
                            end
    
    
                        clear SpikesInTrialTimes SpikesInTrial SpikesInTrialTimesScaled
                    end 
                end
                clear unitTimes Units
            end
        end % end of if (Units~=0 && Units~=255) clause
    end

    % sometimes when some channels are empty and are at the end of
    % channels, there is a posibility that it's not saved in the final data
    % so we need to make sure that the number of channels in inclChans



    allNaN = all(all(all(all(isnan(SpikeData), 1), 2), 3), 4);

    if allNaN
        disp('No spikes found!');
    else
        if (length(SpikeStruct.Spikes.channel)<nChans)
            inclChans = inclChans(1:length(SpikeStruct.Spikes.channel));
        end
    
    
        SpikeStruct.microLabels = microLabels;
        SpikeStruct.balloonTimes = balloonTimes;
        SpikeStruct.ReactionTimes = ReactionTimes; % we need to save reaction time in order to delete trials with RTs more than 10 in the sigle neuron analysis
        SpikeStruct.inclChans = inclChans;
        save([outputFolderName ptID '.spikes.mat'],'SpikeStruct');
    end


    clear SpikeData microLabels balloonTimes nTrials SpikeStruct UnitsNumbers balloonIDs spikes inclChans ReactionTimes
 
end


%%


