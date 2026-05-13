

function [microLabelsElt,microPts,microLabelsNMM,microLocsMNI,generalLabelsElt] = microLabelsBART(ptID)

% MICROLABELSBART outputs a list of patients with microwires and
%   channel labels for microwires
%
% if you just want a list of micro patients, use 'NaP' as the input arg.

% author: EHS20200624
% added atlas labels: EHS20220217

microPts = {'202001','202002','202006u','202007','202009','202011','202014',...
    '202015','202016','202105','202107','202110','202114','202118',... % ,'202117' % issues with this patient's nev file... matbe try to re-sort???
    '202201','202202','202205','202207','202212','202214','202215','202216',... % 202208: micros but no units.
    '202217','202302','202306','202307','202308','202311','202314a','202314b','202401', '202405', '202406'...
    '202407','202409','202413a','202413b'}; %'202309' % was 202309 the pt that seized during BART?

if any(contains(microPts,ptID))
    switch ptID
        case {'202001'}
            microLabelsElt = {'right anterior cingulate','right amygdala'};
            generalLabelsElt = {'left ACC','left MTL'};
        case {'202002'} % according to plexon, it looks like the first 2 banks took channels 97:112, and the 3rd bank is plugged into channels 121:128. This is consistent w/ ChanMap from BART_units
            microLabelsElt = {'left medial orbital gyrus','left anterior cingulate','left amygdala'};
            generalLabelsElt = {'left OFC','left ACC','left MTL'}; 
        case {'202006u'}
            microLabelsElt = {'left dorsal anterior cingulate','left hippocampus'}; % 
            generalLabelsElt = {'left ACC','left MTL'};
        case {'202007'}
            microLabelsElt = {'left subcallosal area','left anterior hippocampus'};
            generalLabelsElt = {'left MFC','left MTL'};
        case {'202009'}
            microLabelsElt = {'right gyrus rectus','right dorsal anterior cingulate'};
            generalLabelsElt = {'right OFC','right ACC'};
        case {'202011'}
            microLabelsElt = {'right gyrus rectus','right parahippocampal gyrus'};
            generalLabelsElt = {'right OFC','right MTL'};
        case {'202014'}
            microLabelsElt = {'right orbitofrontal','right hippocampus'};
            generalLabelsElt = {'right OFC','right MTL'};
        case {'202015'}
            microLabelsElt = {'right orbitofrontal','right hippocampus'};
            generalLabelsElt = {'right OFC','right MTL'};
        case {'202016'}
            microLabelsElt = {'right orbitofrontal','right hippocampus'};
            generalLabelsElt = {'right OFC','right MTL'};
        case {'202105'}
            microLabelsElt = {'left orbitofrontal','right hippocampus'};
            generalLabelsElt = {'left OFC','right MTL'};
        case {'202107'}
            microLabelsElt= {'left anterior cingulate','left subgenual cingulate'}; % 20240830 =TAP flipped these two based on electProj
            generalLabelsElt = {'left ACC','left MFC'};
        case {'202110'}
            microLabelsElt= {'left orbitofrontal','left subgenual cingulate'};
            generalLabelsElt = {'left OFC','left MFC'};
        case {'202114'}
            microLabelsElt= {'right orbitofrontal','right hippocampus'};
            generalLabelsElt = {'right OFC','right MTL'};
        case {'202117'}
            microLabelsElt= {'right orbitofrontal','right hippocampus'};
            generalLabelsElt = {'right OFC','right MTL'};
        case {'202118'}
            microLabelsElt= {'left orbitofrontal','left hippocampus'};
            generalLabelsElt = {'left OFC','left MTL'};
        case {'202201'}
            microLabelsElt= {'left dorsal anterior cingulate','left ventral cingulate','right anterior hippocampus'};
            generalLabelsElt = {'left ACC','left MFC','right MTL'};
        case {'202202'}
            microLabelsElt= {'left orbitofrontal','left ventral cingulate','left dorsal anterior cingulate','right anterior hippocampus'};
            generalLabelsElt = {'left OFC','left MFC','left ACC','right MTL'};
        case {'202205'}
            microLabelsElt= {'left ventral cingulate','left dorsal anterior cingulate','right anterior hippocampus'};
            generalLabelsElt = {'left MFC','left ACC','right MTL'};
        case {'202207'}
            microLabelsElt= {'left orbitofrontal','left dorsal anterior cingulate','right anterior hippocampus'};
            generalLabelsElt = {'left OFC','left ACC','right MTL'};
            %        case {'202208'}
            %           microLabels= {'right OFC','right dorsal Anterior Cingulate','left Anterior hippocampus'};
        case {'202209'}
            microLabelsElt= {'left orbitofrontal','left dorsal anterior cingulate','right entorhinal'};
            generalLabelsElt = {'left OFC','left ACC','right MTL'};
        case {'202212'}
            microLabelsElt= {'left orbitofrontal','left dorsal anterior cingulate','right entorhinal'};
            generalLabelsElt = {'left MFC','left ACC','right MTL'};
        case {'202214'}
            microLabelsElt= {'left hippocampus','left amygdala','right anterior hippocampus'};
            generalLabelsElt = {'left MTL','left MTL','right MTL'};
        case {'202215'}
            microLabelsElt= {'right orbitofrontal','right ventral cingulate'};
            generalLabelsElt = {'right OFC','right ACC'};
        case {'202216'}
            microLabelsElt= {'left dorsal anterior cingulate','left mid cingulate','left anterior hippocampus'};
            generalLabelsElt = {'left ACC','left ACC','left MTL'}; % 20240903 was left dorsal ACC, left anterior hippocampus, left mid cingulate; flipped left mid cingulate and MTL
        case {'202217'}
            microLabelsElt= {'left anterior hippocampus','left mid cingulate','left dorsal anterior cingulate'};
            generalLabelsElt = {'left ACC','left ACC','left MTL'};  % 20240830TAP was left dorsal anterior cingulate, left anterior hippocampus, left mid cingulate. changing to left mid cingulate, left dorsal anterior cingulate, left anterior hippocampus
        case {'202302'}
            microLabelsElt= {'right orbitofrontal','right dorsal cingulate','right anterior hippocampus'};
            generalLabelsElt = {'right OFC','right ACC','right MTL'};
        case {'202306'}
            microLabelsElt= {'left dorsal anterior cingulate','left entorhinal','right entorhinal'};
            generalLabelsElt = {'left ACC','left MTL','right MTL'};
        case {'202307'}
            microLabelsElt= {'right mid cingulate','right anterior hippocampus','left hippocampus'};
            generalLabelsElt = {'right ACC','right MTL','left MTL'};
        case {'202308'}
            microLabelsElt= {'left orbitofrontal','left mid cingulate','right hippocampus'};
            generalLabelsElt = {'left OFC','left ACC','right MTL'};
        case {'202309'} % THIS IS THE PT THAT SEIZED DURING BART
            microLabelsElt= {'right orbitofrontal','right mid cingulate','left hippocampus'};
            generalLabelsElt = {'right OFC','right ACC','left MTL'};
        case {'202311'}
            microLabelsElt= {'right gyrus rectus','right anterior cingulate','left hippocampus'};
            generalLabelsElt= {'right OFC','right ACC','left MTL'}; 
        case {'202314a'}
            microLabelsElt= {'left orbitofrontal cortex','left anterior cingulate', 'right amygdala'};
            generalLabelsElt= {'left OFC','right ACC','right MTL'};
            case {'202314b'}
            microLabelsElt= {'left orbitofrontal cortex','left anterior cingulate', 'right amygdala'};
            generalLabelsElt= {'left OFC','right ACC','right MTL'};
        case {'202401'}
            microLabelsElt= {'left orbitofrontal cortex','left ventral cingulate', 'right amygdala'};
            generalLabelsElt= {'left OFC','left ACC','right MTL'};
        case {'202405'}
            microLabelsElt= {'right orbitofrontal cortex','left anterior cingulate', 'left hippocampus'}; % may be left mid cingulate
            generalLabelsElt= {'right OFC','right ACC','left MTL'};
        case{'202406'}
           microLabelsElt= {'left orbitofrontal cortex','left anterior cingulate', 'right hippocampus'}; % could be left mid cingulate
           generalLabelsElt= {'left OFC','left ACC','right MTL'}; 
       case{'202407'}
           microLabelsElt= {'left orbitofrontal cortex','left anterior cingulate', 'left amygdala'}; % could be left mid cingulate
           generalLabelsElt= {'left OFC','left ACC','left MTL'}; 
       case{'202409'}
           microLabelsElt= {'left orbitofrontal cortex','left anterior cingulate','right hippocampus'}; 
           generalLabelsElt= {'left OFC','left ACC','right MTL'};
       case{'202413a'}
           microLabelsElt= {'right orbitofrontal cortex','right mid cingulate','left anterior hippocampus'}; % on ASANA it says RMACC. is that right mid?
           generalLabelsElt= {'right OFC','right ACC','left MTL'};
       case{'202413b'}
           microLabelsElt= {'right orbitofrontal cortex','right mid cingulate','left anterior hippocampus'}; % on ASANA it says RMACC. is that right mid?
           generalLabelsElt= {'right OFC','right ACC','left MTL'};
    end
    
    % actually finding the micro labels from the channel maps, and their
    % associated locations...
%     load(sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Imaging/Registered/ChannelMap.mat',ptID))
% for ptID == 202002 the channel map from BART_preprocessed is different than channel map from BART_units
    load(sprintf('//155.100.91.44/d/Data/preProcessed/BART_units/%s/Imaging/Registered/ChannelMap.mat',ptID));
    if ~exist('ChanMap','var')
        microChans = ChannelMap1(contains(LabelMap,'m'));
        % TODO:: allow user to pick which atlas to use. so far, just using NMM.
        tmp = ElecAtlasProj(microChans,1);
        % just picking the 4th electrode...
        locIdcs = (0:length(microLabelsElt)*2).*4; % so this creates a series of indices where the final index is equal to the number of micros.
        % i.e. for 3 micros = the last index is 24 (8 microwires x 3 micros = 24)
        microLabelsNMM = tmp(locIdcs(2:2:length(locIdcs)-1))';
        
        % getting micro locations
        microLocsMNI = ElecXYZMNIProj(microChans(locIdcs(2:2:length(locIdcs)-1)),:);
         % 2:2:length(locIdcs)-1 grabs n indices (n being number of micros)
         % then it grabs one value in the middle of the microchannels and indexes ElecXYZMNIProj for that micro.

         % alternatively, we could use the MNI coordinates of distal macro?
         microNames = LabelMap(contains(LabelMap,'m'));
         distChans = microNames(contains(microNames,'1'));
         changem2b = @(str) ['b',str(2:end)];
         distChans = cellfun(changem2b,distChans,'UniformOutput',false);
         if strcmp(ptID,'202207')
             distChans{3} = 'bRAHIP1';
         end
         for ch = 1:length(distChans)
             try
             macroChans(ch) = ChannelMap1(contains(LabelMap,distChans{ch})); % this works unless one of them is NA
             catch
                 if ch == 1
                     removeb = @(str) [str(2:end)];
                     distChans = cellfun(removeb,distChans,'UniformOutput',false);
                 end
                 macroChans(ch) = ChannelMap1(contains(LabelMap,distChans{ch}) & ~contains(LabelMap,'m')); % this returns the macro and micro!
             end
             if ~isnan(macroChans(ch))
                 continue
             else % have to walk back until we get non NaN
                 [r,c] = ind2sub(size(ChannelMap1),find(contains(LabelMap,distChans{ch})));
                 macroChans(ch) = min(ChannelMap1(:,c));
             end
         end

         microLocsMNI = ElecXYZMNIProj(macroChans,:);

    else
        try
            microChans = ChanMap.ChannelMap1(contains(ChanMap.LabelMap,'m'));
            % TODO:: allow user to pick which atlas to use. so far, just using NMM.
            tmp = ChanMap.ElecNMMProj(microChans);
            % just picking the 4th electrode...
            locIdcs = (0:length(microLabelsElt)*2).*4;
            microLabelsNMM = tmp(locIdcs(2:2:length(locIdcs)-1))';
            
            % getting micro locations
            microLocsMNI = ChanMap.ElecXYZMNIProj(microChans(locIdcs(2:2:length(locIdcs)-1)),:);

            % alternatively, we could use the MNI coordinates of distal macro?
            microNames = ChanMap.LabelMap(contains(ChanMap.LabelMap,'m'));
         distChans = microNames(contains(microNames,'1'));
         changem2b = @(str) ['b',str(2:end)];
         distChans = cellfun(changem2b,distChans,'UniformOutput',false);
         if strcmp(ptID,'202207')
             distChans{3} = 'bRAHIP1';
         end
         for ch = 1:length(distChans)
             try
             macroChans(ch) = ChanMap.ChannelMap1(contains(ChanMap.LabelMap,distChans{ch})); % this works unless one of them is NA
             catch
                 if ch == 1
                     removeb = @(str) [str(2:end)];
                     distChans = cellfun(removeb,distChans,'UniformOutput',false);
                 end
                 macroChans(ch) = ChanMap.ChannelMap1(contains(ChanMap.LabelMap,distChans{ch}) & ~contains(ChanMap.LabelMap,'m')); % this returns the macro and micro!
             end
             if ~isnan(macroChans(ch))
                 continue
             else % have to walk back until we get non NaN
                 [r,c] = ind2sub(size(ChanMap.ChannelMap1),find(contains(ChanMap.LabelMap,distChans{ch})));
                 macroChans(ch) = min(ChanMap.ChannelMap1(:,c));
             end
         end

         microLocsMNI = ChanMap.ElecXYZMNIProj(macroChans,:);

            
        catch
            microLabelsNMM = microLabelsElt;
        end
    end
elseif strcmp(ptID,'NaP')
    fprintf('\njust returning list of patients.\n')
    microLabelsElt = {};
    microLabelsNMM = {};
    microLocsMNI = [];
else
    fprintf('\nThis patient may not have had micros...\n')
    microLabelsElt = {};
    microLabelsNMM = {};
    microLocsMNI = [];
end
