function [microLabels,microPts] = microLabelsSTARLING(ptID)

% MICROLABELSANToutputs channel labels for microwires

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};


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


