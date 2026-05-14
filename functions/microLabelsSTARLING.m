function [microLabels,microPts] = microLabelsSTARLING(ptID)
% for this function you just go to the \\155.100.91.44\d\Data\ptID\Imaging\Registered
% and open ChannelMap.mat then open ChannelMap1 and LabelMap and save them
% in order, don't forget to save them in order cause they might appear
% different it always needs to be 97 to .... 

microPts = {'202421', '202511', '202512', '202518', '202521', '202522', '202601'};


if any(contains(microPts,ptID))
    switch ptID
        case{'202421'}
            microLabels= {'right orbitofrontal cortex','right mid cingulate','right anterior hippocampus'}; 

       case{'202511'}
            microLabels= {'right orbitofrontal cortex','right dorsal anterior cingulate','right amygdala'}; 

       case{'202512'}
            microLabels= {'left orbitofrontal corex', 'left anterior cingulate', 'right hippocampus'}; 

       case{'202518'}
            microLabels= {'right orbitofrontal cortex', 'right anterior cingulate', 'right anterior hippocampus'}; 

       case{'202521'}
            microLabels= {'right orbitofrontal cortex', 'right anterior cingulate', 'right anterior hippocampus'}; 

       case{'202522'}
            microLabels= {'right orbitofrontal cortex', 'right anterior cingulate'}; 

       case{'202601'}
            microLabels= {'left orbitofrontal corex', 'left anterior cingulate', 'right hippocampus'}; 
    end

else
	fprintf('\nThis patient may not have had micros...\n')
	microLabels = {};
end


%% help:


% if any(contains(microPts,ptID))
%     switch ptID
%         case {'202001'}
%             microLabels = {'right anterior cingulate','right amygdala'};
%         case {'202002'} % according to plexon, it looks like the first 2 banks took channels 97:112, and the 3rd bank is plugged into channels 121:128. This is consistent w/ ChanMap from BART_units
%             microLabels = {'left medial orbital gyrus','left anterior cingulate','left amygdala'};
%         case {'202006u'}
%             microLabels = {'left dorsal anterior cingulate','left hippocampus'}; % 
%         case {'202007'}
%             microLabels = {'left subcallosal area','left anterior hippocampus'};
%         case {'202009'}
%             microLabels = {'right gyrus rectus','right dorsal anterior cingulate'};
%         case {'202011'}
%             microLabels = {'right gyrus rectus','right parahippocampal gyrus'};
%         case {'202014'}
%             microLabels = {'right orbitofrontal','right hippocampus'};
%         case {'202015'}
%             microLabels = {'right orbitofrontal','right hippocampus'};
%         case {'202016'}
%             microLabels = {'right orbitofrontal','right hippocampus'};
%         case {'202105'}
%             microLabels = {'left orbitofrontal','right hippocampus'};
%         case {'202107'}
%             microLabels= {'left subgenual cingulate','left anterior cingulate'}; % 20240830 =TAP flipped these two based on electProj
%         case {'202110'}
%             microLabels= {'left orbitofrontal','left subgenual cingulate'};
%         case {'202114'}
%             microLabels= {'right orbitofrontal','right hippocampus'};
%         case {'202117'}
%             microLabels= {'right orbitofrontal','right hippocampus'};
%         case {'202118'}
%             microLabels= {'left orbitofrontal','left hippocampus'};
%         case {'202201'}
%             microLabels= {'left ventral cingulate','left dorsal anterior cingulate','right anterior hippocampus'};
%         case {'202202'}
%             microLabels= {'left orbitofrontal','left ventral cingulate','left dorsal anterior cingulate','right anterior hippocampus'};
%         case {'202205'}
%             microLabels= {'left ventral cingulate','left dorsal anterior cingulate','right anterior hippocampus'};
%         case {'202207'}
%             microLabels= {'left orbitofrontal','left dorsal anterior cingulate','right anterior hippocampus'};
%             %        case {'202208'}
%             %           microLabels= {'right OFC','right dorsal Anterior Cingulate','left Anterior hippocampus'};
%         case {'202209'}
%             microLabels= {'left orbitofrontal','left dorsal anterior cingulate','right entorhinal'};
%         case {'202212'}
%             microLabels= {'left orbitofrontal','left dorsal anterior cingulate','right entorhinal'};
%         case {'202214'}
%             microLabels= {'left hippocampus','left amygdala','right anterior hippocampus'};
%         case {'202215'}
%             microLabels= {'right orbitofrontal','right ventral cingulate'};
%         case {'202216'}
%             microLabels= {'left dorsal anterior cingulate','left mid cingulate','left anterior hippocampus'};
%         case {'202217'}
%             microLabels= {'left mid cingulate','left dorsal anterior cingulate','left anterior hippocampus'};
%         case {'202302'}
%             microLabels= {'right orbitofrontal','right dorsal cingulate','right anterior hippocampus'};
%         case {'202306'}
%             microLabels= {'left dorsal anterior cingulate','left entorhinal','right entorhinal'};
%         case {'202307'}
%             microLabels= {'right mid cingulate','right anterior hippocampus','left hippocampus'};
%         case {'202308'}
%             microLabels= {'left orbitofrontal','left mid cingulate','right hippocampus'};
%         case {'202309'} % THIS IS THE PT THAT SEIZED DURING BART
%             microLabels= {'right orbitofrontal','right mid cingulate','left hippocampus'};
%         case {'202311'}
%             microLabels= {'right gyrus rectus','right anterior cingulate','left hippocampus'};
%         case {'202314a'}
%             microLabels= {'left orbitofrontal cortex','left medial frontal cortex', 'right amygdala'};
%         case {'202314b'}
%             microLabels= {'left orbitofrontal cortex','left medial frontal cortex', 'right amygdala'};
%         case {'202401'}
%             microLabels= {'left orbitofrontal cortex','left ventral cingulate', 'right amygdala'};
%         case {'202405'}
%             microLabels= {'right orbitofrontal cortex','left anterior cingulate', 'left hippocampus'}; % may be left mid cingulate
%         case{'202406'}
%            microLabels= {'left orbitofrontal cortex','left anterior cingulate', 'right hippocampus'}; % could be left mid cingulate
%        case{'202407'}
%            microLabels= {'left orbitofrontal cortex','left anterior cingulate', 'left amygdala'}; % could be left mid cingulate
%        case{'202409'}
%            microLabels= {'left orbitofrontal cortex','left anterior cingulate','right hippocampus'}; 
%        case{'202413a'}
%            microLabels= {'right orbitofrontal cortex','right mid cingulate','left anterior hippocampus'}; % on ASANA it says RMACC. is that right mid?
%        case{'202413b'}
%            microLabels= {'right orbitofrontal cortex','right mid cingulate','left anterior hippocampus'}; % on ASANA it says RMACC. is that right mid?
%        case{'202417'}
%            microLabels= {'left orbitofrontal cortex','right mid cingulate','right hippocampus'}; % sheet from OR says MFG ACC.
%        case{'202418'}
%            microLabels= {'left amygdala','left hippocampus','right anterior hippocampus'}; 
%        case{'202418b'}
%            microLabels= {'left amygdala','left hippocampus','right anterior hippocampus'};
%         case{'202421'}
%             microLabels= {'right orbitofrontal cortex','right mid cingulate','right anterior hippocampus'}; 
%        case{'202422'}
%             microLabels= {'left orbitofrontal cortex','left mid cingulate','right hippocampus'}; % may need to update 'mid' cingulate
%        case{'202503'}
%             microLabels= {'left orbitofrontal cortex','right dorsal anterior cingulate','right anterior hippocampus'}; 
%        case{'202504'}
%             microLabels= {'right posterior cingulate','left superior posterior cingulate','left anterior cingulate'}; %on the document from implant, the ACC lead says SFG ACC, but there doesn't appear to be anything else labeled SFG...
%        case{'202505'}
%             microLabels= {'left orbitofrontal cortex','left dorsal anterior cingulate','right anterior hippocampus'};
%        case{'202507'}
%             microLabels= {'left orbitofrontal cortex','left dorsal anterior cingulate','left ventral anterior cingulate'}; 
%         case{'202508'}
%             microLabels= {'left orbitofrontal cortex','left mid cingulate'}; %Cadwell is lablled as MFG ACC
%         case{'202510'}
%             microLabels= {'right orbitofrontal cortex','right anterior cingulate','right posterior cingulate'};
%     end
% 
% else