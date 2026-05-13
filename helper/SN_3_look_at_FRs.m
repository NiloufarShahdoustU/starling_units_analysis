% this code is for the purpose of creating firing rates and take a look at
% raster plots and firing rates. the output is in the outputFolderName

% data.
% AUTHOR: Nill

clear;
clc;
close all;
warning('off','all');
%% loading neural and event data
 
% getting the numbers of different patients

inputFolderName = '\\155.100.91.44\d\Data\Nill\BART\Spikes';
fileList = dir(fullfile(inputFolderName, '*.spikes.mat'));
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\SingleNeuron\output\0_FRs';
PatientsNum = length(fileList);
SpikeDataLength = 6001;
PlotRows = 2;
for pt = 1:PatientsNum
% for pt = 1:1
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    Spikedata = [inputFolderName '\' ptID '.spikes.mat'];
    load(Spikedata);
    disp(['Processing patient ID: ' ptID]);
    Spikes = SpikeStruct.Spikes; % 4d, chan, unit, trial, spikes
    microLabels = SpikeStruct.microLabels;
    inclChans = SpikeStruct.inclChans;
    nChans = length(inclChans);
    nTrials = length(SpikeStruct.balloonTimes);

    
    for ch=1:nChans
        

        % here we check if there is any data for this channel:
        nUnits = length(Spikes.channel(ch).unit);
        if nUnits > 0 
            disp('Accepted channel: '+ string(inclChans(ch)));
            figure('Units', 'normalized', 'Position', [0.1, 0, 0.9, 0.6], 'Visible','off');
            CurrentChanBasedSpike = Spikes.channel(ch);
            
            for un=1:nUnits % for all the units that a channel has I'm gonna have a raster plot and a psth
                TrialsSpikes = nan(nTrials,SpikeDataLength);
                % populating TrialsSpikes
                
                    for trial = 1:nTrials
                        if(length(CurrentChanBasedSpike.unit(un).trial)>0)
                            CurrentChanBasedSpikeTrial = CurrentChanBasedSpike.unit(un).trial(trial).spikes;
        
                   
                            if (length(CurrentChanBasedSpikeTrial)>0)
                                TrialsSpikes(trial,:) = CurrentChanBasedSpikeTrial;
                            end
                        end
                    end
    
                    rowsWithNaN = any(isnan(TrialsSpikes), 2);
                    % Delete these nan trials
                    TrialsSpikes(rowsWithNaN, :) = [];
    
                                    
                    % Raster plot for this unit
                    subplot(PlotRows, nUnits, un);
                    imagesc(TrialsSpikes);
                    colormap(flipud(gray));
                    set(gca, 'YDir', 'reverse'); 
                    yticks = get(gca, 'YTick');
                    set(gca, 'YTickLabel', flipud(get(gca, 'YTickLabel')));
                    xlabel('Time (ms)');
                    xlim([0 6000]);
                    ylabel('Trials');
                    % title(sprintf('Unit %d Raster', un));
                    box off;  
                    
                    hold on;  % Add this to overlay on the imagesc
                    [row, col] = find(TrialsSpikes);
                    plot(col, row, 'k.', 'MarkerSize', 5);  % Plot larger dots
                    hold off;
    
    
    
    
                    
                    % PSTH for this unit
                    subplot(PlotRows, nUnits, nUnits + un); 
                    binWidth = 20;
                    numBins = round(size(TrialsSpikes, 2) / binWidth);
                    binnedSpikes = zeros(1, numBins);
                    for i = 1:numBins
                        binStart = (i-1) * binWidth + 1;
                        binEnd = i * binWidth;
                        binnedSpikes(i) = sum(sum(TrialsSpikes(:, binStart:binEnd)));
                    end
                    psth = binnedSpikes / size(TrialsSpikes, 1);
                    timeVector = linspace(binWidth, size(TrialsSpikes, 2), numBins);
                    plot(timeVector, psth);
                    xlabel('Time (ms)');
                    xlim([0 6000]);
                    ylabel('Average spikes per bin');
                    % title('PSTH');
                    box off;  

    
                clear TrialsSpikes TrialsSpikes rowsWithNaN

            end % unit for
                annotation('textbox', [0.1, 0.9, 0.8, 0.1], 'String', 'patient '+ string(ptID) + ', channel '+ string(inclChans(ch)), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14);

                % for showing brain areas, it depends on the channel number
                % so let's do this: channel numbers start from 104 and end
                % in 120. 
                if inclChans(ch)<= 104 % first micro label
                   annotation('textbox', [0.1, 0.87, 0.8, 0.1], 'String', microLabels(1), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14);
                end
                if (inclChans(ch)>= 105) && (inclChans(ch)<= 112)  % first micro label
                   annotation('textbox', [0.1, 0.87, 0.8, 0.1], 'String', microLabels(2), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14);
                end

                if (inclChans(ch)>= 113) && (inclChans(ch)<= 120) % first micro label
                   annotation('textbox', [0.1, 0.87, 0.8, 0.1], 'String', microLabels(3), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14);
                end

                if (inclChans(ch)>= 121) % first micro label
                   annotation('textbox', [0.1, 0.87, 0.8, 0.1], 'String', microLabels(4), 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14);
                end
               

                set(gca, 'box', 'off', 'tickdir', 'out');
                set(gcf, 'Units', 'inches');
                screenposition = get(gcf, 'Position');
                set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
                filename = string(ptID) + '_ch'+ string(inclChans(ch)) + '_FR';
                saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
        end % if ChannelEmptyOrNot > 0 

        % save each figure based on the patient and channels.
        % I am saving different patients and channels. 

    clear CurrentChanBasedSpike 

    end % channel for     
    clear Spikes microLabels inclChans

end % pt for



%% debug part


% aaaaa = length(Spikes.channel(1).unit(2));
% 
% aaaa = Spikes.channel(2).unit(1);
