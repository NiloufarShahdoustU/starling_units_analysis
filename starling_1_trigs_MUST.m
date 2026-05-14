clc;
clear;
close all;

differentPatients = {'202421', '202509', '202511', '202512'};

input_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\raw');
output_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\spikes\trigs');

d = dir(input_folder);

isub = [d(:).isdir];
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = [];
ptIDs = string(subFolders);

for p = 5:numel(ptIDs)

    ptID = ptIDs{p};
    fprintf('\n--- Processing ptID: %s ---\n', ptID);

    input_folder_pt = fullfile(input_folder, ptID);

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    if any(strcmp(ptID, differentPatients))
        nevList = dir(fullfile(input_folder_pt, '*.nev'));
    else
        nevList = dir(fullfile(input_folder_pt, 'nsp_photodiode_data', '*.nev'));
    end

    if length(nevList) > 1
        error('many nev files available for this patient. Please specify...')
    elseif isempty(nevList)
        error('no nev files found...')
    else
        nevFile = fullfile(nevList.folder, nevList.name);
    end

    [nevPath, nevName, ~] = fileparts(nevFile);
    ns5File = fullfile(nevPath, [nevName '.ns5']);

    NS5 = openNSx(ns5File);

    original_freq = NS5.MetaTags.SamplingFreq;

    lp_cutoff = 60;
    [b, a] = butter(4, lp_cutoff/(original_freq/2), 'low');

    photodiode = double(NS5.Data(2, :));
    photodiode = filtfilt(b, a, photodiode);

    thresh_diff = 500;
    window_ms = round(20 * original_freq / 1000);

    rising_idx = [];
    for i = (window_ms+1):length(photodiode)
        baseline = mean(photodiode(i-window_ms:i-1));

        if (photodiode(i) - baseline) > thresh_diff
            if isempty(rising_idx) || i - rising_idx(end) > window_ms
                rising_idx(end+1) = i;
            end
        end
    end

    falling_idx = nan(size(rising_idx));

    for r = 1:length(rising_idx)
        start_i = rising_idx(r);
        base = mean(photodiode(max(1,start_i-window_ms):start_i));
        drop_thresh = base + thresh_diff/4;

        j = start_i + 1;

        while j <= length(photodiode) && photodiode(j) > drop_thresh
            j = j + 1;
        end

        if j <= length(photodiode)
            falling_idx(r) = j;
        end
    end

    fallRiseBadGap = round(50 * original_freq / 1000);

    valid_mask = (falling_idx - rising_idx) >= fallRiseBadGap;

    rising_idx = rising_idx(valid_mask);
    falling_idx = falling_idx(valid_mask);

    isi = diff(rising_idx);

    doublet_min = round(50 * original_freq / 1000);
    doublet_max = round(200 * original_freq / 1000);

    doublet_idx = find(doublet_min < isi & isi <= doublet_max);
    doublet_rising_idx = rising_idx(doublet_idx);

    n_between_list = nan(1, length(doublet_idx)-1);

    for k = 1:length(doublet_idx)-1
        end_curr = rising_idx(doublet_idx(k)+1);
        start_next = rising_idx(doublet_idx(k+1));

        between_edges = rising_idx(rising_idx > end_curr & rising_idx < start_next);

        n_between_list(k) = numel(between_edges);
    end

    trials_good = 5;
    trials_good_miss = 3;
    trials_message_show_and_space_overlaped = 4;
    trials_message_show_and_space_overlaped_miss = 2;
    consecutive = 0;

    bad_idx = find(n_between_list ~= trials_good & ...
                   n_between_list ~= trials_good_miss & ...
                   n_between_list ~= trials_message_show_and_space_overlaped & ...
                   n_between_list ~= trials_message_show_and_space_overlaped_miss & ...
                   n_between_list ~= consecutive);

    doublet_rising_idx(bad_idx+1) = [];

    n_between_doublets = nan(1, length(doublet_rising_idx));
    spike_times_between = cell(1, length(doublet_rising_idx));

    for k = 1:length(doublet_rising_idx)

        pos = find(rising_idx == doublet_rising_idx(k), 1, 'first');
        end_curr = rising_idx(pos+1);

        if k < length(doublet_rising_idx)
            start_next = doublet_rising_idx(k+1);
        else
            start_next = length(photodiode);
        end

        between_edges = rising_idx(rising_idx > end_curr & rising_idx < start_next);

        n_between_doublets(k) = numel(between_edges);
        spike_times_between{k} = between_edges;
    end

    all_trigs = [];

    for k = 1:length(doublet_rising_idx)

        this_trigs = [];
        this_times = [];

        this_trigs(end+1) = 1;
        this_times(end+1) = double(doublet_rising_idx(k));

        if n_between_doublets(k) >= 5

            codes = [2 3 4 5 6];
            times = double(spike_times_between{k}(1:5));

        elseif n_between_doublets(k) == 4

            codes = [2 3 4 5 6];
            t = double(spike_times_between{k});
            times = [t(1:2), t(2)+round(30 * original_freq / 1000), t(3:4)];

        elseif n_between_doublets(k) == 3

            codes = [2 3 4];
            t = double(spike_times_between{k});
            times = t(1:3);

        elseif n_between_doublets(k) == 2

            codes = [2 3 4];
            t = double(spike_times_between{k});
            times = [t(1), t(2), t(2)+round(30 * original_freq / 1000)];

        else

            codes = [];
            times = [];

        end

        this_trigs = [this_trigs codes];
        this_times = [this_times times];

        all_trigs = [all_trigs [this_trigs; this_times]];
    end

    trigs = all_trigs;

    save(fullfile(output_folder, [ptID '_trigs.mat']), 'trigs');

end