%% task analysis

save_fig_path = 'C:\Users\amarica\Documents\Lab stuff\Random figs\Long_str_stuff';
task_save_fig_path = fullfile(save_fig_path, 'Task');

if ~isfolder(task_save_fig_path)
    mkdir(task_save_fig_path)
end

animal = 'AM026';
workflow = {'stim_wheel_right*'};
recordings = plab.find_recordings(animal, [], workflow);

struct_ephys_data = struct;
missing_ephys = 0;

% for AM021 from 05/04 - end (day 3 - end)
for use_rec=3:length(recordings)
    % use_rec = 1;
%         rec_day = '2023-07-03';
%         use_rec = strcmp(rec_day,{recordings.day});
    % use_rec = length(recordings)-1;

    rec_day = recordings(use_rec).day;
    rec_time = recordings(use_rec).recording{end};

    verbose = true;

    load_parts.behavior = true;
    load_parts.ephys = true; 

    ap.load_recording

    ephys_path = dir(plab.locations.filename('server', animal,rec_day,[],'ephys'));
    if isempty(ephys_path)
        disp(['Day ' rec_day ' does not have ephys'])
        missing_ephys = missing_ephys  + 1;
        continue
    end

    good_templates_idx = find(good_templates);


    if contains(recordings(use_rec).workflow{end}, 'stim_wheel_right')
        trial_stim_values = 90*ones(length(stimOn_times), 1);
    end

    % plot

    %      ap.cellraster(stimOn_times,trial_stim_values);

    % get rewarded trials
    rewarded_trials = logical([trial_events.values.Outcome]);

    % get good reward times
    stimOff_times = photodiode_times(photodiode_values==0);
    [~, good_reward_idx] = min(abs(reward_times - stimOff_times(rewarded_trials)'));
    good_reward_times = reward_times(good_reward_idx);

    % get rewarded completed trials
    completed_trials = ~cellfun(@isempty,{trial_events.values.Outcome});
    completed_rewarded_trials = rewarded_trials;
    completed_rewarded_trials(find(completed_trials==0)) = false;

    %% Get striatum boundaries

%     %%% Get correlation of MUA in sliding sindows
%     depth_corr_window = 100; % MUA window in microns
%     depth_corr_window_spacing = 50; % MUA window spacing in microns
% 
%     max_depths = 3840; % (hardcode, sometimes kilosort2 drops channels)
% 
%     depth_corr_bins = [0:depth_corr_window_spacing:(max_depths-depth_corr_window); ...
%         (0:depth_corr_window_spacing:(max_depths-depth_corr_window))+depth_corr_window];
%     depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;
% 
%     spike_binning_t = 0.01; % seconds
%     spike_binning_t_edges = nanmin(spike_times_timelite):spike_binning_t:nanmax(spike_times_timelite);
% 
%     binned_spikes_depth = zeros(size(depth_corr_bins,2),length(spike_binning_t_edges)-1);
%     for curr_depth = 1:size(depth_corr_bins,2)
%         curr_depth_templates_idx = ...
%             find(template_depths >= depth_corr_bins(1,curr_depth) & ...
%             template_depths < depth_corr_bins(2,curr_depth));
% 
%         binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timelite( ...
%             ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
%     end
% 
%     mua_corr = corrcoef(binned_spikes_depth');
% 
%     figure;
%     imagesc(mua_corr)
%     %%% Estimate start and end depths of striatum
% 
%     % % end of striatum: biggest (smoothed) drop in MUA correlation near end
%     % groups_back = 30;
%     % mua_corr_end = medfilt2(mua_corr(end-groups_back+1:end,end-groups_back+1:end),[3,3]);
%     % mua_corr_end(triu(true(length(mua_corr_end)),0)) = nan;
%     % median_corr = medfilt1(nanmedian(mua_corr_end,2),3);
%     % [x,max_corr_drop] = min(diff(median_corr));
%     % str_end = depth_corr_bin_centers(end-groups_back+max_corr_drop);
% 
%     % (new method)
%     % end of striatum: minimum correlation on dim 1 * dim 2
%     % (to look for the biggest dead space between correlated blocks)
%     groups_back = 20;
%     mua_corr_end = medfilt2(mua_corr(end-groups_back+1:end,end-groups_back+1:end),[3,3]);
%     mua_corr_end(triu(true(length(mua_corr_end)),0)) = nan;
%     mean_corr_dim1 = nanmean(mua_corr_end,2);
%     mean_corr_dim2 = nanmean(mua_corr_end,1);
%     mean_corr_mult = mean_corr_dim1.*mean_corr_dim2';
%     [~,mean_corr_mult_min_idx] = min(mean_corr_mult);
%     str_end = depth_corr_bin_centers(end-groups_back + mean_corr_mult_min_idx - 2); % err early: back up 2 (100 um)
% 
%     % start of striatum: look for ventricle
%     % (by biggest gap between templates)
%     min_gap = 200;
%     sorted_good_template_depths_for_gap = sort([0;template_depths]);
%     %     [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
%     [sort_gap, sort_gap_idx] = sort(diff(sorted_good_template_depths_for_gap));
% 
%     if sort_gap(end) > min_gap
%         str_start = sorted_good_template_depths_for_gap(sort_gap_idx(end)+1)-1;
%         template_str_start = sort_gap_idx(end)+1;
%         if str_start > str_end && sort_gap(end-1) > min_gap
%             str_start = sorted_good_template_depths_for_gap(sort_gap_idx(end-1)+1)-1;
%         end
%     else
%         str_start = sorted_good_template_depths_for_gap(2);
%     end
% 
%     % define str depth
%     str_depth = [str_start,str_end];

        
        [sorted_template_depths, sorted_template_depth_idx] = sort(template_depths);

        % find str start where unit depth distr not linear anymore
        idx_str_start = ischange(sorted_template_depths, 'linear','MaxNumChanges',1);

%         test_x_axis = 1:length(sorted_template_depths);
%         figure;
%         plot(sorted_template_depths, 'o')
%         hold on
%         plot(test_x_axis(test_a), sorted_template_depths(test_a), '*')
%         title(rec_day)


        str_start = sorted_template_depths(idx_str_start);
        str_end = sorted_template_depths(end);
        str_depth = [str_start,str_end];

        % plot
        figure;
        unit_axes = nexttile;
        set(unit_axes,'YDir','reverse');
        hold on;

        norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
        unit_dots = scatter3(norm_spike_n,template_depths(unique(spike_templates)), ...
            unique(spike_templates),20,'k','filled');


        yline(str_depth, 'red')
        xlim(unit_axes,[-0.1,1]);
        ylim([-50, max(template_depths)+50]);
        ylabel('Depth (\mum)')
        xlabel('Normalized log rate')
        title(rec_day)

    %% classify str neuron type

%     % Split striatal/nonstriatal cells
%     str_good_templates = template_depths >= str_depth(1) & template_depths <= str_depth(2);
%     %     non_str_templates = ~str_templates;
% 
%     % Define the window to look for spiking statistics in (spikes go in and
%     % out, so take the bin with the largest firing rate for each cell and work
%     % with that one)
%     % spiking_stat_window = 60*5; % seconds
%     % spiking_stat_bins = min(spike_times_timelite):spiking_stat_window: ...
%     %     max(spike_times_timelite);
% 
%     % % (for whole session)
%     spiking_stat_window = max(spike_times_timelite)-min(spike_times_timelite);
%     spiking_stat_bins = [min(spike_times_timelite),max(spike_times_timelite)];
% 
%     % Get firing rate across the session
%     bin_spikes = nan(size(templates,1), ...
%         length(spiking_stat_bins)-1);
%     for curr_template_idx = 1:length(good_templates_idx)
%         bin_spikes(good_templates_idx(curr_template_idx),:) = ...
%             histcounts(spike_times_timelite(spike_templates == good_templates_idx(curr_template_idx)), ...
%             spiking_stat_bins);
%     end
%     min_spikes = 10;
%     use_spiking_stat_bins = bsxfun(@ge,bin_spikes,prctile(bin_spikes,80,2)) & bin_spikes > min_spikes;
%     spike_rate = sum(bin_spikes.*use_spiking_stat_bins,2)./ ...
%         (sum(use_spiking_stat_bins,2)*spiking_stat_window);
% 
%     % Get proportion of ISI > 2s (Yamin/Cohen 2013) and CV2 (Stalnaker/Schoenbaum 2016)
%     prop_long_isi = nan(size(templates,1),1);
%     cv2 = nan(size(templates,1),1);
%     for curr_template_idx = 1:length(good_templates_idx)
% 
%         long_isi_total = 0;
%         isi_ratios = [];
%         for curr_bin = find(use_spiking_stat_bins(good_templates_idx(curr_template_idx),:))
%             curr_spike_times = spike_times_timelite( ...
%                 spike_times_timelite > spiking_stat_bins(curr_bin) & ...
%                 spike_times_timelite < spiking_stat_bins(curr_bin+1) & ...
%                 spike_templates == good_templates_idx(curr_template_idx));
%             curr_isi = diff(curr_spike_times);
% 
%             long_isi_total = long_isi_total + sum(curr_isi(curr_isi > 2));
% 
%             isi_ratios = [isi_ratios;(2*abs(curr_isi(2:end) - curr_isi(1:end-1)))./ ...
%                 (curr_isi(2:end) + curr_isi(1:end-1))];
%         end
% 
%         prop_long_isi(curr_template_idx) = long_isi_total/ ...
%             (sum(use_spiking_stat_bins(good_templates_idx(curr_template_idx),:))*spiking_stat_window);
%         cv2(good_templates_idx(curr_template_idx)) = nanmean(isi_ratios);
% 
%     end
% 
% 
%     % Cortical classification (like Bartho JNeurophys 2004)
%     waveform_duration_cutoff = 400;
%     %     narrow = non_str_templates & templateDuration_us <= waveform_duration_cutoff;
%     %     wide = non_str_templates & templateDuration_us > waveform_duration_cutoff;
% 
%     % Striatum classification
%     prop_long_isi_cutoff = 0.35;
%     cv2_cutoff = 0.8;
% 
%     msn = str_good_templates & ...
%         templateDuration_us > waveform_duration_cutoff & ...
%         prop_long_isi >= prop_long_isi_cutoff;
% 
%     fsi = str_good_templates & ...
%         templateDuration_us <= waveform_duration_cutoff & ...
%         prop_long_isi < prop_long_isi_cutoff;
% 
%     tan = str_good_templates & ...
%         templateDuration_us > waveform_duration_cutoff & ...
%         prop_long_isi < prop_long_isi_cutoff;
% 
%     uin = str_good_templates & ~msn & ~fsi & ~tan;
% 
%     disp(['Found ' num2str(sum(msn)) ' MSNs' newline ...
%         'Found ' num2str(sum(fsi)) ' FSIs' newline ...
%         'Found ' num2str(sum(tan)) ' TANs' newline ...
%         'Found ' num2str(sum(uin)) ' UINs'])
% 
%     waveform_t = 1e3*((0:size(templates,2)-1)/ephys_sample_rate);
% 
% %         figure; hold on;
% %         p = plot(waveform_t,waveforms(str_good_templates,:)');
% %         set(p(msn(str_good_templates)),'color','m')
% %         set(p(fsi(str_good_templates)),'color','b')
% %         set(p(tan(str_good_templates)),'color','g')
% %         set(p(uin(str_good_templates)),'color','c')
% %         xlabel('Time (ms)')
% %         title('Striatum');
% %         legend([p(find(msn(str_good_templates),1)),p(find(fsi(str_good_templates),1)), ...
% %             p(find(tan(str_good_templates),1)),p(find(uin(str_good_templates),1))],{'MSN','FSI','TAN','UIN'});


    %% psth
    %% - multiunit

    % create time vector around all stim onsets
    bin_window = 0.1;
    bin_edges = -0.5:bin_window:2;
    around_stim_time = stimOn_times(completed_rewarded_trials) + bin_edges;

    %% -- ?? striatum
    str_spikes = spike_depths>str_depth(1) & spike_depths<str_depth(2);
    str_spike_times = spike_times_timelite(str_spikes);

    % spike counts binned for each stim
    str_spikes_in_stim_time = nan(length(find(completed_rewarded_trials)), size(around_stim_time, 2)-1);

    % transpose to get the right shape
    str_spikes_in_stim_time = cell2mat(arrayfun(@(trial_id) histcounts(str_spike_times, around_stim_time(trial_id,:))', ...
        1:size(around_stim_time, 1), 'UniformOutput',false))' / bin_window;

    %% spike data
    % spike counts binned for each stim

    % create time vector around all stim onsets
    bin_window = 0.001;
    bin_edges = -0.5:bin_window:2;
    around_stim_time = stimOn_times(completed_rewarded_trials) + bin_edges;

    % for plot
    bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;

    % task stim
    task_all_spikes_in_stim_time = nan(length(good_templates_idx), length(find(completed_rewarded_trials)), size(around_stim_time, 2)-1);

    for unit_idx=1:length(good_templates_idx)

        unit_spikes = spike_templates == unit_idx;

        unit_spike_times = spike_times_timelite(unit_spikes);

        task_all_spikes_in_stim_time(unit_idx, :, :) = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, around_stim_time(trial_id,:))', ...
            1:size(around_stim_time, 1), 'UniformOutput',false))' / bin_window;
    end

    % spike counts binned for each move after stim
    % create time vector around all move after stim onsets
    bin_window = 0.001;
    bin_edges = -0.5:bin_window:2;
    around_stim_move_time = stim_move_time(completed_rewarded_trials) + bin_edges;

    % for plot
    bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;

    % task stim
    task_all_spikes_in_stim_move_time = nan(length(good_templates_idx), length(find(completed_rewarded_trials)), size(around_stim_move_time, 2)-1);

    for unit_idx=1:length(good_templates_idx)

        unit_spikes = spike_templates == unit_idx;

        unit_spike_times = spike_times_timelite(unit_spikes);

        task_all_spikes_in_stim_move_time(unit_idx, :, :) = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, around_stim_move_time(trial_id,:))', ...
            1:size(around_stim_move_time, 1), 'UniformOutput',false))' / bin_window;
    end


    % spike counts binned for good reward
    % create time vector around reward onsets
    bin_window = 0.001;
    bin_edges = -0.5:bin_window:2;
    around_reward_time = good_reward_times + bin_edges;

    % for plot
    bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;

    % task stim
    task_all_spikes_in_reward_time = nan(length(good_templates_idx), length(good_reward_times), size(around_reward_time, 2)-1);

    for unit_idx=1:length(good_templates_idx)

        unit_spikes = spike_templates == unit_idx;

        unit_spike_times = spike_times_timelite(unit_spikes);

        task_all_spikes_in_reward_time(unit_idx, :, :) = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, around_reward_time(trial_id,:))', ...
            1:size(around_reward_time, 1), 'UniformOutput',false))' / bin_window;
    end


    %% mean

    % stim
    mean_all_spikes_in_stim_time = squeeze(mean(task_all_spikes_in_stim_time, 2));

    % move
    mean_all_spikes_in_stim_move_time = squeeze(mean(task_all_spikes_in_stim_move_time, 2));

    % reward
    mean_all_spikes_in_reward_time = squeeze(mean(task_all_spikes_in_reward_time, 2));

    %% - smooth
    % define gaussian window
    gauss_win = gausswin(51, 3)';

    % stim
    stim_smooth_all_units = filter(gauss_win,sum(gauss_win),mean_all_spikes_in_stim_time, [], 2);

    % stim move
    stim_move_smooth_all_units = filter(gauss_win,sum(gauss_win),mean_all_spikes_in_stim_move_time, [], 2);

    % reward
    reward_smooth_all_units = filter(gauss_win,sum(gauss_win),mean_all_spikes_in_reward_time, [], 2);


    %% - get baseline

    % stim
    stim_smooth_baseline = mean(stim_smooth_all_units(:,bin_centres>-0.2 & bin_centres<0), 2);

    % stim move
    stim_move_smooth_baseline = mean(stim_move_smooth_all_units(:,bin_centres>-0.2 & bin_centres<0), 2);

    % reward
    reward_smooth_baseline = mean(reward_smooth_all_units(:,bin_centres>-0.2 & bin_centres<0), 2);


    %% - normalize

    % stim
    stim_norm_smooth_all_units = (stim_smooth_all_units - stim_smooth_baseline) ...
        ./ (stim_smooth_baseline + std(stim_smooth_baseline));

    % stim move
    stim_move_norm_smooth_all_units = (stim_move_smooth_all_units - stim_move_smooth_baseline) ...
        ./ (stim_move_smooth_baseline + std(stim_move_smooth_baseline));

    % reward
    reward_norm_smooth_all_units = (reward_smooth_all_units - reward_smooth_baseline) ...
        ./ (reward_smooth_baseline + std(reward_smooth_baseline));


    %% - split probe into PSTHs

    %%%%% something is wrong with indexing
    %%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [sorted_template_depths, sorted_template_depth_idx] = sort(template_depths); 

    % define str split
    str_split_step = floor((str_end-str_start)/3);
    str_split = [str_start ...
        str_start+str_split_step ...
        str_start+2*str_split_step ...
        str_end];

    % find str split and get psths
%     a = discretize(template_depths, str_split)
%     b = find(a==1)
%     ismember(b,str_1_templates)

    template_str_1_idx = [find(sorted_template_depths>str_split(1) & sorted_template_depths<=str_split(2))];
    str_1_templates = sorted_template_depth_idx(template_str_1_idx);

    template_str_2_idx = [find(sorted_template_depths>str_split(2) & sorted_template_depths<=str_split(3))];
    str_2_templates = sorted_template_depth_idx(template_str_2_idx);

    template_str_3_idx = [find(sorted_template_depths>str_split(3) & sorted_template_depths<=str_split(4))];
    str_3_templates = sorted_template_depth_idx(template_str_3_idx);

% 
%     a = stim_norm_smooth_all_units(str_1_templates,:)
% 
%     b = stim_norm_smooth_all_units(297,:)

    % psth stim 
    psth_str_1_stim = mean(stim_norm_smooth_all_units(str_1_templates,:), 1);
    psth_str_2_stim = mean(stim_norm_smooth_all_units(str_2_templates,:), 1);
    psth_str_3_stim = mean(stim_norm_smooth_all_units(str_3_templates,:), 1);

    % psth move
    psth_str_1_stim_move = mean(stim_move_norm_smooth_all_units(str_1_templates,:), 1);
    psth_str_2_stim_move = mean(stim_move_norm_smooth_all_units(str_2_templates,:), 1);
    psth_str_3_stim_move = mean(stim_move_norm_smooth_all_units(str_3_templates,:), 1);

    % psth reward
    psth_str_1_reward = mean(reward_norm_smooth_all_units(str_1_templates,:), 1);
    psth_str_2_reward = mean(reward_norm_smooth_all_units(str_2_templates,:), 1);
    psth_str_3_reward = mean(reward_norm_smooth_all_units(str_3_templates,:), 1);

    % get psth for all str and save out
    psths_str_all_stim{use_rec} = mean([psth_str_1_stim; psth_str_2_stim; psth_str_3_stim], 1);
    psths_str_all_stim_move{use_rec} = mean([psth_str_1_stim_move; psth_str_2_stim_move; psth_str_3_stim_move], 1);
    psths_str_all_reward{use_rec} = mean([psth_str_1_reward; psth_str_2_reward; psth_str_3_reward], 1);

    % cortex
    cortex_templates = sorted_template_depth_idx(1:template_str_1_idx(1));
    psth_cortex_stim = mean(stim_norm_smooth_all_units(cortex_templates,:), 1);
    psth_cortex_stim_move = mean(stim_move_norm_smooth_all_units(cortex_templates,:), 1);
    psth_cortex_reward = mean(reward_norm_smooth_all_units(cortex_templates,:), 1);


    % under striatum
    if isempty(template_str_3_idx)
        psth_under_str_stim = nan(1, length(psth_str_1_stim));
        psth_under_str_stim_move = nan(1, length(psth_str_1_stim_move));
        psth_under_str_reward = nan(1, length(psth_str_1_reward));
    else
        under_str_templates = sorted_template_depth_idx(template_str_3_idx(end)+1:end, 1);
        psth_under_str_stim = mean(stim_norm_smooth_all_units(under_str_templates,:), 1);
        psth_under_str_stim_move = mean(stim_move_norm_smooth_all_units(under_str_templates,:), 1);
        psth_under_str_reward = mean(reward_norm_smooth_all_units(under_str_templates,:), 1);

    end

    % append all together
    all_psths_stim = vertcat(psth_cortex_stim, ...
        psth_str_1_stim, ...
        psth_str_2_stim, ...
        psth_str_3_stim, ...
        psth_under_str_stim);

    all_psths_stim_move = vertcat(psth_cortex_stim_move, ...
        psth_str_1_stim_move, ...
        psth_str_2_stim_move, ...
        psth_str_3_stim_move, ...
        psth_under_str_stim_move);

    all_psths_reward = vertcat(psth_cortex_reward, ...
        psth_str_1_reward, ...
        psth_str_2_reward, ...
        psth_str_3_reward, ...
        psth_under_str_reward);


    % plot to check
    figure;
    sgtitle(['Task ' rec_day]);

    stim = subplot(1,3,1);
    title('Stim')
    plot(psth_str_1_stim)
    hold on;
    plot(psth_str_2_stim)
    hold on;
    plot(psth_str_3_stim)
    legend({'1', '2', '3'})

    stim_move = subplot(1,3,2);
    title('Stim move')
    plot(psth_str_1_stim_move)
    hold on;
    plot(psth_str_2_stim_move)
    hold on;
    plot(psth_str_3_stim_move)
    legend({'1', '2', '3'})

    rew = subplot(1,3,3);
    title('Reward')
    plot(psth_str_1_reward)
    hold on;
    plot(psth_str_2_reward)
    hold on;
    plot(psth_str_3_reward)
    legend({'1', '2', '3'})

    % make plot of depth 
    % plot
%     figure;
%     unit_axes = nexttile;
%     set(unit_axes,'YDir','reverse');
%     hold on;
% 
%     norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
%     unit_dots = scatter3(norm_spike_n,template_depths(unique(spike_templates)), ...
%         unique(spike_templates),20,'k','filled');
% 
%     str_1_unit_dots = scatter3(norm_spike_n(str_1_templates),template_depths(str_1_templates), ...
%             str_1_templates,20,'green','filled');
% 
%     str_2_unit_dots = scatter3(norm_spike_n(str_2_templates),template_depths(str_2_templates), ...
%             str_2_templates,20,'blue','filled');
% 
%     str_3_unit_dots = scatter3(norm_spike_n(str_3_templates),template_depths(str_3_templates), ...
%             str_3_templates,20,'magenta','filled');
% 
%     yline(str_split, 'red')
%     xlim(unit_axes,[-0.1,1]);
%     ylim([-50, max(template_depths)+50]);
%     ylabel('Depth (\mum)')
%     xlabel('Normalized log rate')
%     title(rec_day)


    %% plots per day

    task_fig = figure('Position', get(0, 'Screensize')) %, Visible="off");
    sgtitle(['Task ' rec_day]);

    %% - heatmaps
    % -- stim
    upper_caxis = max(abs(stim_norm_smooth_all_units), [],"all");
    stim = subplot(2,3,1);
    imagesc(bin_centres, [], stim_norm_smooth_all_units(sorted_template_depth_idx, :));
    xline(0, 'LineWidth', 1);
    yline([find(sorted_template_depth_idx==str_1_templates(1)) ...
        find(sorted_template_depth_idx==str_2_templates(1)) ...
        find(sorted_template_depth_idx==str_3_templates(1)) ...
        find(sorted_template_depth_idx==str_3_templates(end))], 'green');
    colormap(stim, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;
    ylabel('Heatmaps', 'FontWeight', 'bold', 'FontSize', 14);
    title('Stim', 'FontWeight', 'bold', 'FontSize', 14);

    % -- stim move
    upper_caxis = max(abs(stim_move_norm_smooth_all_units), [],"all");
    move = subplot(2,3,2);
    imagesc(bin_centres, [], stim_move_norm_smooth_all_units(sorted_template_depth_idx, :));
    xline(0, 'LineWidth', 1);
        yline([find(sorted_template_depth_idx==str_1_templates(1)) ...
        find(sorted_template_depth_idx==str_2_templates(1)) ...
        find(sorted_template_depth_idx==str_3_templates(1)) ...
        find(sorted_template_depth_idx==str_3_templates(end))], 'green');
    colormap(move, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;
    title('Move','FontWeight', 'bold', 'FontSize', 14);

    % -- reward
    upper_caxis = max(abs(reward_norm_smooth_all_units), [],"all");
    rew = subplot(2,3,3);
    imagesc(bin_centres, [], reward_norm_smooth_all_units(sorted_template_depth_idx, :));
    xline(0, 'LineWidth', 1);
        yline([find(sorted_template_depth_idx==str_1_templates(1)) ...
        find(sorted_template_depth_idx==str_2_templates(1)) ...
        find(sorted_template_depth_idx==str_3_templates(1)) ...
        find(sorted_template_depth_idx==str_3_templates(end))], 'green');
    colormap(rew, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;
    title('Reward', 'FontWeight', 'bold', 'FontSize', 14);

    %% - PSTHs

    % -- stim
    upper_caxis = max(abs(stim_norm_smooth_all_units), [],"all");
    stim = subplot(2,3,4);
    imagesc(bin_centres, [], all_psths_stim);
    yticks(1:size(all_psths_stim,1));
    yticklabels({'Cortex', 'Str 1', 'Str 2', 'Str 3', 'Under str'});
    xline(0, 'LineWidth', 1);
    colormap(stim, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;
    ylabel('PSTHs', 'FontWeight', 'bold', 'FontSize', 14);

    % -- stim move
    upper_caxis = max(abs(stim_move_norm_smooth_all_units), [],"all");
    move = subplot(2,3,5);
    imagesc(bin_centres, [], all_psths_stim_move);
    yticks(1:size(all_psths_stim,1));
    yticklabels({'Cortex', 'Str 1', 'Str 2', 'Str 3', 'Under str'});
    xline(0, 'LineWidth', 1);
    colormap(move, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;

    % -- reward
    upper_caxis = max(abs(reward_norm_smooth_all_units), [],"all");
    rew = subplot(2,3,6);
    imagesc(bin_centres, [], all_psths_reward);
    yticks(1:size(all_psths_stim,1));
    yticklabels({'Cortex', 'Str 1', 'Str 2', 'Str 3', 'Under str'});
    xline(0, 'LineWidth', 1);
    colormap(rew, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;

    % save this fig
    task_fig_name = [animal '_' rec_day '_task.tif'];
    task_fig_path = fullfile(task_save_fig_path,task_fig_name);
    saveas(task_fig, task_fig_path);


        %% save data in big struct
        struct_ephys_data(use_rec-missing_ephys).rec_day = rec_day;
    
        struct_ephys_data(use_rec-missing_ephys).good_templates = good_templates;
        struct_ephys_data(use_rec-missing_ephys).template_depths = template_depths;
        struct_ephys_data(use_rec-missing_ephys).spike_templates = spike_templates;
    
        struct_ephys_data(use_rec-missing_ephys).spike_times_timelite = spike_times_timelite; 

        struct_ephys_data(use_rec-missing_ephys).stim_to_move = stim_to_move;
        struct_ephys_data(use_rec-missing_ephys).stimOn_times = stimOn_times;
        struct_ephys_data(use_rec-missing_ephys).completed_rewarded_trials = completed_rewarded_trials;

%         struct_ephys_data(use_rec-missing_ephys).msn = msn;
%         struct_ephys_data(use_rec-missing_ephys).fsi = fsi;
%         struct_ephys_data(use_rec-missing_ephys).tan = tan;
%         struct_ephys_data(use_rec-missing_ephys).uin = uin;
    
        struct_ephys_data(use_rec-missing_ephys).bin_edges = bin_edges;
        struct_ephys_data(use_rec-missing_ephys).bin_centres = bin_centres;
        
        struct_ephys_data(use_rec-missing_ephys).str_depth = str_depth;
        struct_ephys_data(use_rec-missing_ephys).str_spikes_in_stim_time = str_spikes_in_stim_time;

        % baseline
        struct_ephys_data(use_rec-missing_ephys).stim_smooth_baseline = stim_smooth_baseline;
        struct_ephys_data(use_rec-missing_ephys).stim_move_smooth_baseline = stim_move_smooth_baseline;
        struct_ephys_data(use_rec-missing_ephys).reward_smooth_baseline = reward_smooth_baseline;

        % norm aligned traces
        struct_ephys_data(use_rec-missing_ephys).stim_norm_smooth_all_units = stim_norm_smooth_all_units;
        struct_ephys_data(use_rec-missing_ephys).stim_move_norm_smooth_all_units = stim_move_norm_smooth_all_units;
        struct_ephys_data(use_rec-missing_ephys).reward_norm_smooth_all_units = reward_norm_smooth_all_units;

        % str templates
        struct_ephys_data(use_rec-missing_ephys).str_1_templates = str_1_templates;
        struct_ephys_data(use_rec-missing_ephys).str_2_templates = str_2_templates;
        struct_ephys_data(use_rec-missing_ephys).str_3_templates = str_3_templates;
        struct_ephys_data(use_rec-missing_ephys).cortex_templates = cortex_templates;

        % stim
        struct_ephys_data(use_rec-missing_ephys).psth_cortex_stim = psth_cortex_stim;
        struct_ephys_data(use_rec-missing_ephys).psth_str_1_stim = psth_str_1_stim;
        struct_ephys_data(use_rec-missing_ephys).psth_str_2_stim = psth_str_2_stim;
        struct_ephys_data(use_rec-missing_ephys).psth_str_3_stim = psth_str_3_stim;
        struct_ephys_data(use_rec-missing_ephys).psth_under_str_stim = psth_under_str_stim;

        % stim move
        struct_ephys_data(use_rec-missing_ephys).psth_cortex_stim_move = psth_cortex_stim_move;
        struct_ephys_data(use_rec-missing_ephys).psth_str_1_stim_move = psth_str_1_stim_move;
        struct_ephys_data(use_rec-missing_ephys).psth_str_2_stim_move = psth_str_2_stim_move;
        struct_ephys_data(use_rec-missing_ephys).psth_str_3_stim_move = psth_str_3_stim_move;
        struct_ephys_data(use_rec-missing_ephys).psth_under_str_stim_move = psth_under_str_stim_move;

        % reward
        struct_ephys_data(use_rec-missing_ephys).psth_cortex_reward = psth_cortex_reward;
        struct_ephys_data(use_rec-missing_ephys).psth_str_1_reward = psth_str_1_reward;
        struct_ephys_data(use_rec-missing_ephys).psth_str_2_reward = psth_str_2_reward;
        struct_ephys_data(use_rec-missing_ephys).psth_str_3_reward = psth_str_3_reward;
        struct_ephys_data(use_rec-missing_ephys).psth_under_str_reward = psth_under_str_reward;

end

task_ephys_data = struct2table(struct_ephys_data);  %convert structure to table
% t_ephys_data = rmmissing(t_ephys_data);
task_ephys_data = task_ephys_data(~cellfun(@isempty, task_ephys_data.(1)), :);


save_name = [animal '_task_str_ephys_data'];
save(save_name, "task_ephys_data", "-v7.3");

clear all;

%% load 

save_fig_path = 'C:\Users\amarica\Documents\Lab stuff\Random figs\Long_str_stuff\Task';

animal = 'AM021';
load([animal '_task_str_ephys_data']);

%% plot psth across days

n_days = height(task_ephys_data);
bin_centres = task_ephys_data.bin_centres{end};


psth_across_days_fig = figure('Position', get(0, 'Screensize'));
sgtitle([animal ' PSTHs across days']);

% stim
cortex_stim = subplot(3, 5, 1);
plot(bin_centres, vertcat(task_ephys_data.psth_cortex_stim{:})');
cortex_stim.ColorOrder = brewermap(n_days,"PuRd");
title('Cortex')
ylabel('Stim', 'FontWeight','bold');
str_1_stim = subplot(3, 5, 2);
plot(bin_centres, vertcat(task_ephys_data.psth_str_1_stim{:})');
str_1_stim.ColorOrder = brewermap(n_days,"PuRd");
title('Str 1')
str_2_stim = subplot(3, 5, 3);
plot(bin_centres, vertcat(task_ephys_data.psth_str_2_stim{:})');
str_2_stim.ColorOrder = brewermap(n_days,"PuRd");
title('Str 2')
str_3_stim = subplot(3, 5, 4);
plot(bin_centres, vertcat(task_ephys_data.psth_str_3_stim{:})');
str_3_stim.ColorOrder = brewermap(n_days,"PuRd");
title('Str 3')
under_str_stim = subplot(3, 5, 5);
plot(bin_centres, vertcat(task_ephys_data.psth_under_str_stim{:})');
under_str_stim.ColorOrder = brewermap(n_days,"PuRd");
title('Under str')

% stim move
cortex_stim_move = subplot(3, 5, 6);
plot(bin_centres, vertcat(task_ephys_data.psth_cortex_stim_move{:})');
cortex_stim_move.ColorOrder = brewermap(n_days,"Greens");
ylabel('Move', 'FontWeight','bold');
str_1_stim_move = subplot(3, 5, 7);
plot(bin_centres, vertcat(task_ephys_data.psth_str_1_stim_move{:})');
str_1_stim_move.ColorOrder = brewermap(n_days,"Greens");
str_2_stim_move = subplot(3, 5, 8);
plot(bin_centres, vertcat(task_ephys_data.psth_str_2_stim_move{:})');
str_2_stim_move.ColorOrder = brewermap(n_days,"Greens");
str_3_stim_move = subplot(3, 5, 9);
plot(bin_centres, vertcat(task_ephys_data.psth_str_3_stim_move{:})');
str_3_stim_move.ColorOrder = brewermap(n_days,"Greens");
under_str_stim_move = subplot(3, 5, 10);
plot(bin_centres, vertcat(task_ephys_data.psth_under_str_stim_move{:})');
under_str_stim_move.ColorOrder = brewermap(n_days,"Greens");

% reward
cortex_reward = subplot(3, 5, 11);
plot(bin_centres, vertcat(task_ephys_data.psth_cortex_reward{:})');
cortex_reward.ColorOrder = brewermap(n_days,"Blues");
ylabel('Reward', 'FontWeight','bold');
str_1_reward = subplot(3, 5, 12);
plot(bin_centres, vertcat(task_ephys_data.psth_str_1_reward{:})');
str_1_reward.ColorOrder = brewermap(n_days,"Blues");
str_2_reward = subplot(3, 5, 13);
plot(bin_centres, vertcat(task_ephys_data.psth_str_2_reward{:})');
str_2_reward.ColorOrder = brewermap(n_days,"Blues");
str_3_reward = subplot(3, 5, 14);
plot(bin_centres, vertcat(task_ephys_data.psth_str_3_reward{:})');
str_3_reward.ColorOrder = brewermap(n_days,"Blues");
under_str_reward = subplot(3, 5, 15);
plot(bin_centres, vertcat(task_ephys_data.psth_under_str_reward{:})');
under_str_reward.ColorOrder = brewermap(n_days,"Blues");

psth_across_days_fig_name = [animal '_psths_across_days.tif'];
psth_across_days_fig_path = fullfile(save_fig_path, psth_across_days_fig_name);
saveas(psth_across_days_fig, psth_across_days_fig_path);

%% depth plots
str_split_fig = figure('Position', get(0, 'Screensize'));
tiledlayout('flow');
sgtitle([animal ' Str split per depth'])
for use_rec=1:height(task_ephys_data)
    unit_axes = nexttile;
    set(unit_axes,'YDir','reverse');
    hold on;

    norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
    unit_dots = scatter3(norm_spike_n,template_depths(unique(spike_templates)), ...
        unique(spike_templates),20,'k','filled');

    str_1_unit_dots = scatter3(norm_spike_n(str_1_templates),template_depths(str_1_templates), ...
            str_1_templates,20,'green','filled');

    str_2_unit_dots = scatter3(norm_spike_n(str_2_templates),template_depths(str_2_templates), ...
            str_2_templates,20,'blue','filled');

    str_3_unit_dots = scatter3(norm_spike_n(str_3_templates),template_depths(str_3_templates), ...
            str_3_templates,20,'magenta','filled');

    yline(str_split, 'red')
    xlim(unit_axes,[-0.1,1]);
    ylim([-50, max(template_depths)+50]);
    ylabel('Depth (\mum)')
    xlabel('Normalized log rate')
    title(rec_day)
end

% legend([sharp_responsive_unit_dots wide_responsive_unit_dots both_responsive_unit_dots], {'Sharp resp', 'Wide resp', 'Both resp'})

str_split_fig_name = [animal '_str_cells_per_depth_Centre.tif'];
str_split_fig_path = fullfile(save_fig_path, str_split_fig_name);
saveas(str_split_fig, str_split_fig_path);


%% TEMP rt sort plot for a unit

rec_day = '2023-07-03';

good_template_depths = template_depths;
[sorted_good_template_depths, sorted_good_template_depth_idx] = sort(good_template_depths);

figure;
for idx=65:75
    plot(stim_norm_smooth_all_units(sorted_good_template_depth_idx(idx), :))
    hold on;
end
legend({num2str([65:75]')})

% 69 and 71 look good - unit 79 and 77
unit_idx = 79;

% get value of reaction time for bin
max_rt = ceil(prctile(stim_to_move(completed_rewarded_trials), 98));

% create time vector around all stim onsets
bin_window = 0.001;
bin_edges = -0.5:bin_window:max_rt;
around_stim_time = stimOn_times(completed_rewarded_trials) + bin_edges;

% for plot
bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;

% stim
unit_spikes = spike_templates == unit_idx;
unit_spike_times = spike_times_timelite(unit_spikes);
unit_spikes_stim = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, around_stim_time(trial_id,:))', ...
    1:size(around_stim_time, 1), 'UniformOutput',false))' / bin_window;

% smooth
gauss_win = gausswin(51, 3)';
stim_smooth_unit = filter(gauss_win,sum(gauss_win),unit_spikes_stim, [], 2);

% baseline and normalize
stim_baseline_smooth_unit = mean(stim_smooth_unit(:,bin_centres>-0.2 & bin_centres<0), 2);
stim_norm_smooth_unit = (stim_smooth_unit - stim_baseline_smooth_unit) ...
    ./ (stim_baseline_smooth_unit + std(stim_baseline_smooth_unit));

%     % check it looks the same
%     figure;
%     plot(mean(stim_norm_smooth_unit, 1))
%
%     figure;
%     idx = 71;
%     plot(stim_norm_smooth_all_units(sorted_good_template_depth_idx(idx), :))

% sort RT for plot
[sorted_rt, sort_rt_idx] = sort(stim_to_move);

% plot
unit_fig = figure('Position', get(0, 'Screensize')); %, Visible="off");
sgtitle(['Good unit ' num2str(unit_idx) ' day ' rec_day]);

upper_caxis = max(abs(stim_norm_smooth_unit), [],"all");
imagesc(bin_centres, [], stim_norm_smooth_unit(sort_rt_idx, :));
xline(0, 'LineWidth', 1);
hold on;
plot(sorted_rt, 1:length(sorted_rt), 'k', 'LineWidth', 1.5);
colormap(unit_fig, AP_colormap('BWR', [], 0.7));
if upper_caxis
    caxis([-upper_caxis upper_caxis]);
end
colorbar;


%% REGRESSION

save_fig_path = 'C:\Users\amarica\Documents\Lab stuff\Random figs\Long_str_stuff\Task';

animal = 'AM021';

load('AM021_swr_bhv_data.mat');
load([animal '_task_str_ephys_data']);
load([animal '_str_ephys_data.mat']);

% get learned day
bhv_days = swr_bhv_data.bhv_days{:};
ephys_days = passive_ephys_data.rec_day';

% remove empty ones
ephys_days = ephys_days(~cellfun('isempty',ephys_days));

%     strcmpi(bhv_days, ephys_days)

%     bhv_and_ephys = ismember(bhv_days, ephys_days);

days_from_learning = swr_bhv_data.days_from_learning{:};
days_from_learning = days_from_learning(ismember(bhv_days, ephys_days));


workflow = {'stim_wheel_right*'};
recordings = plab.find_recordings(animal, [], workflow);

% define struct for regression
regression = struct;
regression.days_from_learning = days_from_learning;


% define struct for regression per neuron
regression_per_neuron = struct;
regression_per_neuron.days_from_learning = days_from_learning;
for use_rec=3:length(recordings)

    rec_day = recordings(use_rec).day;
    rec_time = recordings(use_rec).recording{end};

    verbose = true;
    load_parts.widefield = false;
    load_parts.ephys = true;
    ap.load_recording

    % use_rec=length(recordings)

    this_rec = use_rec - 2;

    %% - get mua trace
    % get str boundaries
    str_depth = task_ephys_data.str_depth{this_rec};

    str_start = str_depth(1);
    str_end = str_depth(2);
    str_split_step = floor((str_end-str_start)/3);
    str_split = [str_start ...
        str_start+str_split_step ...
        str_start+2*str_split_step ...
        str_end];

    str_1_templates = task_ephys_data.str_1_templates{this_rec};
    str_2_templates = task_ephys_data.str_2_templates{this_rec};
    str_3_templates = task_ephys_data.str_3_templates{this_rec};
    cortex_templates = task_ephys_data.cortex_templates{this_rec};

    %% get time bins

    downsampled_time = downsample(timelite.timestamps, 10);

    end_timebin = downsampled_time(end)+mean(diff(downsampled_time));
    neural_downsampled_time = [downsampled_time; end_timebin];

    bin_window = mean(diff(downsampled_time));

    %% psths
    %% - stim

    good_templates_idx = find(good_templates);

    % get rewarded trials
    rewarded_trials = logical([trial_events.values.Outcome]);

    % get good reward times
    stimOff_times = photodiode_times(photodiode_values==0);
    [~, good_reward_idx] = min(abs(reward_times - stimOff_times(rewarded_trials)'));
    good_reward_times = reward_times(good_reward_idx);

    % get completed rew trials
    completed_trials = ~cellfun(@isempty,{trial_events.values.Outcome});
    completed_rewarded_trials = rewarded_trials;
    completed_rewarded_trials(find(completed_trials==0)) = false;

    % create time vector around all stim onsets
    bin_window = 0.01;
    bin_edges = -0.5:bin_window:2;
    around_stim_time = stimOn_times(completed_rewarded_trials) + bin_edges;

    % for plot
    bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;

    % task stim
    spikes_in_stim_time_downsampled = nan(length(good_templates_idx), length(find(completed_rewarded_trials)), size(around_stim_time, 2)-1);

    for unit_idx=1:length(good_templates_idx)

        unit_spikes = spike_templates == unit_idx;

        unit_spike_times = spike_times_timelite(unit_spikes);

        spikes_in_stim_time_downsampled(unit_idx, :, :) = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, around_stim_time(trial_id,:))', ...
            1:size(around_stim_time, 1), 'UniformOutput',false))' / bin_window;
    end

    % get avg across trials
    psth_stim_str_1_all = squeeze(mean(spikes_in_stim_time_downsampled(str_1_templates, :, :), 2));
    psth_stim_str_2_all = squeeze(mean(spikes_in_stim_time_downsampled(str_2_templates, :, :), 2));
    psth_stim_str_3_all = squeeze(mean(spikes_in_stim_time_downsampled(str_3_templates, :, :), 2));

    % baseline
    stim_smooth_baseline = task_ephys_data.stim_smooth_baseline{this_rec};
    str_1_stim_baseline = stim_smooth_baseline(str_1_templates);
    str_2_stim_baseline = stim_smooth_baseline(str_2_templates);
    str_3_stim_baseline = stim_smooth_baseline(str_3_templates);

    % baseline subtract
    psth_stim_str_1_all_baseline_sub = psth_stim_str_1_all - repmat(str_1_stim_baseline,1,size(psth_stim_str_1_all,2));
    psth_stim_str_1_mean_baseline_sub = mean(psth_stim_str_1_all_baseline_sub, 1)';

    psth_stim_str_2_all_baseline_sub = psth_stim_str_2_all - repmat(str_2_stim_baseline,1,size(psth_stim_str_2_all,2));
    psth_stim_str_2_mean_baseline_sub = mean(psth_stim_str_2_all_baseline_sub, 1)';

    psth_stim_str_3_all_baseline_sub = psth_stim_str_3_all - repmat(str_3_stim_baseline,1,size(psth_stim_str_3_all,2));
    psth_stim_str_3_mean_baseline_sub = mean(psth_stim_str_3_all_baseline_sub, 1)';

    % % plot to check
    % figure;
    % plot(psth_stim_str_1_mean_baseline_sub)


    %% long neural trace
    % LEFT HERE

    % multiunit striatum 1
    str_1_all_trace = cell2mat(arrayfun(@(unit_idx) histcounts(spike_times_timelite(spike_templates == unit_idx), neural_downsampled_time), ...
        str_1_templates, 'UniformOutput',false)) / bin_window;
    str_1_mua_trace = mean(str_1_all_trace, 1)';

    % multiunit striatum 2
    str_2_all_trace = cell2mat(arrayfun(@(unit_idx) histcounts(spike_times_timelite(spike_templates == unit_idx), neural_downsampled_time), ...
        str_2_templates, 'UniformOutput',false)) / bin_window;
    str_2_mua_trace = mean(str_2_all_trace, 1)';

    % multiunit striatum 3
    str_3_all_trace = cell2mat(arrayfun(@(unit_idx) histcounts(spike_times_timelite(spike_templates == unit_idx), neural_downsampled_time), ...
        str_3_templates, 'UniformOutput',false)) / bin_window;
    str_3_mua_trace = mean(str_3_all_trace, 1)';

    % multiunit cortex
    cortex_all_trace = cell2mat(arrayfun(@(unit_idx) histcounts(spike_times_timelite(spike_templates == unit_idx), neural_downsampled_time), ...
        cortex_templates, 'UniformOutput',false)) / bin_window;
    cortex_mua_trace = mean(cortex_all_trace, 1)';

    % subtract baseline
    stim_smooth_baseline = task_ephys_data.stim_smooth_baseline{this_rec};
    str_1_stim_baseline = stim_smooth_baseline(str_1_templates);
    str_1_all_trace_baseline_sub = str_1_all_trace - repmat(str_1_stim_baseline,1,size(str_1_all_trace,2));
    str_1_mua_trace_baseline_sub = mean(str_1_all_trace_baseline_sub, 1)';

    str_2_stim_baseline = stim_smooth_baseline(str_2_templates);
    str_2_all_trace_baseline_sub = str_2_all_trace - repmat(str_2_stim_baseline,1,size(str_2_all_trace,2));
    str_2_mua_trace_baseline_sub = mean(str_2_all_trace_baseline_sub, 1)';

    str_3_stim_baseline = stim_smooth_baseline(str_3_templates);
    str_3_all_trace_baseline_sub = str_3_all_trace - repmat(str_3_stim_baseline,1,size(str_3_all_trace,2));
    str_3_mua_trace_baseline_sub = mean(str_3_all_trace_baseline_sub, 1)';

    cortex_stim_baseline = stim_smooth_baseline(cortex_templates);
    cortex_all_trace_baseline_sub = cortex_all_trace - repmat(cortex_stim_baseline,1,size(cortex_all_trace,2));
    cortex_mua_trace_baseline_sub = mean(cortex_all_trace_baseline_sub, 1)';

    % figure;
    % imagesc(str_1_all_trace)
    % figure;
    % imagesc(str_1_all_trace_baseline_sub)
    %
    % figure;
    % plot(str_1_mua_trace)
    % hold on
    % plot(str_1_mua_trace_baseline_sub)
    %
    % figure;
    % plot(str_2_mua_trace)
    % hold on
    % plot(str_2_mua_trace_baseline_sub)
    %
    % figure;
    % plot(str_3_mua_trace)
    % hold on
    % plot(str_3_mua_trace_baseline_sub)

    % % % test_str_1_baseline_mean = mean(str_1_stim_baseline);
    % % % test_str_1_mua_trace_baseline_sub = str_1_mua_trace - test_str_1_baseline_mean;
    % % %
    % % % figure;
    % % % plot(str_1_mua_trace)
    % % % hold on
    % % % plot(test_str_1_mua_trace_baseline_sub)


    %% - regressors

    % get move times
    move_on = find([0; diff(wheel_move)] == 1);
    moveOn_times = timelite.timestamps(move_on);

    stim_move_times = stimOn_times(completed_trials) + stim_to_move;
    iti_move_times = setdiff(moveOn_times, stim_move_times);

    % define time shifts 
    stim_t_interval_for_shift = [0, 0.5];
    stim_move_t_interval_for_shift = [-0.5 1];
    iti_move_t_interval_for_shift = [-0.5, 1];
    reward_t_interval_for_shift = [0, 0.5];

%     %% - do cross-val regression from baseline sub mua trace
% 
%     str_1_var_expl = struct;
%     str_2_var_expl = struct;
%     str_3_var_expl = struct;
% 
% %     % ----- move only ---------------------------------------------
% %     % str 1
% %     [str_1_move_avg_coeff, str_1_move_predict_mua_trace, str_1_var_expl.move] = ...
% %         ephys.run_regression(str_1_mua_trace_baseline_sub, ...
% %         {stim_move_times, iti_move_times}, ...
% %         {stim_move_t_interval_for_shift, iti_move_t_interval_for_shift}, ...
% %         downsampled_time);
% %     % str 2
% %     [str_2_move_avg_coeff, str_2_move_predict_mua_trace, str_2_var_expl.move] = ...
% %         ephys.run_regression(str_2_mua_trace_baseline_sub, ...
% %         {stim_move_times, iti_move_times}, ...
% %         {stim_move_t_interval_for_shift, iti_move_t_interval_for_shift}, ...
% %         downsampled_time);
% %     % str 3
% %     [str_3_move_avg_coeff, str_3_move_predict_mua_trace, str_3_var_expl.move] = ...
% %         ephys.run_regression(str_3_mua_trace_baseline_sub, ...
% %         {stim_move_times, iti_move_times}, ...
% %         {stim_move_t_interval_for_shift, iti_move_t_interval_for_shift}, ...
% %         downsampled_time);
% % 
% %     % ------- reward only -------------------------------------------
% %     % str 1
% %     [str_1_rew_avg_coeff, str_1_rew_predict_mua_trace, str_1_var_expl.rew] = ephys.run_regression(str_1_mua_trace_baseline_sub, ...
% %         {reward_times}, {reward_t_interval_for_shift}, ...
% %         downsampled_time);
% %     % str 2
% %     [str_2_rew_avg_coeff, str_2_rew_predict_mua_trace, str_2_var_expl.rew] = ephys.run_regression(str_2_mua_trace_baseline_sub, ...
% %         {reward_times}, {reward_t_interval_for_shift}, ...
% %         downsampled_time);
% %     % str 3
% %     [str_3_rew_avg_coeff, str_3_rew_predict_mua_trace, str_3_var_expl.rew] = ephys.run_regression(str_3_mua_trace_baseline_sub, ...
% %         {reward_times}, {reward_t_interval_for_shift}, ...
% %         downsampled_time);
% % 
% %     % ------- stim only ----------------------------------------------
% %     % str 1
% %     [str_1_stim_avg_coeff, str_1_stim_predict_mua_trace, str_1_var_expl.stim] = ...
% %         ephys.run_regression(str_1_mua_trace_baseline_sub, ...
% %         {stimOn_times}, {stim_t_interval_for_shift}, ...
% %         downsampled_time);
% %     % str 2
% %     [str_2_stim_avg_coeff, str_2_stim_predict_mua_trace, str_2_var_expl.stim] = ...
% %         ephys.run_regression(str_2_mua_trace_baseline_sub, ...
% %         {stimOn_times}, {stim_t_interval_for_shift}, ...
% %         downsampled_time);
% %     % str 3
% %     [str_3_stim_avg_coeff, str_3_stim_predict_mua_trace, str_3_var_expl.stim] = ...
% %         ephys.run_regression(str_3_mua_trace_baseline_sub, ...
% %         {stimOn_times}, {stim_t_interval_for_shift}, ...
% %         downsampled_time);
% 
% %     % ---------- stim and move -----------------------------------------
% %     % str 1
% %     [str_1_stim_move_avg_coeff, str_1_stim_move_predict_mua_trace, str_1_var_expl.stim_move] = ... 
% %         ephys.run_regression(str_1_mua_trace_baseline_sub, ...
% %         {stimOn_times, stim_move_times, iti_move_times}, ...
% %         {stim_t_interval_for_shift, stim_move_t_interval_for_shift, iti_move_t_interval_for_shift}, ...
% %         downsampled_time);
% %     % str 2
% %     [str_2_stim_move_avg_coeff, str_2_stim_move_predict_mua_trace, str_2_var_expl.stim_move] = ... 
% %         ephys.run_regression(str_2_mua_trace_baseline_sub, ...
% %         {stimOn_times, stim_move_times, iti_move_times}, ...
% %         {stim_t_interval_for_shift, stim_move_t_interval_for_shift, iti_move_t_interval_for_shift}, ...
% %         downsampled_time);
% %     % str 3
% %     [str_3_stim_move_avg_coeff, str_3_stim_move_predict_mua_trace, str_3_var_expl.stim_move] = ... 
% %         ephys.run_regression(str_3_mua_trace_baseline_sub, ...
% %         {stimOn_times, stim_move_times, iti_move_times}, ...
% %         {stim_t_interval_for_shift, stim_move_t_interval_for_shift, iti_move_t_interval_for_shift}, ...
% %         downsampled_time);
% 
%     % -------------- move and reward ----------------------------------
%     % str 1
%     [str_1_move_rew_avg_coeff, str_1_move_rew_predict_mua_trace, str_1_var_expl.move_rew] = ephys.run_regression(str_1_mua_trace_baseline_sub, ...
%         {stim_move_times, iti_move_times, reward_times}, ...
%         {stim_move_t_interval_for_shift, iti_move_t_interval_for_shift, ...
%         reward_t_interval_for_shift}, ...
%         downsampled_time);
%     % str 2
%     [str_2_move_rew_avg_coeff, str_2_move_rew_predict_mua_trace, str_2_var_expl.move_rew] = ephys.run_regression(str_2_mua_trace_baseline_sub, ...
%         {stim_move_times, iti_move_times, reward_times}, ...
%         {stim_move_t_interval_for_shift, iti_move_t_interval_for_shift, ...
%         reward_t_interval_for_shift}, ...
%         downsampled_time);
%     % str 3
%     [str_3_move_rew_avg_coeff, str_3_move_rew_predict_mua_trace, str_3_var_expl.move_rew] = ephys.run_regression(str_3_mua_trace_baseline_sub, ...
%         {stim_move_times, iti_move_times, reward_times}, ...
%         {stim_move_t_interval_for_shift, iti_move_t_interval_for_shift, ...
%         reward_t_interval_for_shift}, ...
%         downsampled_time);
% 
%     % ---------- stim move rew  ---------------------------------------------
%     % str 1
%     [str_1_stim_move_rew_avg_coeff, str_1_stim_move_rew_predict_mua_trace, str_1_var_expl.stim_move_rew] = ephys.run_regression(str_1_mua_trace_baseline_sub, ...
%         {stimOn_times, stim_move_times, iti_move_times, reward_times}, ...
%         {stim_t_interval_for_shift, stim_move_t_interval_for_shift, ...
%         iti_move_t_interval_for_shift, reward_t_interval_for_shift}, ...
%         downsampled_time);
%     % str 2
%     [str_2_stim_move_rew_avg_coeff, str_2_stim_move_rew_predict_mua_trace, str_2_var_expl.stim_move_rew] = ephys.run_regression(str_2_mua_trace_baseline_sub, ...
%         {stimOn_times, stim_move_times, iti_move_times, reward_times}, ...
%         {stim_t_interval_for_shift, stim_move_t_interval_for_shift, ...
%         iti_move_t_interval_for_shift, reward_t_interval_for_shift}, ...
%         downsampled_time);
%     % str 3
%     [str_3_stim_move_rew_avg_coeff, str_3_stim_move_rew_predict_mua_trace, str_3_var_expl.stim_move_rew] = ephys.run_regression(str_3_mua_trace_baseline_sub, ...
%         {stimOn_times, stim_move_times, iti_move_times, reward_times}, ...
%         {stim_t_interval_for_shift, stim_move_t_interval_for_shift, ...
%         iti_move_t_interval_for_shift, reward_t_interval_for_shift}, ...
%         downsampled_time);
% 
% 
%     disp('Str 1:')
%     disp(str_1_var_expl)
%     disp('Str 2:')
%     disp(str_2_var_expl)
%     disp('Str 3:')
%     disp(str_3_var_expl)
% 
% %     figure;
% %     plot(str_1_mua_trace_baseline_sub);
% %     hold on;
% %     plot(str_1_stim_move_rew_predict_mua_trace)
% %     title([rec_day ' '])
% % 
% % 
% %     figure = tiledlayout(3, 1);
% %     nexttile
% %     plot(str_1_stim_move_rew_avg_coeff(1:50))
% %     hold on;
% %     plot(str_1_stim_move_rew_avg_coeff(51:200))
% %     hold on;
% %     plot(str_1_stim_move_rew_avg_coeff(201:350))
% %     hold on;
% %     plot(str_1_stim_move_rew_avg_coeff(351:400))
% %     legend({'Stim', 'Stim Move', 'ITI Move', 'Reward'})
% %     title([rec_day ' Kernels for mua trace str 1 FUNCTION'])
% % 
% %     nexttile
% %     plot(str_2_stim_move_rew_avg_coeff(1:50))
% %     hold on;
% %     plot(str_2_stim_move_rew_avg_coeff(51:200))
% %     hold on;
% %     plot(str_2_stim_move_rew_avg_coeff(201:350))
% %     hold on;
% %     plot(str_2_stim_move_rew_avg_coeff(351:400))
% %     legend({'Stim', 'Stim Move', 'ITI Move', 'Reward'})
% %     title([rec_day ' Kernels for mua trace str 2 FUNCTION'])
% % 
% %     nexttile
% %     plot(str_3_stim_move_rew_avg_coeff(1:50))
% %     hold on;
% %     plot(str_3_stim_move_rew_avg_coeff(51:200))
% %     hold on;
% %     plot(str_3_stim_move_rew_avg_coeff(201:350))
% %     hold on;
% %     plot(str_3_stim_move_rew_avg_coeff(351:400))
% %     legend({'Stim', 'Stim Move', 'ITI Move', 'Reward'})
% %     title([rec_day ' Kernels for mua trace str 3 FUNCTION'])
% 
%     % do reg from residuals
% 
%     str_1_mua_trace_residual = str_1_mua_trace_baseline_sub - str_1_move_rew_predict_mua_trace;
%     [str_1_residual_stim_avg_coeff, str_1_residual_stim_predict_mua_trace, ...
%         str_1_var_expl.residual_stim] = ephys.run_regression(str_1_mua_trace_residual, ...
%         {stimOn_times}, {stim_t_interval_for_shift}, ...
%         downsampled_time);
%     str_2_mua_trace_residual = str_2_mua_trace_baseline_sub - str_2_move_rew_predict_mua_trace;
%     [str_2_residual_stim_avg_coeff, str_2_residual_stim_predict_mua_trace, ...
%         str_2_var_expl.residual_stim] = ephys.run_regression(str_2_mua_trace_residual, ...
%         {stimOn_times}, {stim_t_interval_for_shift}, ...
%         downsampled_time);
%     str_3_mua_trace_residual = str_3_mua_trace_baseline_sub - str_3_move_rew_predict_mua_trace;
%     [str_3_residual_stim_avg_coeff, str_3_residual_stim_predict_mua_trace, ...
%         str_3_var_expl.residual_stim] = ephys.run_regression(str_3_mua_trace_residual, ...
%         {stimOn_times}, {stim_t_interval_for_shift}, ...
%         downsampled_time);
% 
%     figure
%     plot(str_1_residual_stim_avg_coeff)
% 
% 
%     % plot residuals aligned to stim
%     test1 = interp1(downsampled_time, str_1_mua_trace_residual, around_stim_time);
%     mean_test1 = mean(test1, 1);
% 
%     test2 = interp1(downsampled_time, str_2_mua_trace_residual, around_stim_time);
%     mean_test2 = mean(test2, 1);
% 
%     test3 = interp1(downsampled_time, str_3_mua_trace_residual, around_stim_time);
%     mean_test3 = mean(test3, 1);
% 
% 
%     % plot predicted aligned to stim
%     pred_test1 = interp1(downsampled_time, str_1_residual_stim_predict_mua_trace, around_stim_time);
%     mean_pred_test1 = mean(pred_test1, 1);
% 
%     pred_test2 = interp1(downsampled_time, str_2_residual_stim_predict_mua_trace, around_stim_time);
%     mean_pred_test2 = mean(pred_test2, 1);
% 
%     pred_test3 = interp1(downsampled_time, str_3_residual_stim_predict_mua_trace, around_stim_time);
%     mean_pred_test3 = mean(pred_test3, 1);
% 
% 
%     figure = tiledlayout(3,1)
%     nexttile
%     plot(mean_test1);
%     hold on;
%     plot(mean_pred_test1)
%     title('Str 1')
% 
%     nexttile
%     plot(mean_test2);
%     hold on;
%     plot(mean_pred_test2)
%     title('Str 2')
% 
%     nexttile
%     plot(mean_test3);
%     hold on;
%     plot(mean_pred_test3)
%     title('Str 3')


    %% do regression on each neuron
% TEMP COM OUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     str_all_templates = vertcat(str_1_templates, str_2_templates, str_3_templates);
%     str_all_all_trace_baseline_sub = vertcat(str_1_all_trace_baseline_sub, ...
%         str_2_all_trace_baseline_sub, ...
%         str_3_all_trace_baseline_sub);
%     mean_neuron_stim_align_str_all_trace_baseline_sub = cell(1, length(str_all_templates));
%     mean_neuron_str_all_trace_residual_stim_align = cell(1, length(str_all_templates));
% 
%     % move rew
%     move_rew_design_matrix = ephys.make_design_matrix(...
%         {stim_move_times, iti_move_times, reward_times}, ...
%         {stim_move_t_interval_for_shift, iti_move_t_interval_for_shift, ...
%         reward_t_interval_for_shift}, ...
%         downsampled_time);
%     % TAKE EACH NEURON
%     neuron_str_all_move_rew_avg_coeff = cell(1, length(str_all_templates));
%     neuron_str_all_move_rew_predict_trace = cell(1, length(str_all_templates));
%     neuron_str_all_trace_residual = cell(1, length(str_all_templates));
%     neuron_str_all_trace_residual_stim_align = cell(1, length(str_all_templates));
% 
%     for neuron_idx = 1:length(str_all_templates)
%         neuron_str_all_trace_baseline_sub = str_all_all_trace_baseline_sub(neuron_idx, :)';
%         neuron_stim_align_str_all_trace_baseline_sub = interp1(downsampled_time, neuron_str_all_trace_baseline_sub, around_stim_time);
%         mean_neuron_stim_align_str_all_trace_baseline_sub{neuron_idx} = mean(neuron_stim_align_str_all_trace_baseline_sub, 1);
%         %     figure; imagesc(neuron_stim_align_str_all_trace_baseline_sub)
%         %     figure; plot(mean_neuron_stim_align_str_all_trace_baseline_sub)
% 
%         % move rew
%         [neuron_str_all_move_rew_avg_coeff{neuron_idx}, neuron_str_all_move_rew_predict_trace{neuron_idx}, neuron_str_all_var_expl(neuron_idx).move_rew] = ...
%             ephys.run_regression(neuron_str_all_trace_baseline_sub, ...
%             move_rew_design_matrix, ...
%             downsampled_time);
% 
% %         [k,predicted_signals,explained_var,predicted_signals_reduced] = ...
% %             AP_regresskernel(regressors,neuron_str_all_trace_baseline_sub,t_shifts)
% 
%         % residuals
%         neuron_str_all_trace_residual{neuron_idx} = neuron_str_all_trace_baseline_sub - neuron_str_all_move_rew_predict_trace{neuron_idx};
% 
%         % signed rank test
%         neuron_str_all_trace_residual_stim_align{neuron_idx} = interp1(downsampled_time, neuron_str_all_trace_residual{neuron_idx}, around_stim_time);
%         mean_neuron_str_all_trace_residual_stim_align{neuron_idx} = mean(neuron_str_all_trace_residual_stim_align{neuron_idx}, 1);
% 
% %         pre_stim = mean(neuron_str_all_trace_residual_stim_align(:, bin_edges>=-bin_window_for_test&bin_edges<0), 2);
% %         post_stim = mean(neuron_str_all_trace_residual_stim_align(:, bin_edges<=bin_window_for_test&bin_edges>0), 2);
% % 
% %         str_all_sign_rank_test_p_val(neuron_idx) = signrank(post_stim, pre_stim, 'tail', 'right');
% %             figure;
% %             plot(pre_stim)
% %             hold on
% %             plot(post_stim)
% % %         
% % %             figure
% % %             plot(post_stim-pre_stim)
% % % % 
% % % %             figure;
% % % %             plot(abs(post_stim)-abs(pre_stim))
% % % % 
% % %             figure;
% % %             plot(bin_edges, mean_neuron_stim_align_str_all_trace_baseline_sub{neuron_idx})
% % %            xline([-bin_window_for_test bin_window_for_test])
% % %         %     xline(-)
% % %         
% 
%             upper_caxis = max(max(neuron_str_all_trace_residual_stim_align{neuron_idx}, [], 'all'), max(neuron_stim_align_str_all_trace_baseline_sub, [], 'all'));
%             figure; 
%             tiledlayout(2,2);
%             nexttile;
%             imagesc(bin_centres, [], neuron_str_all_trace_residual_stim_align{neuron_idx})
%             colormap(AP_colormap('BWR', [], 0.7));
%             clim([-upper_caxis upper_caxis]);
%             ylabel('Residual', 'FontWeight','bold')
%             nexttile;
%             plot(bin_edges, mean_neuron_str_all_trace_residual_stim_align{neuron_idx})
%             xline(0)
%             
% 
%             raw_subfig = nexttile;
%             imagesc(bin_centres, [], neuron_stim_align_str_all_trace_baseline_sub)
%             ylabel('Raw', 'FontWeight','bold')
%             colormap(AP_colormap('BWR', [], 0.7));
%             clim([-upper_caxis upper_caxis]);
%             nexttile;
%             plot(bin_edges, mean_neuron_stim_align_str_all_trace_baseline_sub{neuron_idx})
%             xline(0)
%     end
% 
% 
% 
%     % try different windows
% 
%     bin_window_for_test = 0.3;
%     str_all_sign_rank_test_p_val = nan(1, length(str_all_templates));
%     pre_stim = cell(1, length(str_all_templates));
%     post_stim = cell(1, length(str_all_templates));
%     for neuron_idx = 1:length(str_all_templates)
%         neuron_str_all_trace_residual_stim_align = interp1(downsampled_time, neuron_str_all_trace_residual{neuron_idx}, around_stim_time);
%         pre_stim{neuron_idx} = mean(neuron_str_all_trace_residual_stim_align(:, bin_edges>=-bin_window_for_test&bin_edges<0), 2);
%         post_stim{neuron_idx} = mean(neuron_str_all_trace_residual_stim_align(:, bin_edges<=bin_window_for_test&bin_edges>0), 2);
%         max_post_stim{neuron_idx} = max(neuron_str_all_trace_residual_stim_align(:, bin_edges<=bin_window_for_test&bin_edges>0), 2);
%         str_all_sign_rank_test_p_val(neuron_idx) = signrank(post_stim{neuron_idx}, pre_stim{neuron_idx}, 'tail', 'right');
%         str_all_max_sign_rank_test_p_val(neuron_idx) = signrank(max_post_stim{neuron_idx}, pre_stim{neuron_idx}, 'tail', 'right');
%     end
% figure;
%     tiledlayout('flow');
%     for neuron_idx = 1:length(str_all_templates)
%         if str_all_sign_rank_test_p_val(neuron_idx)<0.05
%             nexttile
%             plot(bin_edges, mean_neuron_str_all_trace_residual_stim_align{neuron_idx})
%             xline([-bin_window_for_test bin_window_for_test])
%             title(['Neuron ' num2str(neuron_idx) ' test val ' num2str(str_all_sign_rank_test_p_val(neuron_idx))])
%         end         
%     end
%     sgtitle(['Normal - Residual trace Responsive bin ' num2str(bin_window_for_test)])
%     
% 
%     figure;
%     tiledlayout('flow');
%     for neuron_idx = 1:length(str_all_templates)
%         if str_all_sign_rank_test_p_val(neuron_idx)>=0.05
%             nexttile
%             plot(bin_edges, mean_neuron_str_all_trace_residual_stim_align{neuron_idx})
%             xline([-bin_window_for_test bin_window_for_test])
%             title(['Neuron ' num2str(neuron_idx) ' test val ' num2str(str_all_sign_rank_test_p_val(neuron_idx))])
%         end
%     end
%     sgtitle(['Normal - Residual trace Non-Responsive bin ' num2str(bin_window_for_test)])
% 
% 
% 
%     %     figure;
% %     tiledlayout('flow');
% %     for neuron_idx = 1:10%length(str_1_templates)
% %         if str_1_sign_rank_test_p_val(neuron_idx)>=0.05
% %             nexttile
% %             plot(bin_edges, mean_neuron_stim_align_str_1_trace_baseline_sub{neuron_idx})
% %             xline([-bin_window_for_test bin_window_for_test])
% %             title(['Neuron ' num2str(neuron_idx) ' test val ' num2str(str_1_sign_rank_test_p_val(neuron_idx))])
% %         end
% %     end
% %     sgtitle('Full trace Non-Responsive')
% 
% %     figure;
% %     tiledlayout('flow');
% %     for neuron_idx = 1:10%length(str_1_templates)
% %         if str_1_sign_rank_test_p_val(neuron_idx)<0.05
% %             nexttile
% %             plot(bin_edges, mean_neuron_stim_align_str_1_trace_baseline_sub{neuron_idx})
% %             xline([-bin_window_for_test bin_window_for_test])
% %             title(['Neuron ' num2str(neuron_idx) ' test val ' num2str(str_1_sign_rank_test_p_val(neuron_idx))])
% %         end         
% %     end
% %     sgtitle('Full trace Responsive')
% TEMP COM OUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% TEMP %%%%%%%%%%%%%

    mean_neuron_stim_align_str_1_trace_baseline_sub = cell(1, length(str_1_templates));
    mean_neuron_str_1_trace_residual_stim_align = cell(1, length(str_1_templates));

    % move rew
    move_rew_design_matrix = ephys.make_design_matrix(...
        {stim_move_times, iti_move_times, reward_times}, ...
        {stim_move_t_interval_for_shift, iti_move_t_interval_for_shift, ...
        reward_t_interval_for_shift}, ...
        downsampled_time);
 % TAKE EACH NEURON
    neuron_str_1_move_rew_avg_coeff = cell(1, length(str_1_templates));
    neuron_str_1_move_rew_predict_trace = cell(1, length(str_1_templates));
    neuron_str_1_trace_residual = cell(1, length(str_1_templates));
    neuron_str_1_trace_residual_stim_align = cell(1, length(str_1_templates));

    for neuron_idx = 1:length(str_1_templates)
        tic
        neuron_str_1_trace_baseline_sub = str_1_all_trace_baseline_sub(neuron_idx, :)';
        neuron_stim_align_str_1_trace_baseline_sub = interp1(downsampled_time, neuron_str_1_trace_baseline_sub, around_stim_time);
        mean_neuron_stim_align_str_1_trace_baseline_sub{neuron_idx} = mean(neuron_stim_align_str_1_trace_baseline_sub, 1);
        %     figure; imagesc(neuron_stim_align_str_1_trace_baseline_sub)
        %     figure; plot(mean_neuron_stim_align_str_1_trace_baseline_sub)

        % move rew
        [neuron_str_1_move_rew_avg_coeff{neuron_idx}, neuron_str_1_move_rew_predict_trace{neuron_idx}, neuron_str_1_var_expl(neuron_idx).move_rew] = ...
            ephys.run_regression(neuron_str_1_trace_baseline_sub, ...
            move_rew_design_matrix, ...
            downsampled_time);

%         [k,predicted_signals,explained_var,predicted_signals_reduced] = ...
%             AP_regresskernel(regressors,neuron_str_1_trace_baseline_sub,t_shifts)

        % residuals
        neuron_str_1_trace_residual{neuron_idx} = neuron_str_1_trace_baseline_sub - neuron_str_1_move_rew_predict_trace{neuron_idx};

        % signed rank test
        neuron_str_1_trace_residual_stim_align{neuron_idx} = interp1(downsampled_time, neuron_str_1_trace_residual{neuron_idx}, around_stim_time);
        mean_neuron_str_1_trace_residual_stim_align{neuron_idx} = mean(neuron_str_1_trace_residual_stim_align{neuron_idx}, 1);

%         pre_stim = mean(neuron_str_1_trace_residual_stim_align(:, bin_edges>=-bin_window_for_test&bin_edges<0), 2);
%         post_stim = mean(neuron_str_1_trace_residual_stim_align(:, bin_edges<=bin_window_for_test&bin_edges>0), 2);
% 
%         str_1_sign_rank_test_p_val(neuron_idx) = signrank(post_stim, pre_stim, 'tail', 'right');
%             figure;
%             plot(pre_stim)
%             hold on
%             plot(post_stim)
% %         
% %             figure
% %             plot(post_stim-pre_stim)
% % % 
% % %             figure;
% % %             plot(abs(post_stim)-abs(pre_stim))
% % % 
% %             figure;
% %             plot(bin_edges, mean_neuron_stim_align_str_1_trace_baseline_sub{neuron_idx})
% %            xline([-bin_window_for_test bin_window_for_test])
% %         %     xline(-)
% %         

%             upper_caxis = max(max(neuron_str_1_trace_residual_stim_align{neuron_idx}, [], 'all'), max(neuron_stim_align_str_1_trace_baseline_sub, [], 'all'));
%             figure; 
%             tiledlayout(2,2);
%             nexttile;
%             imagesc(bin_centres, [], neuron_str_1_trace_residual_stim_align{neuron_idx})
%             colormap(AP_colormap('BWR', [], 0.7));
%             clim([-upper_caxis upper_caxis]);
%             ylabel('Residual', 'FontWeight','bold')
%             nexttile;
%             plot(bin_edges, mean_neuron_str_1_trace_residual_stim_align{neuron_idx})
%             xline(0)
%             
% 
%             raw_subfig = nexttile;
%             imagesc(bin_centres, [], neuron_stim_align_str_1_trace_baseline_sub)
%             ylabel('Raw', 'FontWeight','bold')
%             colormap(AP_colormap('BWR', [], 0.7));
%             clim([-upper_caxis upper_caxis]);
%             nexttile;
%             plot(bin_edges, mean_neuron_stim_align_str_1_trace_baseline_sub{neuron_idx})
%             xline(0)

                time_loop = toc

    end


    % try different windows

    bin_window_for_test = 0.3;
    str_1_sign_rank_test_p_val = nan(1, length(str_1_templates));
    pre_stim = cell(1, length(str_1_templates));
    post_stim = cell(1, length(str_1_templates));
    for neuron_idx = 1:length(str_1_templates)
        neuron_str_1_trace_residual_stim_align = interp1(downsampled_time, neuron_str_1_trace_residual{neuron_idx}, around_stim_time);

        % zero everything below 0
        neuron_str_1_trace_residual_stim_align(neuron_str_1_trace_residual_stim_align<0) = 0;
        pre_stim{neuron_idx} = mean(neuron_str_1_trace_residual_stim_align(:, bin_edges>=-bin_window_for_test&bin_edges<0), 2);
        post_stim{neuron_idx} = mean(neuron_str_1_trace_residual_stim_align(:, bin_edges<=bin_window_for_test&bin_edges>0), 2);
        max_post_stim{neuron_idx} = max(neuron_str_1_trace_residual_stim_align(:, bin_edges<=bin_window_for_test&bin_edges>0), [], 2);
        max_pre_stim{neuron_idx} = max(neuron_str_1_trace_residual_stim_align(:, bin_edges>=-bin_window_for_test&bin_edges<0),[], 2);
        str_1_sign_rank_test_p_val(neuron_idx) = signrank(post_stim{neuron_idx}, pre_stim{neuron_idx}, 'tail', 'right');
        str_1_max_sign_rank_test_p_val(neuron_idx) = signrank(max_post_stim{neuron_idx}, max_pre_stim{neuron_idx}, 'tail', 'right');
    end


    figure;
    tiledlayout('flow');
    for neuron_idx = 1:length(str_1_templates)
        if str_1_sign_rank_test_p_val(neuron_idx)<0.05
            nexttile
            plot(bin_edges, mean_neuron_str_1_trace_residual_stim_align{neuron_idx})
            xline([-bin_window_for_test bin_window_for_test])
            title(['Neuron ' num2str(neuron_idx) ' test val ' num2str(str_1_sign_rank_test_p_val(neuron_idx))])
        end         
    end
    sgtitle(['Normal - Residual trace Responsive bin ' num2str(bin_window_for_test)])
    

    figure;
    tiledlayout('flow');
    for neuron_idx = 1:length(str_1_templates)
        if str_1_sign_rank_test_p_val(neuron_idx)>=0.05
            nexttile
            plot(bin_edges, mean_neuron_str_1_trace_residual_stim_align{neuron_idx})
            xline([-bin_window_for_test bin_window_for_test])
            title(['Neuron ' num2str(neuron_idx) ' test val ' num2str(str_1_sign_rank_test_p_val(neuron_idx))])
        end
    end
    sgtitle(['Normal - Residual trace Non-Responsive bin ' num2str(bin_window_for_test)])
    

    figure;
    tiledlayout('flow');
    for neuron_idx = 1:length(str_1_templates)
        if str_1_max_sign_rank_test_p_val(neuron_idx)<0.05
            nexttile
            plot(bin_edges, mean_neuron_str_1_trace_residual_stim_align{neuron_idx})
            xline([-bin_window_for_test bin_window_for_test])
            title(['Neuron ' num2str(neuron_idx) ' test val ' num2str(str_1_max_sign_rank_test_p_val(neuron_idx))])
        end         
    end
    sgtitle(['Max - Residual trace Responsive bin ' num2str(bin_window_for_test)])
    

    figure;
    tiledlayout('flow');
    for neuron_idx = 1:length(str_1_templates)
        if str_1_max_sign_rank_test_p_val(neuron_idx)>=0.05
            nexttile
            plot(bin_edges, mean_neuron_str_1_trace_residual_stim_align{neuron_idx})
            xline([-bin_window_for_test bin_window_for_test])
            title(['Neuron ' num2str(neuron_idx) ' test val ' num2str(str_1_max_sign_rank_test_p_val(neuron_idx))])
        end
    end
    sgtitle(['Max - Residual trace Non-Responsive bin ' num2str(bin_window_for_test)])


    neuron_idx = 79;
    figure;
    plot(max_post_stim{neuron_idx}, 'o')
    hold on;
    plot(max_pre_stim{neuron_idx}, 'o')

    figure;
    plot(max_post_stim{neuron_idx}-max_pre_stim{neuron_idx})
    median(max_post_stim{neuron_idx}-max_pre_stim{neuron_idx})

    neuron_idx = 75;
    figure;
    plot(max_post_stim{neuron_idx}, 'o')
    hold on;
    plot(max_pre_stim{neuron_idx}, 'o')
    title('75')


    figure;
    plot(max_post_stim{neuron_idx}-max_pre_stim{neuron_idx})
    median(max_post_stim{neuron_idx}-max_pre_stim{neuron_idx})
    title('75')

    figure;
    imagesc(neuron_stim_align_str_1_trace_baseline_sub)
    title('75')

    neuron_str_1_trace = str_1_all_trace(neuron_idx, :)';
    neuron_stim_align_str_1_trace = interp1(downsampled_time, neuron_str_1_trace, around_stim_time);
       
    figure;
    imagesc(neuron_stim_align_str_1_trace)
    title('75 raw')
%%%% END TEMP %%%%%%%%%%%%%%
    

    %% compare to psths
    %% - baseline sub kernels and baseline sub psth
%     figure;
%     plot(str_1_stim_move_rew_avg_coeff(1:50))
%     hold on;
%     plot(psth_stim_str_1_mean_baseline_sub(51:100))
%     legend({'Kernel', 'PSTH'})
%     title([rec_day ' Str 1 PSTH and reg kernel comparison for stim'])

    %% - aligned predict trace and baseline sub psth
%     % align predict trace
%     predict_mua_trace_stim_align_stim_move_rew = interp1(downsampled_time, str_1_stim_move_rew_predict_mua_trace, around_stim_time);
%     predict_mua_trace_stim_align_mean = mean(predict_mua_trace_stim_align_stim_move_rew, 1);
% 
%     figure;
%     plot(predict_mua_trace_stim_align_mean)
%     hold on;
%     plot(psth_stim_str_1_mean_baseline_sub)
%     legend({'Predict trace aligned', 'PSTH'})
%     title([rec_day ' Str 1 PSTH and predict trace comparison for stim'])

    %% - kernel and passive psth
%     passive_rec_day = find(ismember(passive_ephys_data.rec_day, rec_day));
%     passive_stimOn_times = passive_ephys_data.stimOn_times{passive_rec_day};
%     passive_spike_times_timelite = passive_ephys_data.spike_times_timelite{passive_rec_day};
%     passive_spike_templates = passive_ephys_data.spike_templates{passive_rec_day};
%     passive_good_trials = passive_ephys_data.contra_good_trials{passive_rec_day};
% 
%     % create time vector around all stim onsets
%     bin_window = 0.01;
%     bin_edges = -0.5:bin_window:2;
%     passive_around_stim_time = passive_stimOn_times(passive_good_trials) + bin_edges;
% 
%     % for plot
%     bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;
% 
%     % task stim
%     passive_spikes_in_stim_time_downsampled = nan(length(good_templates_idx), length(passive_stimOn_times(passive_good_trials)), size(passive_around_stim_time, 2)-1);
% 
%     for unit_idx=1:length(good_templates_idx)
% 
%         unit_spikes = passive_spike_templates == unit_idx;
% 
%         unit_spike_times = passive_spike_times_timelite(unit_spikes);
% 
%         passive_spikes_in_stim_time_downsampled(unit_idx, :, :) = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, passive_around_stim_time(trial_id,:))', ...
%             1:size(passive_around_stim_time, 1), 'UniformOutput',false))' / bin_window;
%     end
% 
%     % get avg across trials
%     passive_psth_stim_str_1_all = squeeze(mean(passive_spikes_in_stim_time_downsampled(str_1_templates, :, :), 2));
%     passive_psth_stim_str_2_all = squeeze(mean(passive_spikes_in_stim_time_downsampled(str_2_templates, :, :), 2));
%     passive_psth_stim_str_3_all = squeeze(mean(passive_spikes_in_stim_time_downsampled(str_3_templates, :, :), 2));
% 
%     % baseline
%     passive_str_1_stim_baseline = mean(passive_psth_stim_str_1_all(bin_edges<0), 2);
%     passive_str_2_stim_baseline = mean(passive_psth_stim_str_2_all(bin_edges<0), 2);
%     passive_str_3_stim_baseline = mean(passive_psth_stim_str_3_all(bin_edges<0), 2);
% 
%     % baseline subtract
%     passive_psth_stim_str_1_all_baseline_sub = passive_psth_stim_str_1_all - repmat(passive_str_1_stim_baseline,1,size(passive_psth_stim_str_1_all,2));
%     passive_psth_stim_str_1_mean_baseline_sub = mean(passive_psth_stim_str_1_all_baseline_sub, 1)';
% 
%     passive_psth_stim_str_2_all_baseline_sub = passive_psth_stim_str_2_all - repmat(passive_str_2_stim_baseline,1,size(passive_psth_stim_str_2_all,2));
%     passive_psth_stim_str_2_mean_baseline_sub = mean(passive_psth_stim_str_2_all_baseline_sub, 1)';
% 
%     passive_psth_stim_str_3_all_baseline_sub = passive_psth_stim_str_3_all - repmat(passive_str_3_stim_baseline,1,size(passive_psth_stim_str_3_all,2));
%     passive_psth_stim_str_3_mean_baseline_sub = mean(passive_psth_stim_str_3_all_baseline_sub, 1)';
% 
%     % figure to compare to kernel
%     figure;
%     plot(passive_psth_stim_str_1_mean_baseline_sub(51:100));
%     hold on
%     plot(str_1_stim_move_rew_avg_coeff(1:50));
%     legend({'Passive stim psth', 'Task stim kernel'})
%     title([animal ' ' rec_day ' STIM Passive psth vs task kernel'])


    %% save for plot
%     
%     regression.str_1_var_expl{this_rec} = str_1_var_expl;
%     regression.str_1_stim_move_rew_avg_coeff_stim{this_rec} = str_1_stim_move_rew_avg_coeff(1:50);
%     regression.str_1_stim_move_rew_avg_coeff_stim_move{this_rec} = str_1_stim_move_rew_avg_coeff(51:200);
%     regression.str_1_stim_move_rew_avg_coeff_iti_move{this_rec} = str_1_stim_move_rew_avg_coeff(201:350);
%     regression.str_1_stim_move_rew_avg_coeff_reward{this_rec} = str_1_stim_move_rew_avg_coeff(351:400);
% 
%     regression.str_2_var_expl{this_rec} = str_2_var_expl;
%     regression.str_2_stim_move_rew_avg_coeff_stim{this_rec} = str_2_stim_move_rew_avg_coeff(1:50);
%     regression.str_2_stim_move_rew_avg_coeff_stim_move{this_rec} = str_2_stim_move_rew_avg_coeff(51:200);
%     regression.str_2_stim_move_rew_avg_coeff_iti_move{this_rec} = str_2_stim_move_rew_avg_coeff(201:350);
%     regression.str_2_stim_move_rew_avg_coeff_reward{this_rec} = str_2_stim_move_rew_avg_coeff(351:400);
% 
%     regression.str_3_var_expl{this_rec} = str_3_var_expl;
%     regression.str_3_stim_move_rew_avg_coeff_stim{this_rec} = str_3_stim_move_rew_avg_coeff(1:50);
%     regression.str_3_stim_move_rew_avg_coeff_stim_move{this_rec} = str_3_stim_move_rew_avg_coeff(51:200);
%     regression.str_3_stim_move_rew_avg_coeff_iti_move{this_rec} = str_3_stim_move_rew_avg_coeff(201:350);
%     regression.str_3_stim_move_rew_avg_coeff_reward{this_rec} = str_3_stim_move_rew_avg_coeff(351:400);
% 

    %%%% GO HERE %%%%%%%%%%
    regression_per_neuron.str_all_templates{this_rec} = str_all_templates;
    regression_per_neuron.bin_window_for_test = bin_window_for_test;
    regression_per_neuron.mean_neuron_stim_align_str_all_trace_baseline_sub{this_rec} = ...
        mean_neuron_stim_align_str_all_trace_baseline_sub;

    regression_per_neuron.neuron_str_all_move_rew_avg_coeff{this_rec} = ...
        neuron_str_all_move_rew_avg_coeff;
    regression_per_neuron.neuron_str_all_move_rew_predict_trace{this_rec} = ...
        neuron_str_all_move_rew_predict_trace;
    regression_per_neuron.neuron_str_all_var_expl_move_rew{this_rec} = ...
        neuron_str_all_var_expl.move_rew;

    regression_per_neuron.neuron_str_all_trace_residual{this_rec} = ...
        neuron_str_all_trace_residual;
    regression_per_neuron.mean_neuron_str_all_trace_residual_stim_align{this_rec} = ...
        mean_neuron_str_all_trace_residual_stim_align;

    regression_per_neuron.pre_stim{this_rec} = ...
        pre_stim;
    regression_per_neuron.post_stim{this_rec} = ...
        post_stim;
    regression_per_neuron.max_post_stim{this_rec} = ...
        max_post_stim;
    regression_per_neuron.str_all_sign_rank_test_p_val{this_rec} = ...
        str_all_sign_rank_test_p_val;


%     regression.stim_move_rew_predict_mua_trace_stim_align_mean{this_rec} = predict_mua_trace_stim_align_mean;
%     regression.stim_move_rew_psth_stim_str_all_mean_baseline_sub{this_rec} = psth_stim_str_1_mean_baseline_sub;
%     regression.stim_move_rew_passive_psth_stim_str_1_mean_baseline_sub{this_rec} = passive_psth_stim_str_1_mean_baseline_sub;

end

regression_per_neuron = rmfield(regression_per_neuron, 'str_all_templates');
for use_rec=3:length(recordings)

    this_rec = use_rec - 2;

    str_1_templates = task_ephys_data.str_1_templates{this_rec};
    str_2_templates = task_ephys_data.str_2_templates{this_rec};
    str_3_templates = task_ephys_data.str_3_templates{this_rec};
    str_all_templates = vertcat(str_1_templates, str_2_templates, str_3_templates);

    regression_per_neuron.str_all_templates{this_rec} = str_all_templates;

    regression_per_neuron.bin_edges{this_rec} = task_ephys_data.bin_edges{this_rec};
    regression_per_neuron.bin_centres{this_rec} = task_ephys_data.bin_centres{this_rec};
end

% fieldnames(regression_per_neuron)

fields_order = {'days_from_learning',
    'str_all_templates',
    'bin_edges',
    'bin_centres',
    'mean_neuron_stim_align_str_all_trace_baseline_sub',
    'neuron_str_all_move_rew_avg_coeff',
    'neuron_str_all_move_rew_predict_trace',
    'neuron_str_all_var_expl_move_rew',
    'neuron_str_all_trace_residual',
    'mean_neuron_str_all_trace_residual_stim_align',
    'bin_window_for_test',
    'pre_stim',
    'post_stim',
    'str_all_sign_rank_test_p_val'};

regression_per_neuron = orderfields(regression_per_neuron, fields_order);
% save out?  
save_name = [animal '_per_neuron_task_regression_data'];
save(save_name, "regression_per_neuron", "-v7.3");

%% load stuff

load([animal '_task_regression_data']);
load([animal '_per_neuron_task_regression_data']);

% LEFT HERE %%%%%

%% across days expl variance
stim_diff = struct;
for day_idx=1:length(regression.days_from_learning)
    stim_diff(day_idx).str_1 = regression.str_1_var_expl{day_idx}.stim_move_rew - regression.str_1_var_expl{day_idx}.move_rew;
    stim_diff(day_idx).str_2 = regression.str_2_var_expl{day_idx}.stim_move_rew - regression.str_2_var_expl{day_idx}.move_rew;
    stim_diff(day_idx).str_3 = regression.str_3_var_expl{day_idx}.stim_move_rew - regression.str_3_var_expl{day_idx}.move_rew;
end

figure;
plot(regression.days_from_learning, [stim_diff.str_1], '-o')
hold on;
plot(regression.days_from_learning, [stim_diff.str_2],  '-o')
hold on;
plot(regression.days_from_learning, [stim_diff.str_3],  '-o')
yline(0, 'LineWidth', 2)
xlabel('Days from learning')
legend({'Str 1', 'Str 2', 'Str 3', 'Threshold'})
title('Adding stim to regression diff in var expl')

%% str diff regions stim kernels
figure = tiledlayout('flow');
for day_idx=1:length(regression.days_from_learning)
    nexttile
    plot(regression.str_1_stim_move_rew_avg_coeff_stim{day_idx})
    hold on;
    plot(regression.str_2_stim_move_rew_avg_coeff_stim{day_idx})
    hold on;
    plot(regression.str_3_stim_move_rew_avg_coeff_stim{day_idx})
    legend({'Str 1', 'Str 2', 'Str 3'})
    title(['Day ' num2str(regression.days_from_learning(day_idx))])
end
sgtitle([animal 'Stim kernels across days']

%% plot across days with all kernels

animal = 'AM021';

load([animal '_task_regression_data']);

kernels_across_days_fig = figure('Position', get(0, 'Screensize'));
tiledlayout('flow')
sgtitle([animal ' Kernels across days'], 'FontSize', 14, 'FontWeight', 'bold')

kernel_stim = nexttile;
plot(vertcat(regression.stim_move_rew_avg_coeff_stim{regression.days_from_learning<0})', '--', 'linewidth', 2);
hold on
plot(vertcat(regression.stim_move_rew_avg_coeff_stim{regression.days_from_learning>=0})', 'linewidth', 2);
kernel_stim.ColorOrder = brewermap(length(regression.days_from_learning),"PuRd");
title('Stim')
xlim([0 50])

kernel_stim_move = nexttile;
plot(vertcat(regression.stim_move_rew_avg_coeff_stim_move{regression.days_from_learning<0})', '--', 'linewidth', 2);
hold on
plot(vertcat(regression.stim_move_rew_avg_coeff_stim_move{regression.days_from_learning>=0})', 'linewidth', 2);
kernel_stim_move.ColorOrder = brewermap(length(regression.days_from_learning),"Blues");
title('Stim move')
xlim([0 150])

kernel_iti_move = nexttile;
plot(vertcat(regression.stim_move_rew_avg_coeff_iti_move{regression.days_from_learning<0})', '--', 'linewidth', 2);
hold on
plot(vertcat(regression.stim_move_rew_avg_coeff_iti_move{regression.days_from_learning>=0})', 'linewidth', 2);
kernel_iti_move.ColorOrder = brewermap(length(regression.days_from_learning),"Greens");
title('ITI Move')
xlim([0 150])

kernel_reward = nexttile;
plot(vertcat(regression.stim_move_rew_avg_coeff_reward{regression.days_from_learning<0})', '--', 'linewidth', 2);
hold on
plot(vertcat(regression.stim_move_rew_avg_coeff_reward{regression.days_from_learning>=0})', 'linewidth', 2);
kernel_reward.ColorOrder = brewermap(length(regression.days_from_learning),"Purples");
title('Reward')
xlim([0 50])

kernels_across_days_fig_name = [animal '_kernels_across_days.tif'];
kernels_across_days_fig_path = fullfile(save_fig_path, kernels_across_days_fig_name);
saveas(kernels_across_days_fig, kernels_across_days_fig_path);

%% plot check kernel and passive + psth and predict

passive_vs_stim_kernel_fig = figure('Position', get(0, 'Screensize'));
tiledlayout('flow');
sgtitle([animal ' STIM Passive psth vs task stim kernel'], 'FontSize', 14, 'FontWeight', 'bold');
for this_rec=1:length(regression.days_from_learning)
    nexttile;
    passive_psth_stim_str_1_mean_baseline_sub = regression.all_passive_psth_stim_str_1_mean_baseline_sub{this_rec};
    plot(passive_psth_stim_str_1_mean_baseline_sub(51:100));
    hold on;
    plot(regression.all_avg_coeff_stim{this_rec});
    legend({'Passive stim psth', 'Task stim stim kernel'});
    title([rec_day ' day from learning: ' num2str(regression.days_from_learning(this_rec))]);
end

passive_vs_stim_kernel_fig_name = [animal '_passive_vs_stim_kernels.tif'];
passive_vs_stim_kernel_fig_path = fullfile(save_fig_path, passive_vs_stim_kernel_fig_name);
saveas(passive_vs_stim_kernel_fig, passive_vs_stim_kernel_fig_path);

% psth vs predict mua trace
psth_vs_predict_mua_fig = figure('Position', get(0, 'Screensize'));
tiledlayout('flow');
sgtitle([animal ' Task psth vs predicted mua trace'], 'FontSize', 14, 'FontWeight', 'bold');
for this_rec=1:length(regression.days_from_learning)
    nexttile;
    plot(regression.all_psth_stim_str_1_mean_baseline_sub{this_rec});
    hold on;
    plot(regression.all_predict_mua_trace_stim_align_mean{this_rec})
    legend({'Task stim psth', 'Predict mua trace'});
    title([rec_day ' day from learning: ' num2str(regression.days_from_learning(this_rec))]);
end

psth_vs_predict_mua_fig_name = [animal '_task_psth_vs_predict_mua.tif'];
psth_vs_predict_mua_fig_path = fullfile(save_fig_path, psth_vs_predict_mua_fig_name);
saveas(psth_vs_predict_mua_fig, psth_vs_predict_mua_fig_path);

%% TEST CORTEX on a day with stim response

%% - do cross-val regression from baseline sub mua trace
cortex_var_expl_baseline_sub = struct;
%% all

% append to get design matrix
design_matrix = horzcat(stim_regressor_matrix, ...
    move_regressor_matrix, ...
    reward_regressor_matrix); % timepoints x lags

[str_1_all_avg_coeff, str_1_all_predict_mua_trace, str_1_all_var_expl_baseline] = ephys.run_regression(cortex_mua_trace_baseline_sub, ...
    regressors, t_shifts, downsampled_time);


figure;
plot(cortex_mua_trace_baseline_sub);
% hold on;
% plot(smooth_cortex_mua_trace_baseline_sub);
hold on;
plot(str_1_all_predict_mua_trace)

figure;
plot(str_1_all_avg_coeff(1:50))
hold on;
plot(str_1_all_avg_coeff(51:200))
hold on;
plot(str_1_all_avg_coeff(201:250))
legend({'Stim', 'Move', 'Reward'})
title([rec_day ' Kernels for baseline sub mua trace str 1'])

% get ssms
cortex_var_expl_baseline_sub.all = 1 - (var(cortex_mua_trace_baseline_sub - str_1_all_predict_mua_trace)) / (var(cortex_mua_trace_baseline_sub));

%% TOY DATASET 

% toy_data_event_onset = randi([0 1], length(downsampled_time), 1);
% toy_data_kernel_1 = gausswin(11, 3)/sum(gausswin(11, 3));
% toy_data_event_trace_full = conv(toy_data_event_onset, toy_data_kernel_1, 'full');
% toy_data_event_trace = toy_data_event_trace_full((size(toy_data_kernel_1,1)-1)/2+1:...
%     end-(size(toy_data_kernel_1,1)-1)/2);
% 
% 
% toy_data_regressor = toy_data_event_onset;
% 
% toy_data_baseline = mean(toy_data_event_trace);
% toy_data_event_trace_baseline_sub = toy_data_event_trace - toy_data_baseline;
% % figure;
% % plot(toy_data_event_trace_baseline_sub)
% 
% % get idx for stim time shift
% timestep = mean(diff(downsampled_time));
% toy_data_t_shift = -0.5:timestep:0.5;
% toy_data_t_shift_idx = double(int32(toy_data_t_shift * 1/timestep));
% 
% % calculate design matrix for stim regressor 
% toy_data_regressor_matrix = lagmatrix(toy_data_regressor,toy_data_t_shift_idx);
% toy_data_regressor_matrix(isnan(toy_data_regressor_matrix)) = 0;
% 
% % append to get design matrix
% toy_dataset_design_matrix = horzcat(toy_data_regressor_matrix); % timepoints x lags
% 
% [toy_avg_coeff, toy_predict_mua_trace, toy_var_expl.all] = ephys.run_regression(toy_data_event_trace_baseline_sub, ...
%     toy_dataset_design_matrix, downsampled_time);
% 
% figure;
% plot(toy_data_event_trace);
% hold on;
% plot(toy_data_predict_mua_trace)
% 
% figure;
% plot(toy_data_avg_coeff)
% legend({'Test 1'})
% title('Toy dataset kernels')
% 
% figure; 
% plot(toy_data_kernel_1)
