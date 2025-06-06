%% Load dataset

load('all_swr_bhv_data.mat');

animals = ["AM011"; "AM012"; "AM014"; "AM015"; "AM016"; ...
    "AM017"; "AM018"; "AM019"; "AM021"; "AM022"; "AM026"]; %; "AM029"];

for animal_idx=1:length(animals)
    animal = animals{animal_idx};
    disp(['Start ' animal])
    workflow_passive = {'lcr_passive'};
    recordings_passive = plab.find_recordings(animal, [], workflow_passive);
    workflow_task = {'stim_wheel_right*'};
    recordings_task = plab.find_recordings(animal, [], workflow_task);
    training_days = ismember({recordings_passive.day}, {recordings_task.day});
    train_rec_passive = recordings_passive(training_days);

    bhv_days = swr_bhv_data.bhv_days{animal_idx};
    days_from_learning = swr_bhv_data.days_from_learning{animal_idx};

    ctx_maps_to_str = table;
    for use_rec=1:length(train_rec_passive)
        rec_day = train_rec_passive(use_rec).day;
        rec_time = train_rec_passive(use_rec).recording{end};
        ap.load_recording;

        this_day_from_learning = days_from_learning(ismember(bhv_days, rec_day));

        ctx_maps_to_str.rec_day(use_rec) = {rec_day};
        ctx_maps_to_str.days_from_learning(use_rec) = {this_day_from_learning};

        %% Get striatum MUA in bins of set lengths
        % Set MUA bin length (microns)
        mua_length = 200;

        % Get striatum start = lowest unit density, end = end of probe
        % (this probably isn't robust/optimal)
        unit_density_bins = 0:100:3840;
        unit_density = histcounts(template_depths,unit_density_bins);
        [~,unit_density_min_bottom_idx] = min(fliplr(unit_density));
        unit_density_min_idx = length(unit_density_bins) - unit_density_min_bottom_idx;
        template_depths_sorted = sort(template_depths);
        str_start =  template_depths_sorted(find(template_depths_sorted >= ...
            unit_density_bins(unit_density_min_idx+1),1));
        str_end = max(channel_positions(:,2));

%         %%% ANDRADA
%         [sorted_template_depths, sorted_template_depth_idx] = sort(template_depths);
% 
%         % find str start where unit depth distr not linear anymore
%         idx_str_start = ischange(sorted_template_depths, 'linear','MaxNumChanges',1);
% 
%         %         test_x_axis = 1:length(sorted_template_depths);
%         %         figure;
%         %         plot(sorted_template_depths, 'o')
%         %         hold on
%         %         plot(test_x_axis(test_a), sorted_template_depths(test_a), '*')
%         %         title(rec_day)
% 
%         str_start = sorted_template_depths(idx_str_start);
%         str_end = sorted_template_depths(end);
%         str_depth = [str_start,str_end];
%         %%%%%%%%%%%

        % Discretize spikes by depth
        depth_group_edges = str_start:mua_length:str_end;
        if length(depth_group_edges)<2
            ctx_maps_to_str.depth_group_edges(use_rec) = {depth_group_edges};
            ctx_maps_to_str.no_depths(use_rec) = {0};
            warning([animal ' ' rec_day ' bad striatum edges'])
            continue
        end
        depth_group = discretize(spike_depths,depth_group_edges);

        % Get time bins corresponding to widefield frame exposures
        % (skip the beginning and end of the recording to avoid artifacts)
        sample_rate = (1/mean(diff(wf_t)));
        skip_seconds = 60;
        time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

        % Bin spikes in depth and time
        binned_spikes = zeros(max(depth_group),length(time_bins)-1);
        for curr_depth = 1:max(depth_group)
            curr_spike_times = spike_times_timelite(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end

        %% Regress ctx fluorescence to striatal MUA

        % Normalize MUA by standard devation
        % (to have kernels in comparable units)
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);

        % Set parameters for regression
        % (widefield components to use - too many can overload memory or add noise,
        % too few gives a less detailed map)
        use_svs = 1:500;
        % (frame lags to use: for maps rather than prediction, can just do zero-lag
        % regression with no temporal component)
        kernel_t = [0,0];
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        % (ridge regression normalizing factor lambda: too low gives just noise,
        % higher removes map detail)
        lambda = 5;
        % (whether to z-score regressors or signals: not here, since everything
        % already normalized)
        zs = [false,false];
        % (cross validation: not needed for maps, only for prediction)
        cvfold = 1;
        % (whether to use (yes) and return as ouput (no) a constant in regression)
        use_constant = true;
        return_constant = false;

        % Interpolate the widefield to the center of each time bin to match the
        % MUA timepoints
        fVdf_deconv_resample = interp1(wf_t,wf_V(use_svs,:)',time_bin_centers)';

        % Do regression from cortex to spikes
        [cortex_kernel,predicted_spikes,explained_var] = ...
            ap.regresskernel(fVdf_deconv_resample, ...
            binned_spikes_std,kernel_frames, ...
            lambda,zs,cvfold,return_constant,use_constant);

        % (note: cortex_kernel is components x time x MUA - so in this example is
        % 500 used components x 1 time point x 12 MUA signals)

%         % Convert ctx kernal V into pixels, plot
         cortex_kernel_px = squeeze(plab.wf.svd2px(wf_U(:,:,use_svs),cortex_kernel));
%         ap.imscroll(cortex_kernel_px);
%         axis image;
%         clim(max(abs(clim)).*[-1,1]*0.7);
%         ap.wf_draw('ccf','k');
%         colormap(ap.colormap('PWG'));


        %% ADD RAW TRACE ALIGNED TO STIM
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        contra_stim = 90;
        centre_stim = 0;
        ipsi_stim = -90;

        trial_stim_values = vertcat(trial_events.values.TrialStimX);
        trial_stim_values = trial_stim_values(1:length(stimOn_times));


        % stim aligned wheel move
        % create matrix of times for stim onset
        timestep = 0.01;
        start_time = -0.5;
        end_time = 1;  %%%%%%% can change to [0, 0.5]
        timevec = start_time:timestep:end_time;
        stim_frame = (-start_time)*(1/timestep)+1;
        time_stimulus = stimOn_times+timevec;
        % stim aligned wheel move
        t = timelite.timestamps;
        stim_wheel_move = interp1(t,+wheel_move,time_stimulus);
        no_move_trials = sum(stim_wheel_move(:,stim_frame:end),2)==0;

        contra_good_trials = (trial_stim_values == contra_stim) & no_move_trials;
        centre_good_trials = (trial_stim_values == centre_stim) & no_move_trials;
        ipsi_good_trials = (trial_stim_values == ipsi_stim) & no_move_trials;

        psth_opts.window = [-0.5,1];
        psth_opts.bin_size = 0.001;

        binned_spikes_stim_align = ap.ephys_psth(spike_times_timelite(~isnan(depth_group)), ...
            num2cell(stimOn_times),  ...
            depth_group(~isnan(depth_group)));

        % around stim time
        bin_edges = psth_opts.window(1):psth_opts.bin_size:psth_opts.window(2);
        bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;

%         around_stim_time = stimOn_times + bin_edges;
%         binned_spikes_stim_align = zeros([max(depth_group) [size(around_stim_time)]]);
%         for curr_depth = 1:max(depth_group)
%             binned_spikes_stim_align(curr_depth,:,:) = interp1(time_bin_centers, binned_spikes(curr_depth, :), around_stim_time);
%         end

        % contra stim
        contra_no_move_binned_spikes_stim_align = binned_spikes_stim_align(:,:, contra_good_trials);
        contra_no_move_psth_stim_align = squeeze(nanmean(contra_no_move_binned_spikes_stim_align, 3));
        % - normalize
        contra_stim_baseline = mean(contra_no_move_psth_stim_align(:,bin_centres>-0.2 & bin_centres<0), 2);
        contra_norm_no_move_psth_stim_align = (contra_no_move_psth_stim_align - contra_stim_baseline) ...
                ./ (contra_stim_baseline + std(contra_stim_baseline));

        % centre stim
        centre_no_move_binned_spikes_stim_align = binned_spikes_stim_align(:,:, centre_good_trials);
        centre_no_move_psth_stim_align = squeeze(nanmean(centre_no_move_binned_spikes_stim_align, 3));
        % - normalize
        centre_stim_baseline = mean(centre_no_move_psth_stim_align(:,bin_centres>-0.2 & bin_centres<0), 2);
        centre_norm_no_move_psth_stim_align = (centre_no_move_psth_stim_align - centre_stim_baseline) ...
                ./ (centre_stim_baseline + std(centre_stim_baseline));
        
        % ipsi stim
        ipsi_no_move_binned_spikes_stim_align = binned_spikes_stim_align(:,:, ipsi_good_trials);
        ipsi_no_move_psth_stim_align = squeeze(nanmean(ipsi_no_move_binned_spikes_stim_align, 3));
        % - normalize
        ipsi_stim_baseline = mean(ipsi_no_move_psth_stim_align(:,bin_centres>-0.2 & bin_centres<0), 2);
        ipsi_norm_no_move_psth_stim_align = (ipsi_no_move_psth_stim_align - ipsi_stim_baseline) ...
                ./ (ipsi_stim_baseline + std(ipsi_stim_baseline));

%         %% responsive cells
%         % contra
%         % create pre and post stim onsets
% 
%         non_nan_depth_group = depth_group(~isnan(depth_group));
%         unique_depth_group = unique(non_nan_depth_group);
% 
% %         bin_window_for_sharp = 0.1;
% %         pre_stim_time = stimOn_times - [bin_window_for_sharp 0];
% %         post_stim_time = stimOn_times + [0.05 bin_window_for_sharp+0.05];
% % 
% %         unit_spikes_small_pre_stim = nan(length(good_templates_idx), length(find(contra_good_trials)));
% %         unit_spikes_small_post_stim = nan(length(good_templates_idx), length(find(contra_good_trials)));
% %         contra_sharp_p_units = nan(length(good_templates_idx), 1);
%         for curr_depth=1:length(unique_depth_group)
%             depth_spike_times = spike_times_timelite(non_nan_depth_group == curr_depth);
% 
%             % spike counts binned pre/post stim
%             this_stim_time = pre_stim_time(contra_good_trials, :);
%             depth_spikes_small_pre_stim(depth_idx, :)  = cell2mat(arrayfun(@(trial_id) histcounts(depth_spike_times, this_stim_time(trial_id,:))', ...
%                 1:size(this_stim_time, 1), 'UniformOutput',false))' / bin_window_for_sharp;
% 
%             this_stim_time = post_stim_time(contra_good_trials, :);
%             depth_spikes_small_post_stim(depth_idx, :) = cell2mat(arrayfun(@(trial_id) histcounts(depth_spike_times, this_stim_time(trial_id,:))', ...
%                 1:size(this_stim_time, 1), 'UniformOutput',false))' / bin_window_for_sharp;
% 
%             % signed rank test
%             contra_sharp_p_depths(depth_idx) = signrank(depth_spikes_small_post_stim(depth_idx, :), depth_spikes_small_pre_stim(depth_idx, :), 'tail', 'right');
%         end

        %% save everything
        ctx_maps_to_str.no_depths(use_rec) = {max(depth_group)};
        ctx_maps_to_str.depth_group_edges(use_rec) = {depth_group_edges};
        %         ctx_maps_to_str.use_svs(use_rec) = {use_svs};
        %         ctx_maps_to_str.cortex_kernel(use_rec) = {cortex_kernel};
        ctx_maps_to_str.cortex_kernel_px(use_rec) = {cortex_kernel_px};
        ctx_maps_to_str.explained_var(use_rec) = {explained_var.total};
        
        ctx_maps_to_str.bin_edges(use_rec) = {bin_edges};
        ctx_maps_to_str.bin_centres(use_rec) = {bin_centres};
        ctx_maps_to_str.contra_no_move_binned_spikes_stim_align(use_rec) = {contra_no_move_binned_spikes_stim_align};
        ctx_maps_to_str.contra_no_move_psth_stim_align(use_rec) = {contra_no_move_psth_stim_align};
        ctx_maps_to_str.contra_norm_no_move_psth_stim_align(use_rec) = {contra_norm_no_move_psth_stim_align};
        
        ctx_maps_to_str.centre_no_move_binned_spikes_stim_align(use_rec) = {centre_no_move_binned_spikes_stim_align};
        ctx_maps_to_str.centre_no_move_psth_stim_align(use_rec) = {centre_no_move_psth_stim_align};
        ctx_maps_to_str.centre_norm_no_move_psth_stim_align(use_rec) = {centre_norm_no_move_psth_stim_align};
        
        ctx_maps_to_str.ipsi_no_move_binned_spikes_stim_align(use_rec) = {ipsi_no_move_binned_spikes_stim_align};
        ctx_maps_to_str.ipsi_no_move_psth_stim_align(use_rec) = {ipsi_no_move_psth_stim_align};
        ctx_maps_to_str.ipsi_norm_no_move_psth_stim_align(use_rec) = {ipsi_norm_no_move_psth_stim_align};

        disp(['Done day ' num2str(use_rec)])
        
    end
    all_ctx_maps_to_str_cell{animal_idx} = ctx_maps_to_str;
    disp(['Done ' animal])
end

all_ctx_maps_to_str = table();
all_ctx_maps_to_str.animals = animals;
all_ctx_maps_to_str.recording_data = all_ctx_maps_to_str_cell';

save_name = 'new_all_ctx_maps_and_passive_str_ephys_data';
save(save_name, "all_ctx_maps_to_str", "-v7.3");

