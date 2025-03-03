%% Load experiments and save in struct

%% load bhv, save path and dataset
save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Users\Andrada-Maria_Marica\long_str_ctx_data';

load(fullfile(save_path, 'swr_bhv.mat'));

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

%% - go through each animal
for animal_idx=1:length(animals)
    task_ephys_animal = table;

    animal = animals{animal_idx};
    workflow_passive = {'lcr_passive'};
    recordings_passive = plab.find_recordings(animal, [], workflow_passive);
    workflow_task = {'stim_wheel_right*'};
    recordings_task = plab.find_recordings(animal, [], workflow_task);
    training_days = ismember({recordings_passive.day}, {recordings_task.day});
    train_rec_passive = recordings_passive(training_days);
    bhv_days = {train_rec_passive.day};

    % take out after
    days_from_learning = bhv.days_from_learning(strcmp(bhv.animal, animal));
    figure;
    tiledlayout('flow');
    sgtitle([animal ' resp cells'])
    % take out after

    for use_rec=1:length(recordings_task)

        rec_day = recordings_task(use_rec).day;
        rec_time = recordings_task(use_rec).recording{end};
        verbose = true;
        load_parts.behavior = true;
        load_parts.ephys = true;
        ap.load_recording

        %% Get striatum boundaries - Just skip missing depth ones

        AP_longstriatum_find_striatum_depth
        mua_length = 200;
        depth_group_edges = striatum_depth(1):mua_length:striatum_depth(end);
        try
            depth_group = discretize(spike_depths,depth_group_edges);
        catch ME
            warning(['Undefined str depth ' animal ' ' rec_day])
            %             depth_group = nan;
            task_ephys_animal.animal(use_rec) = {animal};
            task_ephys_animal.rec_day(use_rec) = {rec_day};
            continue
        end

        unit_depth_group = discretize(template_depths, depth_group_edges);

        %% single unit labels and trial values
        single_unit_idx = strcmp(template_qc_labels, 'singleunit');

        trial_stim_values = ones(length(n_trials), 1);

        %% cell type labels
        AP_longstriatum_classify_striatal_units
        str_tan_idx = striatum_celltypes.tan;
        str_fsi_idx = striatum_celltypes.fsi;
        str_msn_idx = striatum_celltypes.msn;

        %% psth
        [~,binned_spikes_stim_align] = ap.psth(spike_times_timelite, ...
            stimOn_times(1:n_trials),  ...
            depth_group);

        %% EDIT: unit responsivenes
        %% -- contra stim
        % sharp
        % create pre and post stim onsets
        bin_window_for_pre = 0.2;
        bin_window_for_post = 0.2;
        pre_stim_time = stimOn_times(1:n_trials) - [bin_window_for_pre 0];
        post_stim_time = stimOn_times(1:n_trials) + [0 bin_window_for_post];

        % group stim times into cell arrays before use arrayfun
        % example:
        use_align = arrayfun(@(x) stimOn_times(1:n_trials), ...
            unique(trial_stim_values),'uni',false);
        
        [event_psths,~] = ap.psth(spike_times_timelite,use_align,spike_templates);

        window_center = [-0.1,0.1];
        bin_size = 0.2;
        [~,event_spikes] = ap.psth(spike_times_timelite,use_align,spike_templates, ...
            'window',window_center,'bin_size',bin_size);

        % added so same code from passive works
        event_spikes = {event_spikes};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real_event_avg_spikes = arrayfun(@(x) squeeze(mean(diff(event_spikes{x}, [], 2), 1)), ...
            1:length(unique(trial_stim_values)),'uni',false);

        % get mean pre and post
        mean_pre_stim = arrayfun(@(stim_idx) ...
            squeeze(mean(event_spikes{stim_idx}(:,1, :), 1)), ...
            1:length(unique(trial_stim_values)), ...
            'uni', false);
        mean_post_stim = arrayfun(@(stim_idx) ...
            squeeze(mean(event_spikes{stim_idx}(:,2, :), 1)), ...
            1:length(unique(trial_stim_values)), ...
            'uni', false);
        std_post_stim = arrayfun(@(stim_idx) ...
            squeeze(std(event_spikes{stim_idx}(:,2, :), [], 1)), ...
            1:length(unique(trial_stim_values)), ...
            'uni', false);

        % shuffle test for responsiveness
        num_shuffles = 10000;
        shuffle_event_avg_spikes = cell(length(unique(trial_stim_values)), 1);
        for stim_idx=1:length(unique(trial_stim_values))
            stim_shuffle_event_avg_spikes = nan(length(find(good_templates)), num_shuffles);
            for shuffle=1:num_shuffles
                this_stim_event_spikes = event_spikes{stim_idx};
                stim_shuffle_event_avg_spikes(:, shuffle) = squeeze( ...
                    mean(diff(ap.shake(this_stim_event_spikes, 2), [], 2), 1));
            end
            shuffle_event_avg_spikes{stim_idx} = stim_shuffle_event_avg_spikes;
        end

        all_event_avg_spikes = arrayfun(@(x) ...
            horzcat(real_event_avg_spikes{x}, shuffle_event_avg_spikes{x}), ...
            1:length(unique(trial_stim_values)),'uni',false);

        unit_rank = cell(length(unique(trial_stim_values)), 1);
        unit_resp_p_value = cell(length(unique(trial_stim_values)), 1);
        for stim_idx=1:length(unique(trial_stim_values))
            [unit_rank{stim_idx}, ~] = arrayfun(@(template_idx) ...
                tiedrank(all_event_avg_spikes{stim_idx}(template_idx,:)), ...
                1:length(find(good_templates)), ...
                'uni',false);
            unit_resp_p_value{stim_idx} = arrayfun(@(x) ...
                1 - unit_rank{stim_idx}{x}(1) / length(unit_rank{stim_idx}{x}), ...
                1:length(find(good_templates)), ...
                'uni',false);
        end
        resp_units = arrayfun(@(stim_idx) [unit_resp_p_value{stim_idx}{:}] < 0.01, ...
            1:length(unique(trial_stim_values)),'uni',false);

        % for plot
        time_vector = -0.5:0.001:1;
        baseline_idx = time_vector >= -0.2 & time_vector <= 0;

        % Function to normalize and smooth
        gauss_win = gausswin(51, 3)';
        normalize_and_smooth = @(psth) ...
            filter(gauss_win, sum(gauss_win), psth - mean(psth(baseline_idx), 2), [], 2);

        % Apply normalization and smoothing to all units and stim IDs
        smooth_event_psths = arrayfun(@(stim) ...
            cell2mat(arrayfun(@(unit) ...
            normalize_and_smooth(squeeze(event_psths(unit, :, 1))), ...
            (1:length(unique(spike_templates)))', ...
            'UniformOutput', false)), ...
            1:length(unique(trial_stim_values)), ...
            'UniformOutput', false);

        time_for_plot = -0.2:0.001:0.5;
        time_for_plot_indices = find(-0.5:0.001:1 >= time_for_plot(1) & ...
            -0.5:0.001:1 <= time_for_plot(end));

        for stim_idx=1:length(unique(trial_stim_values))

            % Sort by response
            time_range_sort = 0:0.001:0.2;
            time_indices_sort = find(-0.5:0.001:1 >= time_range_sort(1) & ...
                -0.5:0.001:1 <= time_range_sort(end));
            avg_response = mean(smooth_event_psths{stim_idx}(resp_units{stim_idx}, time_indices_sort), 2);
            [~, sorted_indices] = sort(avg_response, 'descend');

            % psth to plot
            psths_plot = smooth_event_psths{stim_idx}(resp_units{stim_idx}, time_for_plot_indices);
            sorted_psths_plot = psths_plot(sorted_indices, :);

            nexttile;
            imagesc(time_for_plot, [], sorted_psths_plot)
            colormap(AP_colormap('BWR'))
            max_abs = max(abs(smooth_event_psths{stim_idx}(resp_units{stim_idx}, :)), [], 'all');
            if isempty(max_abs)
                continue
            end
            clim([-max_abs, max_abs]);
            %         clim([min(smooth_event_psths(resp_units{1}, :, 1), [], 'all'), ...
            %             max(smooth_event_psths(resp_units{1}, :, 1), [], 'all')]);
            title([rec_day ' LD ' num2str(days_from_learning(use_rec))])
            drawnow;
        end

        %% save data in mouse table
        task_ephys_animal.animal(use_rec) = {animal};
        task_ephys_animal.rec_day(use_rec) = {rec_day};

        task_ephys_animal.depth_group_edges(use_rec) = {depth_group_edges};
        task_ephys_animal.unit_depth_group(use_rec) = {unit_depth_group};

        task_ephys_animal.binned_spikes_stim_align(use_rec) = {binned_spikes_stim_align};

        task_ephys_animal.bin_window_for_pre(use_rec) = {bin_window_for_pre};
        task_ephys_animal.bin_window_for_post(use_rec) = {bin_window_for_post};
        task_ephys_animal.unit_resp_p_value(use_rec) = {unit_resp_p_value};

        task_ephys_animal.mean_pre_stim(use_rec) = {mean_pre_stim};
        task_ephys_animal.mean_post_stim(use_rec) = {mean_post_stim};
        task_ephys_animal.std_post_stim(use_rec) = {std_post_stim};

        task_ephys_animal.unit_smooth_event_psths(use_rec) = {smooth_event_psths}; 

        task_ephys_animal.single_unit_idx(use_rec) = {single_unit_idx}; 

        task_ephys_animal.str_tan_idx(use_rec) = {str_tan_idx}; 
        task_ephys_animal.str_fsi_idx(use_rec) = {str_fsi_idx}; 
        task_ephys_animal.str_msn_idx(use_rec) = {str_msn_idx}; 

        disp(['Done day ' num2str(use_rec)])

    end
    all_task_ephys_save_cell{animal_idx} = task_ephys_animal;
    disp(['Done ' animal]);
end

task_ephys = vertcat(all_task_ephys_save_cell{:});
save_name = fullfile(save_path, 'task_ephys');
save(save_name, "task_ephys", "-v7.3");

