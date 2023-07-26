%% Load experiment


animal = 'AP003';

% use_workflow = {'lcr_passive'};
use_workflow = {'lcr_passive_fullscreen'};
% use_workflow = {'lcr_passive','lcr_passive_fullscreen'};
% use_workflow = {'stim_wheel_right_stage1','stim_wheel_right_stage2'};
% use_workflow = 'sparse_noise';

recordings = ap.find_recordings(animal,use_workflow);

% use_rec = 1;
rec_day = '2023-06-07';
use_rec = strcmp(rec_day,{recordings.day});
% use_rec = length(recordings)-1;

rec_day = recordings(use_rec).day;
rec_time = recordings(use_rec).protocol{end};

verbose = true;

ap.load_experiment

qMetrics_path = plab.locations.make_server_filename(animal,rec_day,[],'ephys','qMetrics');
if isfolder(qMetrics_path)
    qMetricsExist = 0;
    qMetrics_probe_folders = dir(qMetrics_path);
    qMetrics_probe_folders = qMetrics_probe_folders(3:end);

    if  isempty(find([qMetrics_probe_folders.isdir]))
        % Load Julie's quality metrics
        qMetricsExist = ~isempty(dir(fullfile(qMetrics_path, 'qMetric*.mat'))) || ~isempty(dir(fullfile(qMetrics_path, 'templates._bc_qMetrics.parquet')));

        if qMetricsExist
            [param, qMetric] = bc_loadSavedMetrics(qMetrics_path);
%             unitType = bc_getQualityUnitType(param, qMetric);
            load(fullfile(qMetrics_path, 'am_bc_unit_type.mat'))
            disp(['Loaded qMetrics']);
           
            % Define good units from labels
            bc_good_templates_idx_mask = am_bc_units == 1 | am_bc_units == 2;
            bc_good_templates_idx = find(bc_good_templates_idx_mask)';
            bc_good_templates_identity = bc_good_templates_idx - 1;

            bc_axons_idx_mask = am_bc_units == 3;
            bc_axons_idx = find(bc_axons_idx_mask)';
            bc_axons_identity = bc_axons_idx - 1;

            bc_noise_idx_mask = am_bc_units == 0;
            bc_noise_idx = find(bc_noise_idx_mask)';
            bc_noise_identity = bc_noise_idx - 1;

            bc_inc_axons_idx_mask = bc_good_templates_idx_mask | bc_axons_idx_mask;
            bc_inc_axons_idx = find(bc_inc_axons_idx_mask)';
            bc_inc_axons_identity = bc_inc_axons_idx - 1;

            disp(['Bombcell found ' num2str(length(bc_good_templates_idx)) ' good units']);

        else
            qMetrics_probe_paths = strcat(qMetrics_path, filesep, qMetrics_probe_folders);

            for probe_idx=1:length(qMetrics_probe_paths)

                qMetrics_probe_path = qMetrics_probe_paths{probe_idx};

                % Load Julie's quality metrics
                qMetricsExist = ~isempty(dir(fullfile(qMetrics_probe_path, 'qMetric*.mat'))) || ~isempty(dir(fullfile(qMetrics_probe_path, 'templates._bc_qMetrics.parquet')));
                if qMetricsExist
                    [param, qMetric] = bc_loadSavedMetrics(qMetrics_probe_path);
                    %             unitType = bc_getQualityUnitType(param, qMetric);
                    load(fullfile(qMetrics_probe_path, 'am_bc_unit_type.mat'))
                    disp(['Loaded qMetrics']);

                    % Define good units from labels
                    bc_good_templates_idx_mask = am_bc_units == 1 | am_bc_units == 2;
                    bc_good_templates_idx = find(bc_good_templates_idx_mask)';
                    bc_good_templates_identity = bc_good_templates_idx - 1;

                    bc_axons_idx_mask = am_bc_units == 3;
                    bc_axons_idx = find(bc_axons_idx_mask)';
                    bc_axons_identity = bc_axons_idx - 1;

                    bc_noise_idx_mask = am_bc_units == 0;
                    bc_noise_idx = find(bc_noise_idx_mask)';
                    bc_noise_identity = bc_noise_idx - 1;

                    bc_inc_axons_idx_mask = bc_good_templates_idx_mask | bc_axons_idx_mask;
                    bc_inc_axons_idx = find(bc_inc_axons_idx_mask)';
                    bc_inc_axons_identity = bc_inc_axons_idx - 1;

                    disp(['Bombcell found ' num2str(length(bc_good_templates_idx)) ' good units']);

                end
            end
        end
    end
end

% Throw out all non-good template data
templates = templates(bc_good_templates_idx,:,:);
template_depths = template_depths(bc_good_templates_idx);
waveforms = waveforms(bc_good_templates_idx,:);
templateDuration = templateDuration(bc_good_templates_idx);
templateDuration_us = templateDuration_us(bc_good_templates_idx);

% Throw out all non-good spike data
good_spike_idx = ismember(spike_templates_0idx,bc_good_templates_identity);
spike_times_openephys = spike_times_openephys(good_spike_idx);
spike_templates_0idx = spike_templates_0idx(good_spike_idx);
template_amplitudes = template_amplitudes(good_spike_idx);
spike_depths = spike_depths(good_spike_idx);
spike_times_timeline = spike_times_timeline(good_spike_idx);

% Rename the spike templates according to the remaining templates
% (and make 1-indexed from 0-indexed)
new_spike_idx = nan(max(spike_templates_0idx)+1,1);
new_spike_idx(bc_good_templates_idx) = 1:length(bc_good_templates_idx);
spike_templates = new_spike_idx(spike_templates_0idx+1);


%% plot

trial_stim_values = vertcat(trial_events.values.TrialStimX);
% stimOn_times = photodiode_times(1:2:end);

AP_cellraster(stimOn_times,trialStim_values);

%% psth 
%% - multiunit

% create time vector around all stim onsets
bin_window = 0.1;
timevec = -0.5:bin_window:2;
around_stim_time = stimOn_times + timevec; 

% possible stims
possible_stim = unique(trial_stim_values);

%% -- hippocampus
hipp_depth = [2300 3000];
hipp_spikes = spike_depths>hipp_depth(1) & spike_depths<hipp_depth(2);
hipp_spike_times = spike_times_timeline(hipp_spikes);

% spike counts binned for each stim
hipp_spikes_in_stim_time = nan(size(around_stim_time, 1)/length(possible_stim), size(around_stim_time, 2)-1, length(possible_stim));
for stim_idx=1:length(possible_stim)
    this_stim = possible_stim(stim_idx);
    this_stim_time = around_stim_time(trial_stim_values == this_stim, :);
    % transpose to get the right shape
    hipp_spikes_in_stim_time(:, :, stim_idx) = cell2mat(arrayfun(@(trial_id) histcounts(hipp_spike_times, this_stim_time(trial_id,:))', ...
        1:size(this_stim_time, 1), 'UniformOutput',false))';
end

figure;
title('Hippocampus psth')
plot(mean(hipp_spikes_in_stim_time(:, :, 1)))
hold on;
plot(mean(hipp_spikes_in_stim_time(:, :, 2)))
hold on; 
plot(mean(hipp_spikes_in_stim_time(:, :, 3)))


%% -- visual cortex
vis_cortex_depth = [1500 2300];
vis_cortex_spikes = spike_depths>vis_cortex_depth(1) & spike_depths<vis_cortex_depth(2);
vis_cortex_spike_times = spike_times_timeline(vis_cortex_spikes);

% spike counts binned for each stim
vis_cortex_spikes_in_stim_time = nan(size(around_stim_time, 1)/length(possible_stim), size(around_stim_time, 2)-1, length(possible_stim));

for stim_idx=1:length(possible_stim)
    this_stim = possible_stim(stim_idx);
    this_stim_time = around_stim_time(trial_stim_values == this_stim, :);
    % transpose to get the right shape
    vis_cortex_spikes_in_stim_time(:, :, stim_idx) = cell2mat(arrayfun(@(trial_id) histcounts(vis_cortex_spike_times, this_stim_time(trial_id,:))', ...
        1:size(this_stim_time, 1), 'UniformOutput',false))';
end

figure;
title('Vis cortex psth')
plot(mean(vis_cortex_spikes_in_stim_time(:, :, 1)))
hold on;
plot(mean(vis_cortex_spikes_in_stim_time(:, :, 2)))
hold on; 
plot(mean(vis_cortex_spikes_in_stim_time(:, :, 3)))

%% - single unit

% create time vector around all stim onsets
bin_window = 0.1;
timevec = -0.5:bin_window:2;
around_stim_time = stimOn_times + timevec; 

% possible stims
possible_stim = unique(trial_stim_values);

%% -- visual cortex

good_unit_idx = 80;

good_vis_unit = (spike_templates == good_unit_idx) & vis_cortex_spikes;
good_vis_unit_spike_times = spike_times_timeline(good_vis_unit);

% spike counts binned for each stim
good_vis_unit_spikes_in_stim_time = nan(length(trial_stim_values)/length(possible_stim), size(around_stim_time, 2)-1, length(possible_stim));

for stim_idx=1:length(possible_stim)
    this_stim = possible_stim(stim_idx);
    this_stim_time = around_stim_time(trial_stim_values == this_stim, :);
    % transpose to get the right shape
    good_vis_unit_spikes_in_stim_time(:, :, stim_idx) = cell2mat(arrayfun(@(trial_id) histcounts(good_vis_unit_spike_times, this_stim_time(trial_id,:))', ...
        1:size(this_stim_time, 1), 'UniformOutput',false))';
end

figure;
title('Vis cortex psth')
plot(mean(good_vis_unit_spikes_in_stim_time(:, :, 1)))
hold on;
plot(mean(good_vis_unit_spikes_in_stim_time(:, :, 2)))
hold on; 
plot(mean(good_vis_unit_spikes_in_stim_time(:, :, 3)))

%% unit responsivenes 

% create pre and post stim onsets
pre_stim_time = stimOn_times - [0.1 0];
post_stim_time = stimOn_times + [0 0.1];

% possible stims
possible_stim = unique(trial_stim_values);
contra_stim = possible_stim(end);

unit_spikes_pre_stim = nan(length(bc_good_templates_idx), length(trial_stim_values)/length(possible_stim));
unit_spikes_post_stim = nan(length(bc_good_templates_idx), length(trial_stim_values)/length(possible_stim));
p_units = nan(length(bc_good_templates_idx), 1);
for unit_idx=1:length(bc_good_templates_idx)
    
    unit_spikes = (spike_templates == unit_idx) & vis_cortex_spikes;

    if sum(unit_spikes) == 0
        continue
    end

    unit_spike_times = spike_times_timeline(unit_spikes);

    % spike counts binned pre/post stim
    this_stim_time = pre_stim_time(trial_stim_values == contra_stim, :);
    unit_spikes_pre_stim(unit_idx, :)  = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, this_stim_time(trial_id,:))', ...
        1:size(this_stim_time, 1), 'UniformOutput',false))';

    this_stim_time = post_stim_time(trial_stim_values == contra_stim, :);
    unit_spikes_post_stim(unit_idx, :) = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, this_stim_time(trial_id,:))', ...
        1:size(this_stim_time, 1), 'UniformOutput',false))';

    % signed rank test
    p_units(unit_idx) = signrank(unit_spikes_pre_stim(unit_idx, :), unit_spikes_post_stim(unit_idx, :));
end

responsive_units = find(p_units<0.05);
unresponsive_units = find(p_units>0.05);

unit_spikes_diff_stim = mean(unit_spikes_post_stim - unit_spikes_pre_stim, 2);
[sorted_unit_spikes_diff_stim, sorted_units] = sort(unit_spikes_diff_stim, 'descend');

responsive_unit_spikes_diff_stim = mean(unit_spikes_post_stim(responsive_units) - unit_spikes_pre_stim(responsive_units), 2);
[~, sorted_responsive_units_idx] = sort(responsive_unit_spikes_diff_stim, 'descend');
sorted_responsive_units = responsive_units(sorted_responsive_units_idx);

unresponsive_unit_spikes_diff_stim = mean(unit_spikes_post_stim(unresponsive_units) - unit_spikes_pre_stim(unresponsive_units), 2);
[~, sorted_unresponsive_units_idx] = sort(unresponsive_unit_spikes_diff_stim, 'descend');
sorted_unresponsive_units = unresponsive_units(sorted_unresponsive_units_idx);

%% raw data
% spike counts binned for each stim
% create time vector around all stim onsets
bin_window = 0.1;
timevec = -0.5:bin_window:2;
around_stim_time = stimOn_times + timevec; 

all_spikes_in_stim_time = nan(length(bc_good_templates_idx), length(trial_stim_values)/length(possible_stim), size(around_stim_time, 2)-1);

for unit_idx=1:length(bc_good_templates_idx)
   
    unit_spikes = (spike_templates == unit_idx) & vis_cortex_spikes;

    if sum(unit_spikes) == 0
        continue
    end

    unit_spike_times = spike_times_timeline(unit_spikes);

    this_stim_time = around_stim_time(trial_stim_values == contra_stim, :);

    all_spikes_in_stim_time(unit_idx, :, :) = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, this_stim_time(trial_id,:))', ...
        1:size(this_stim_time, 1), 'UniformOutput',false))';
end

mean_all_spikes_in_stim_time = squeeze(mean(all_spikes_in_stim_time, 2));

sorted_mean_all_spikes_in_stim_time = mean_all_spikes_in_stim_time(sorted_units, :);

responsive_mean_all_spikes_in_stim_time = mean_all_spikes_in_stim_time(sorted_responsive_units, :);

unresponsive_mean_all_spikes_in_stim_time = mean_all_spikes_in_stim_time(sorted_unresponsive_units, :);


% raw data
% all cells
figure;
imagesc(sorted_mean_all_spikes_in_stim_time);
title('All cells')

% responsive cells
figure;
imagesc(responsive_mean_all_spikes_in_stim_time);
title('Responsive cells')

% unresponsive cells
figure;
imagesc(unresponsive_mean_all_spikes_in_stim_time);
title('Unresponsive cells')


%% normalized data

% create pre and post stim onsets
baseline_stim_time = stimOn_times - [0.5 0];

% possible stims
possible_stim = unique(trial_stim_values);
contra_stim = possible_stim(end);

unit_spikes_baseline_stim = nan(length(bc_good_templates_idx), length(trial_stim_values)/length(possible_stim));
for unit_idx=1:length(bc_good_templates_idx)
    
    unit_spikes = (spike_templates == unit_idx) & vis_cortex_spikes;

    if sum(unit_spikes) == 0
        continue
    end

    unit_spike_times = spike_times_timeline(unit_spikes);

    % spike counts binned pre/post stim
    this_stim_time = baseline_stim_time(trial_stim_values == contra_stim, :);
    unit_spikes_baseline_stim(unit_idx, :)  = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, this_stim_time(trial_id,:))', ...
        1:size(this_stim_time, 1), 'UniformOutput',false))';

end

%%%%%%%%%%% LEFT HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseline_all = mean(unit_spikes_baseline_stim(sorted_units), 2);

% normalized data
% all cells
figure;
imagesc(norm_sorted_mean_all_spikes_in_stim_time);
title('Norm All cells')

% responsive cells
figure;
imagesc(norm_responsive_mean_all_spikes_in_stim_time);
title('Norm Responsive cells')

% unresponsive cells
figure;
imagesc(norm_unresponsive_mean_all_spikes_in_stim_time);
title('Norm Unresponsive cells')