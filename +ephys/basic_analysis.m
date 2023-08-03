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

AP_cellraster(stimOn_times,trial_stim_values);

%% psth 
%% - multiunit

% create time vector around all stim onsets
bin_window = 0.1;
bin_edges = -0.5:bin_window:2;
around_stim_time = stimOn_times + bin_edges; 

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
        1:size(this_stim_time, 1), 'UniformOutput',false))' / bin_window;
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
        1:size(this_stim_time, 1), 'UniformOutput',false))' / bin_window;
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
bin_window = 0.05;
bin_edges = -0.5:bin_window:2;
around_stim_time = stimOn_times + bin_edges; 

% possible stims
possible_stim = unique(trial_stim_values);

% for plot
bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;

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
        1:size(this_stim_time, 1), 'UniformOutput',false))' / bin_window;
end

figure;
title('Vis cortex psth')
plot(bin_centres, mean(good_vis_unit_spikes_in_stim_time(:, :, 1)))
hold on;
plot(bin_centres, mean(good_vis_unit_spikes_in_stim_time(:, :, 2)))
hold on; 
plot(bin_centres, mean(good_vis_unit_spikes_in_stim_time(:, :, 3)))
xline(0.1)

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
    
    unit_spikes = spike_templates == unit_idx; 

    unit_spike_times = spike_times_timeline(unit_spikes);

    % spike counts binned pre/post stim
    this_stim_time = pre_stim_time(trial_stim_values == contra_stim, :);
    unit_spikes_pre_stim(unit_idx, :)  = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, this_stim_time(trial_id,:))', ...
        1:size(this_stim_time, 1), 'UniformOutput',false))' / bin_window;

    this_stim_time = post_stim_time(trial_stim_values == contra_stim, :);
    unit_spikes_post_stim(unit_idx, :) = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, this_stim_time(trial_id,:))', ...
        1:size(this_stim_time, 1), 'UniformOutput',false))' / bin_window;

    % signed rank test
    p_units(unit_idx) = signrank(unit_spikes_pre_stim(unit_idx, :), unit_spikes_post_stim(unit_idx, :));
end

responsive_units = find(p_units<0.05);
unresponsive_units = find(p_units>0.05);

%% raw data
% spike counts binned for each stim
% create time vector around all stim onsets
bin_window = 0.001;
bin_edges = -0.5:bin_window:2;
around_stim_time = stimOn_times + bin_edges; 

% for plot
bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;

all_spikes_in_stim_time = nan(length(bc_good_templates_idx), length(trial_stim_values)/length(possible_stim), size(around_stim_time, 2)-1);

for unit_idx=1:length(bc_good_templates_idx)
   
    unit_spikes = spike_templates == unit_idx; 

    unit_spike_times = spike_times_timeline(unit_spikes);

    this_stim_time = around_stim_time(trial_stim_values == contra_stim, :);

    all_spikes_in_stim_time(unit_idx, :, :) = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, this_stim_time(trial_id,:))', ...
        1:size(this_stim_time, 1), 'UniformOutput',false))' / bin_window;
end

mean_all_spikes_in_stim_time = squeeze(mean(all_spikes_in_stim_time, 2));

%% get baseline

% create pre stim onsets
bin_window = 0.2;
baseline_stim_time = stimOn_times - [bin_window 0];

% possible stims
possible_stim = unique(trial_stim_values);
contra_stim = possible_stim(end);

unit_spikes_baseline_stim = nan(length(bc_good_templates_idx), length(trial_stim_values)/length(possible_stim));
for unit_idx=1:length(bc_good_templates_idx)
    
    unit_spikes = spike_templates == unit_idx; 

    unit_spike_times = spike_times_timeline(unit_spikes);

    % spike counts binned pre/post stim
    this_stim_time = baseline_stim_time(trial_stim_values == contra_stim, :);
    unit_spikes_baseline_stim(unit_idx, :)  = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, this_stim_time(trial_id,:))', ...
        1:size(this_stim_time, 1), 'UniformOutput',false))' / bin_window;

end

baseline_all = mean(unit_spikes_baseline_stim, 2);

%% normalize

norm_mean_all_spikes_in_stim_time = (mean_all_spikes_in_stim_time - baseline_all) ...
    ./ (baseline_all + std(baseline_all));

norm_mean_responsive_spikes_in_stim_time = norm_mean_all_spikes_in_stim_time(responsive_units, :);

norm_mean_unresponsive_spikes_in_stim_time = norm_mean_all_spikes_in_stim_time(unresponsive_units, :);


%% sorting units based on normalized values

% get spikes 100ms post stim
post_stim_time = [0 0.1];
post_stim_spikes = mean(norm_mean_all_spikes_in_stim_time(:,bin_centres>post_stim_time(1) & bin_centres<post_stim_time(2)), 2);

% sort for plotting
[sorted_post_stim_spikes, sorted_units] = sort(post_stim_spikes, 'descend');
sorted_mean_all_spikes_in_stim_time = mean_all_spikes_in_stim_time(sorted_units, :);
sorted_norm_mean_all_spikes_in_stim_time = norm_mean_all_spikes_in_stim_time(sorted_units, :);

% get responsive sorted units for plotting
[~, sorted_responsive_units_idx] = sort(post_stim_spikes(responsive_units), 'descend');
sorted_responsive_units = responsive_units(sorted_responsive_units_idx);
sorted_mean_responsive_spikes_in_stim_time = mean_all_spikes_in_stim_time(sorted_responsive_units, :);
sorted_norm_mean_responsive_spikes_in_stim_time = norm_mean_responsive_spikes_in_stim_time(sorted_responsive_units_idx, :);

% get unresponsive sorted units for plotting
[~, sorted_unresponsive_units_idx] = sort(post_stim_spikes(unresponsive_units), 'descend');
sorted_unresponsive_units = unresponsive_units(sorted_unresponsive_units_idx);
sorted_mean_unresponsive_spikes_in_stim_time = mean_all_spikes_in_stim_time(sorted_unresponsive_units, :);
sorted_norm_mean_unresponsive_spikes_in_stim_time = norm_mean_unresponsive_spikes_in_stim_time(sorted_unresponsive_units_idx, :);

%% get mean trace for resp vs unresp

% define and normalize gaussian window
gauss_win = gausswin(51, 3)';

% all units norm
smooth_norm_all_units = filter(gauss_win,sum(gauss_win),norm_mean_all_spikes_in_stim_time, [], 2);
mean_smooth_norm_all_units = mean(smooth_norm_all_units, 1);

% responsive norm
smooth_norm_responsive_units = filter(gauss_win,sum(gauss_win),norm_mean_responsive_spikes_in_stim_time, [], 2);
mean_smooth_norm_responsive_units = mean(smooth_norm_responsive_units, 1);

% unresponsive norm
smooth_norm_unresponsive_units = filter(gauss_win,sum(gauss_win),norm_mean_unresponsive_spikes_in_stim_time, [], 2);
mean_smooth_norm_unresponsive_units = mean(smooth_norm_unresponsive_units, 1);

% figure;
% tiledlayout(2,2)
% nexttile
% imagesc(smooth_norm_responsive_units(sorted_responsive_units_idx, :));
% title('Responsive')
% nexttile
% plot(mean_smooth_norm_responsive_units);
% title('Responsive')
% nexttile
% imagesc(smooth_norm_unresponsive_units(sorted_unresponsive_units_idx, :));
% title('Unresponsive')
% nexttile
% plot(mean_smooth_norm_unresponsive_units);
% title('Unresponsive')


%% plots

%% - raw data
% all cells
% subplot(1,3,1)
% imagesc(bin_centres, [], sorted_mean_all_spikes_in_stim_time);                      
% % colormap(brewermap([], 'RdBu'));
% % caxis([-4500 4500])
% colorbar;
% title('All cells')

%% - normalized data
figure;
sgtitle('Normalized data') 
upper_caxis = max(max(abs(smooth_norm_responsive_units), [],"all"), ...
    max(abs(smooth_norm_unresponsive_units), [],"all"));

% responsive cells
resp = subplot(1,3,1)
imagesc(bin_centres, [], smooth_norm_responsive_units(sorted_responsive_units_idx, :));
xline(0, 'LineWidth', 1);
xline(0.5, 'LineWidth', 1);
colormap(resp, AP_colormap('BWR', [], 0.7));
caxis([-upper_caxis upper_caxis]);
colorbar;
title('Responsive cells');

% unresponsive cells
unresp = subplot(1,3,2);
imagesc(bin_centres, [], smooth_norm_unresponsive_units(sorted_unresponsive_units_idx, :));
xline(0, 'LineWidth', 1);
xline(0.5, 'LineWidth', 1);
colormap(unresp, AP_colormap('BWR', [], 0.7));
caxis([-upper_caxis upper_caxis]);
colorbar;
title('Unresponsive cells');

% mean traces
subplot(1,3,3);
plot(bin_centres, mean_smooth_norm_responsive_units);
hold on;
plot(bin_centres, mean_smooth_norm_unresponsive_units);
hold on;
xline(0, 'LineWidth', 1);
xline(0.5, 'LineWidth', 1);
legend({'Responsive', 'Unresponsive'});


% % all cells
% figure;
% imagesc(bin_centres, [], sorted_norm_mean_all_spikes_in_stim_time);
% % colormap(brewermap([], 'RdBu'));
% % caxis([-4500 4500])
% colorbar;
% title('Norm All cells')

% % zscore
% zscore_sorted_mean_all_spikes_in_stim_time = zscore(sorted_mean_all_spikes_in_stim_time, [], 2);
% figure;
% imagesc(bin_centres, [], zscore_sorted_mean_all_spikes_in_stim_time);
% % colormap(brewermap([], 'RdBu'));
% % caxis([-4500 4500])
% colorbar;
% title('Zscore All cells')


% % check individual cells
% unit_idx = 1;
% figure;
% yyaxis left
% plot(sorted_mean_all_spikes_in_stim_time(unit_idx,:))
% ylim([-1 max(sorted_mean_all_spikes_in_stim_time(unit_idx,:))])
% hold on;
% yyaxis right
% plot(norm_sorted_mean_all_spikes_in_stim_time(unit_idx,:))
% title(['Unit ' num2str(unit_idx)])
% legend({'Raw', 'Norm'})
