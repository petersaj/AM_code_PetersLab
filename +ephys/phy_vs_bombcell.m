%% Compare phy and bombcell

%% Load experiment

% for bombcell 
rerun = 1;

animal = 'AP009';

use_workflow = 'lcr_passive';
% use_workflow = {'lcr_passive_fullscreen'};
% use_workflow = {'lcr_passive','lcr_passive_fullscreen'};
% use_workflow = {'stim_wheel_right_stage1','stim_wheel_right_stage2'};
% use_workflow = 'sparse_noise';

recordings = ap.find_recordings(animal,[], use_workflow);

% use_rec = 1;
rec_day = '2023-07-14';
use_rec = strcmp(rec_day,{recordings.day});
% use_rec = length(recordings)-1;

rec_day = recordings(use_rec).day;
rec_time = recordings(use_rec).recording{end};

verbose = true;

ap.load_experiment


trialStim_values = vertcat(trial_events.values.TrialStimX);
% stimOn_times = photodiode_times(1:2:end);

AP_cellraster(stimOn_times,trialStim_values);

%% - Phy section


ephys_quality_control = true;
if ephys_quality_control && exist('cluster_groups','var') % && qMetricsExist == 0
    % If there's a manual classification
    if verbose; disp('Keeping manually labelled good units...'); end

    % Check that all used spike templates have a label
    spike_templates_0idx_unique = unique(spike_templates_0idx);
    if ~all(ismember(spike_templates_0idx_unique,uint32(cluster_groups{1}))) || ...
            ~all(ismember(cluster_groups{2},{'good','mua','noise'}))
        warning([animal ' ' rec_day ': not all templates labeled']);
    end

    % Define good units from labels
    good_templates_idx = uint32(cluster_groups{1}( ...
        strcmp(cluster_groups{2},'good') | strcmp(cluster_groups{2},'mua')));
    good_templates = ismember(0:size(templates,1)-1,good_templates_idx);

    % temp by me
    disp(['Found ' num2str(length(good_templates_idx)) ' good units']);
else
    % If no cluster groups at all, keep all
    warning([animal ' ' rec_day ' - no ephys quality control']);
    disp('No ephys quality control, keeping all and re-indexing');
    good_templates_idx = unique(spike_templates_0idx);
    good_templates = ismember(0:size(templates,1)-1,good_templates_idx);
end

% % only keep good units
% % Throw out all non-good template data
% templates = templates(good_templates,:,:);
% template_depths = template_depths(good_templates);
% waveforms = waveforms(good_templates,:);
% templateDuration = templateDuration(good_templates);
% templateDuration_us = templateDuration_us(good_templates);
% 
% % Throw out all non-good spike data
% good_spike_idx = ismember(spike_templates_0idx,good_templates_idx);
% spike_times_openephys = spike_times_openephys(good_spike_idx);
% spike_templates_0idx = spike_templates_0idx(good_spike_idx);
% template_amplitudes = template_amplitudes(good_spike_idx);
% spike_depths = spike_depths(good_spike_idx);
% spike_times_timeline = spike_times_timeline(good_spike_idx);
% 
% % Rename the spike templates according to the remaining templates
% % (and make 1-indexed from 0-indexed)
% new_spike_idx = nan(max(spike_templates_0idx)+1,1);
% new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
% spike_templates = new_spike_idx(spike_templates_0idx+1);

% % choose my own axons
% curr_waveform = waveforms(find(good_templates_idx==108),:);
%
% % Create the original knot points.
% lengthX = size(waveforms,2);
% x = 1:lengthX;
%
% % Plot it and show how the line has sharp bends.
% plot(x,curr_waveform, '-sr', 'LineWidth', 2);
% % Use splines to interpolate a smoother curve,
% % with 10 times as many points,
% % that goes exactly through the same data points.
% samplingRateIncrease = 10;
% newXSamplePoints = linspace(1, lengthX, lengthX * samplingRateIncrease);
% smoothed_waveforms = spline(x, curr_waveform, newXSamplePoints);
% % Plot smoothedY and show how the line is
% % smooth, and has no sharp bends.
% hold on; % Don't destroy the first curve we plotted.
% plot(newXSamplePoints, smoothed_waveforms, '-ob');
% title('Spline Interpolation Demo', 'FontSize', 20);
% legend('Original Points', 'Spline Points');
%
% % slopes
% slopes = [0, diff(smoothed_waveforms)];
% plot(newXSamplePoints, slopes, 'k-', 'LineWidth', 3);
% % Draw x axis
% line(xlim, [0,0], 'Color', 'k', 'LineWidth', 2);
% grid on;
% legend('Original Points', 'Spline Points', 'Slope');


%% - Bombcell section

% animal = 'AP009';
%
% use_workflow = {'lcr_passive'};
% % use_workflow = {'lcr_passive_fullscreen'};
% % use_workflow = {'lcr_passive','lcr_passive_fullscreen'};
% % use_workflow = {'stim_wheel_right_stage1','stim_wheel_right_stage2'};
% % use_workflow = 'sparse_noise';
%
% recordings = ap.find_recordings(animal,use_workflow);
%
% % use_rec = 1;
% rec_day = '2023-06-30';
% use_rec = strcmp(rec_day,{recordings.day});
% % use_rec = length(recordings)-1;
%
% rec_day = recordings(use_rec).day;
% rec_time = recordings(use_rec).protocol{end};
%
% verbose = true;
%
% ap.load_experiment

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


dif_good_templates = find(bc_good_templates_idx_mask' ~= good_templates);

same_good_templates = find(bc_good_templates_idx_mask' & good_templates);

same_noise_templates = find(bc_noise_idx_mask' & ~good_templates);

% without axons
bc_good_not_phy = find(bc_good_templates_idx_mask' & (bc_good_templates_idx_mask' ~= good_templates));
phy_good_not_bc = find(good_templates & (bc_good_templates_idx_mask' ~= good_templates));


% with axons
bc_good_not_phy = find(bc_inc_axons_idx_mask' & (bc_inc_axons_idx_mask' ~= good_templates));
phy_good_not_bc = find(good_templates & (bc_inc_axons_idx_mask' ~= good_templates));

% take out <100 spikes
new_bc_good_not_phy = bc_good_not_phy(qMetric.nSpikes(bc_good_not_phy)>100);

% new_bc_noise_idx = bc_noise_idx(qMetric.nSpikes(bc_noise_idx))

% plots for extra metrics
figure;
histogram(extra_metric_values.all_units_corr_max_chan,50)
hold on;
xline(extra_metric_values.cutoff_all_units_corr_max_chan, 'r');
title('Corr max channels')


figure;
histogram(extra_metric_values.all_unit_non_zero_corr,50)
hold on;
xline(extra_metric_values.cutoff_all_unit_non_zero_corr, 'r');
title('Corr non-zero channels')


figure;
histogram(extra_metric_values.all_units_amp_decay,50)
hold on;
xline(extra_metric_values.cutoff_all_units_amp_decay, 'r');
title('Amp decay')





% not_good_spatialDecay = qMetric.spatialDecaySlope(bc_good_not_phy);
% new_noise = bc_good_not_phy(not_good_spatialDecay>-0.1); % 15 units, all noise
% extra_cell = bc_good_not_phy(not_good_spatialDecay<=-0.1);% 37 units, still some noise
% qMetric.spatialDecaySlope(extra_cell)'

%
% qMetric.clusterID(qMetric.spatialDecaySlope<=-0.1)

% enough_spikes_units = bc_good_not_phy(qMetric.nSpikes(bc_good_not_phy)>300);




% axons and good phy
% axons_bc_phy = bc_axons & good_templates;

% % bc not phy
% % only keep good units
% % Throw out all non-good template data
% templates = templates(good_templates,:,:);
% template_depths = template_depths(good_templates);
% waveforms = waveforms(good_templates,:);
% templateDuration = templateDuration(good_templates);
% templateDuration_us = templateDuration_us(good_templates);
%
% % Throw out all non-good spike data
% good_spike_idx = ismember(spike_templates_0idx,good_templates_idx);
% spike_times_openephys = spike_times_openephys(good_spike_idx);
% spike_templates_0idx = spike_templates_0idx(good_spike_idx);
% template_amplitudes = template_amplitudes(good_spike_idx);
% spike_depths = spike_depths(good_spike_idx);
% spike_times_timeline = spike_times_timeline(good_spike_idx);
%
% % Rename the spike templates according to the remaining templates
% % (and make 1-indexed from 0-indexed)
% new_spike_idx = nan(max(spike_templates_0idx)+1,1);
% new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
% spike_templates = new_spike_idx(spike_templates_0idx+1);


%% -- keep only axons and plot (run bombcell section first)
% bc_axons_idx = find(bc_axons);

% Throw out all non-good template data
templates = templates(bc_axons_idx,:,:);
template_depths = template_depths(bc_axons_idx);
waveforms = waveforms(bc_axons_idx,:);
templateDuration = templateDuration(bc_axons_idx);
templateDuration_us = templateDuration_us(bc_axons_idx);

% Throw out all non-good spike data
good_spike_idx = ismember(spike_templates_0idx,bc_axons_identity);
spike_times_openephys = spike_times_openephys(good_spike_idx);
spike_templates_0idx = spike_templates_0idx(good_spike_idx);
template_amplitudes = template_amplitudes(good_spike_idx);
spike_depths = spike_depths(good_spike_idx);
spike_times_timeline = spike_times_timeline(good_spike_idx);

% Rename the spike templates according to the remaining templates
% (and make 1-indexed from 0-indexed)
new_spike_idx = nan(max(spike_templates_0idx)+1,1);
new_spike_idx(bc_axons_idx) = 1:length(bc_axons_idx);
spike_templates = new_spike_idx(spike_templates_0idx+1);

trialStim_values = vertcat(trial_events.values.TrialStimX);
% stimOn_times = photodiode_times(1:2:end);

AP_cellraster(stimOn_times,trialStim_values);

figure;
for i=1:size(waveforms,1)
    subplot(11,6,i)
    plot(1:size(waveforms,2),waveforms(i,:))
    title(num2str(i))
end


% questionable_axons = bc_axons_idx(:,[8 9])-1

%% -- keep only neurons and plot (run bombcell section first)
% bc_good_templates_idx = find(bc_good_templates);

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

trialStim_values = vertcat(trial_events.values.TrialStimX);
% stimOn_times = photodiode_times(1:2:end);

AP_cellraster(stimOn_times,trialStim_values);

figure;
for i=1:size(waveforms,1)
    subplot(11,6,i)
    plot(1:size(waveforms,2),waveforms(i,:))
    title(num2str(i))
end

% figure;
% plot(unique(spike_templates), ones(1, size(unique(spike_templates), 1)), '-o');
% hold on;
% plot(bc_good_templates_idx, zeros(1, size(bc_good_templates_idx, 1)), '-o')

%% -- keep only noise and plot (run bombcell section first)
% bc_axons_idx = find(bc_axons);

% Throw out all non-good template data
templates = templates(bc_noise_idx,:,:);
template_depths = template_depths(bc_noise_idx);
waveforms = waveforms(bc_noise_idx,:);
templateDuration = templateDuration(bc_noise_idx);
templateDuration_us = templateDuration_us(bc_noise_idx);

% Throw out all non-good spike data
good_spike_idx = ismember(spike_templates_0idx,bc_noise_identity);
spike_times_openephys = spike_times_openephys(good_spike_idx);
spike_templates_0idx = spike_templates_0idx(good_spike_idx);
template_amplitudes = template_amplitudes(good_spike_idx);
spike_depths = spike_depths(good_spike_idx);
spike_times_timeline = spike_times_timeline(good_spike_idx);

% Rename the spike templates according to the remaining templates
% (and make 1-indexed from 0-indexed)
new_spike_idx = nan(max(spike_templates_0idx)+1,1);
new_spike_idx(bc_noise_idx) = 1:length(bc_noise_idx);
spike_templates = new_spike_idx(spike_templates_0idx+1);

trialStim_values = vertcat(trial_events.values.TrialStimX);
% stimOn_times = photodiode_times(1:2:end);

AP_cellraster(stimOn_times,trialStim_values);

figure;
for i=1:size(waveforms,1)
    subplot(11,6,i)
    plot(1:size(waveforms,2),waveforms(i,:))
    title(num2str(i))
end


