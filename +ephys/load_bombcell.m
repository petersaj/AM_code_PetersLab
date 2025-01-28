ap.load_recording


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
        end

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


% plot to check
%     trial_stim_values = vertcat(trial_events.values.TrialStimX);
% 
%      AP_cellraster(stimOn_times,trial_stim_values);