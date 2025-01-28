% striatum boundary fix

animal = 'AM014';
rec_day = '2024-01-23';
workflow = 'passive';

% analysis.rerun_bombcell(animal, rec_day)

load_parts.ephys = true;

ap.load_recording

% ap.run_bombcell('P:\Data\AM014\2024-01-23\ephys\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA-AP\continuous.dat', ...
%     'P:\Data\AM014\2024-01-23\ephys\kilosort4', ...
%     'P:\Data\AM014\2024-01-23\ephys\experiment1\recording1\structure.oebin', ...
%     'kilosort4')

trial_stim_values = vertcat(trial_events.values.TrialStimX);
trial_stim_values = trial_stim_values(1:length(stimOn_times));

ap.cellraster(stimOn_times,trial_stim_values);

template_qc_labels(752)


AP_longstriatum_find_striatum_depth


%% depth plot
contra_stim_fig = figure('Position', [970   300   520   550]);
tiledlayout('flow');
str_depth = [str_start str_end];

% plot
unit_axes = nexttile;
set(unit_axes,'YDir','reverse');
hold on;

norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
unit_dots = scatter3(norm_spike_n,template_depths(unique(spike_templates)), ...
    unique(spike_templates),20,'k','filled');

curr_spike_templates = spike_templates(depth_group == curr_depth);
% curr_depth_unit_dots = scatter3(norm_spike_n(curr_spike_templates),template_depths(curr_spike_templates), ...
%     curr_spike_templates,20,'filled', 'MarkerFaceColor', [0.7 0 0]);

x_limits = xlim(unit_axes); % get x-axis limits
fill_x = [x_limits(1), x_limits(2), x_limits(2), x_limits(1)];
fill_y_all = [str_depth(1), str_depth(1), str_depth(2), str_depth(2)];
fill(fill_x, fill_y_all, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Create shaded box

xlim(unit_axes,[-0.1,1]);
ylim([-50, max(template_depths)+50]);
ylabel('Depth (\mum)', 'FontSize', 38)
xlabel('Normalized log rate', 'FontSize', 38)

ax = get(gca);
% Customize tick labels and spacing
ax.XAxis.FontSize = 30;  % Set X-axis tick label font size
ax.YAxis.FontSize = 30;  % Set Y-axis tick label font size

% Set specific intervals for ticks to make them more sparse
xticks([0, 0.5, 1]); % 5 ticks on X-axis
yticks([0, 1725, 3500]); % Customize Y-axis tick intervals