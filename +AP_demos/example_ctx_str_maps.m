
%% Load dataset

animal = 'AM021';
% rec_day = '2024-07-23';
% rec_time = '1340';
workflow = {'lcr_passive'};
recordings = plab.find_recordings(animal, [], workflow);

use_rec = length(recordings)-1;
rec_day = recordings(use_rec).day;
rec_time = recordings(use_rec).recording{end};
verbose = true;
load_parts.behavior = true;
load_parts.ephys = true;
load_parts.widefield = true;
ap.load_recording;


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

% Discretize spikes by depth
depth_group_edges = str_start:mua_length:str_end;
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

%% Regress cortical fluorescence to striatal MUA

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

% Convert cortical kernal V into pixels, plot
cortex_kernel_px = squeeze(plab.wf.svd2px(wf_U(:,:,use_svs),cortex_kernel));
ap.imscroll(cortex_kernel_px);
axis image;
clim(max(abs(clim)).*[-1,1]*0.7);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG'));

%% FOR TRANSFER
%% - plot example map
curr_depth = 1;

figure;
imagesc(cortex_kernel_px(:,:,curr_depth))
axis image;
axis off;
clim(max(abs(clim)).*[-1,1]*0.7);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG'));
title('Example map');
colorbar

%% - plot example mua trace

figure;
plot(time_bin_centers, binned_spikes_std(curr_depth,:))
hold on
plot(time_bin_centers, predicted_spikes(curr_depth,:))
xlabel('Time (s)', 'FontSize', 14)
ylabel('Spikes', 'FontSize', 14)
legend({'MUA trace', 'Predicted MUA trace'})

%% - plot depth units

contra_stim_fig = figure('Position', get(0, 'Screensize'));
tiledlayout('flow');
sgtitle('Units in example map')

str_depth = [str_start str_end];

% plot
unit_axes = nexttile;
set(unit_axes,'YDir','reverse');
hold on;

norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
unit_dots = scatter3(norm_spike_n,template_depths(unique(spike_templates)), ...
    unique(spike_templates),20,'k','filled');

% curr_spike_templates = spike_templates(depth_group == curr_depth);
% curr_depth_unit_dots = scatter3(norm_spike_n(curr_spike_templates),template_depths(curr_spike_templates), ...
%     curr_spike_templates,20,'magenta','filled');

yline(str_depth, 'red')
xlim(unit_axes,[-0.1,1]);
ylim([-50, max(template_depths)+50]);
ylabel('Depth (\mum)')
xlabel('Normalized log rate')

%% - plot mua and stim onsets

% curr depth
downsampled_time = downsample(timelite.timestamps, 10);
end_timebin = downsampled_time(end)+mean(diff(downsampled_time));
neural_downsampled_time = [downsampled_time; end_timebin];
bin_window = mean(diff(downsampled_time));

unique_curr_spike_templates = unique(curr_spike_templates);
str_1_all_trace = cell2mat(arrayfun(@(unit_idx) histcounts(spike_times_timelite(spike_templates == unit_idx), neural_downsampled_time), ...
    unique_curr_spike_templates, 'UniformOutput',false)) / bin_window;
str_1_mua_trace = mean(str_1_all_trace, 1)';

gauss_win = gausswin(51, 3)';
smooth_str_1_mua_trace = filter(gauss_win,sum(gauss_win),str_1_mua_trace, [], 1);

figure; plot(smooth_str_1_mua_trace)

figure;
plot(downsampled_time, smooth_str_1_mua_trace)
hold on;
xline(stimOn_times)

% all
str_all_trace = histcounts(spike_times_timelite, neural_downsampled_time) / bin_window;
str_mua_trace = mean(str_all_trace, 1)';

smooth_str_mua_trace = filter(gauss_win,sum(gauss_win),str_mua_trace, [], 1);

figure; plot(smooth_str_mua_trace)

figure;
plot(downsampled_time, smooth_str_mua_trace, 'k')
hold on;
xline(stimOn_times(contra_good_trials), 'r')
xlim([20 120])

%% - plot avg across trials

contra_stim = 90;
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

psth_opts.window = [-0.5,1];
psth_opts.bin_size = 0.001;

binned_spikes_stim_align = ap.ephys_psth(spike_times_timelite(~isnan(depth_group)), ...
    num2cell(stimOn_times),  ...
    depth_group(~isnan(depth_group)));

% around stim time
bin_edges = psth_opts.window(1):psth_opts.bin_size:psth_opts.window(2);
bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;

% contra stim
contra_no_move_binned_spikes_stim_align = binned_spikes_stim_align(:,:, contra_good_trials);
contra_no_move_psth_stim_align = squeeze(nanmean(contra_no_move_binned_spikes_stim_align, 3));
% - normalize
contra_stim_baseline = mean(contra_no_move_psth_stim_align(:,bin_centres>-0.2 & bin_centres<0), 2);
contra_norm_no_move_psth_stim_align = (contra_no_move_psth_stim_align - contra_stim_baseline) ...
    ./ (contra_stim_baseline + std(contra_stim_baseline));

contra_smooth_norm_no_move_psth_stim_align = filter(gauss_win,sum(gauss_win),contra_norm_no_move_psth_stim_align, [], 2);

all_depths_contra_smooth_norm_stim_align = mean(contra_smooth_norm_no_move_psth_stim_align, 1);

figure;
plot(bin_edges, all_depths_contra_smooth_norm_stim_align, 'k')
xline(0, 'r', 'LineWidth', 2)
xlabel('Time from stim onset (s)', 'FontSize', 14)
ylabel('Increase in firing rate (spks/s)', 'FontSize', 14)
title('Striatum psth eg day')

%% - wf all trials over vis 

figure; colormap(gray);
imagesc(wf_avg); 
axis image off
clim([0, 20000])
ap.wf_draw('ccf', 'y')
roi_poly = drawpolygon;
vis_roi_mask = createMask(roi_poly);
% figure; imagesc(vis_roi_mask); axis image;
% title('ROI mask');

vis_roi_trace = ap.wf_roi(wf_U,wf_V,wf_avg,[],vis_roi_mask);
figure;plot(wf_t, vis_roi_trace, 'b');
title('Visual ROI fluorescence');
xline(stimOn_times(contra_good_trials), 'r')
ylim([-0.025, 0.025])
xlim([20 120])

%% - wf avg trials
roi_stim_align = interp1(wf_t,vis_roi_trace,time_stimulus)';
avg_contra_roi_stim_align = mean(roi_stim_align(:, contra_good_trials), 2);
figure;
plot(timevec, avg_contra_roi_stim_align, 'b')
xline(0, 'r', 'LineWidth', 2)
xline(0.5, 'r', 'LineWidth', 2)
xlabel('Time from stim onset (s)', 'FontSize', 14)
ylabel('{\Delta}F/F', 'FontSize', 14)
title('Vis cortex ROI eg day')

