%% load
save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Users\Andrada-Maria_Marica\long_str_ctx_data';

bhv_data_path = fullfile(save_path, "swr_bhv.mat");
% wf_data_path = fullfile(save_path, "ctx_wf.mat");
task_ephys_data_path = fullfile(save_path, "task_ephys.mat");
ctx_str_maps_data_path = fullfile(save_path, 'ctx_maps_to_str.mat');

load(bhv_data_path)
% load(wf_data_path)
load(task_ephys_data_path)
load(ctx_str_maps_data_path)
% 
% master_U_fn = fullfile(save_path,'U_master.mat');
% load(master_U_fn, 'U_master');

%% kmeans ids

% get number of things from total MUA - after re-run save script
all_cortex_kernel_px = cat(3, all_ctx_maps_to_str.cortex_kernel_px{:});
all_flattened_cortex_kernel_px = reshape(all_cortex_kernel_px, [], ...
    size(all_cortex_kernel_px, 3))';

% run kmeans
rng('default');
rng(0);
num_clusters = 4;
[cluster_ids, centroids, sumd] = kmeans(double(all_flattened_cortex_kernel_px), num_clusters, 'Distance', 'correlation',  'Replicates',5);


% check who the masters are
centroid_images = reshape(centroids, [num_clusters, [size(all_cortex_kernel_px, 1) size(all_cortex_kernel_px, 2)]]);
figure;
tiledlayout('flow');
for idx=1:num_clusters
    nexttile;
    imagesc(squeeze(centroid_images(idx,:,:)))
    axis image;
    axis off;
    clim(max(abs(clim)).*[-1,1]*0.7);
    ap.wf_draw('ccf','k');
    colormap(ap.colormap('PWG'));
    title(['Cluster ', num2str(idx)]);
end

%% heatmaps

%% - stuff
num_recordings = height(task_ephys);
% unique_stims_nan = unique(vertcat(task_ephys.trial_stim_values{:})); 
% unique_stims = unique_stims_nan(~isnan(unique_stims_nan));

%% - make big vectors of days from learning and mouse id

%% !!!!!!!!!!! LEFT HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

n_depths = arrayfun(@(rec_idx) ...
    size(task_ephys.binned_spikes_stim_align{rec_idx}, 3) * ~isempty(task_ephys.binned_spikes_stim_align{rec_idx}), ...
    1:num_recordings);

maps_n_depths = arrayfun(@(rec_idx) ...
    size(all_ctx_maps_to_str.cortex_kernel_px{rec_idx}, 3) * ~isempty(all_ctx_maps_to_str.cortex_kernel_px{rec_idx}), ...
    1:height(all_ctx_maps_to_str));

find(maps_n_depths - n_depths)

num_trials = arrayfun(@(rec_idx) ...
    size(task_ephys.binned_spikes_stim_align{rec_idx}, 1) * ~isempty(task_ephys.binned_spikes_stim_align{rec_idx}), ...
    1:num_recordings);

per_rec_cluster_ids = mat2cell(cluster_ids, n_depths);


for_psth_cluster_ids_cell = arrayfun(@(rec_idx) repelem(1:num_trials(rec_idx), n_depths(rec_idx))', ...
    1:num_recordings, 'uni', false);
for_psth_cluster_ids_all = vertcat(for_psth_cluster_ids_cell{:});

for_psth_n_depths = repelem(n_depths, num_trials);
size(for_psth_n_depths)
sum(n_depths.*num_trials)

expected_size = sum(num_trials)

% Repeat trial indices for each depth
trial_idx_cell = arrayfun(@(rec_idx) repelem(1:num_trials(rec_idx), n_depths(rec_idx))', ...
    1:num_recordings, 'uni', false);
trial_idx_all = vertcat(trial_idx_cell{:});

% Repeat depth indices for each trial
depth_idx_cell = arrayfun(@(rec_idx) repelem(1:n_depths(rec_idx), num_trials(rec_idx))', ...
    1:num_recordings, 'uni', false);
depth_idx_all = vertcat(depth_idx_cell{:});

% Combine the trial and depth indices into a single matrix
trial_depth_index = [trial_idx_all, depth_idx_all];

% for_psth_days_from_learning = repelem(bhv.days_from_learning, n_depths);
% [~, ~, for_psth_animal_ids] = unique(repelem(bhv.animal, n_depths));


% try to get psths 
test_transp = task_ephys.binned_spikes_stim_align.';
test_non_empty = test_transp(~cellfun(@isempty, test_transp)); % Keeps only non-empty cells

all_stim_binned_spikes = cat(1, test_non_empty{:});


flattened_stim_binned_spikes = arrayfun(@(rec_idx) ...
    reshape(task_ephys.binned_spikes_stim_align{rec_idx}, ...
    n_depths(rec_idx)*num_trials(rec_idx), ...
    size(task_ephys.binned_spikes_stim_align{rec_idx}, 2)), ...
    1:num_recordings, 'uni', false);

cat_flattened_stim_binned_spikes = cat(1, flattened_stim_binned_spikes{:});



A_back = reshape(flattened_stim_binned_spikes{2}, num_trials(2), 1501, n_depths(2));

% figure;
% imagesc(task_ephys.binned_spikes_stim_align{2}(:,:,1))
% 
% figure;
% imagesc(A_back(:,:,1))

% %% PASSIVE STUFF
% %% - get psths
% 
% all_avg_stim_grouped_binned_spikes = cell(numel(unique_stims), 1);
% for stim_idx = 1:numel(unique_stims)
%     stim_grouped_binned_spikes = arrayfun(@(rec_idx) ...
%         task_ephys.binned_spikes_stim_align{rec_idx}(ephys.trial_stim_values{rec_idx} == unique_stims(stim_idx), :, :), ...
%         1:num_recordings, 'UniformOutput', false);
%     avg_stim_grouped_binned_spikes = cellfun(@(rec_spikes) ...
%         squeeze(mean(rec_spikes, 1)), ...
%         stim_grouped_binned_spikes, 'UniformOutput', false);
%     all_avg_stim_grouped_binned_spikes {stim_idx} = horzcat(avg_stim_grouped_binned_spikes{:});
% end

%% - group psths according to animal, cluster and day from learning
% only use nonnan learning days
use_days = ~isnan(for_psth_days_from_learning);

[unique_cluster_ids, ~, ~] = unique(cluster_ids(use_days));
[unique_days_from_learning, ~, ~] = unique(for_psth_days_from_learning(use_days));
[unique_animal_ids, ~, ~] = unique(for_psth_animal_ids(use_days));

% Combine cluster_ids and days into a single matrix for indexing
group_indices = [for_psth_animal_ids(use_days), cluster_ids(use_days), for_psth_days_from_learning(use_days)];
[group_indices_unique_clusters, ~, group_clusters_indices] = unique(group_indices, 'rows');

all_sum_avg_stim_grouped_spikes = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    all_sum_avg_stim_grouped_spikes{stim_idx} = ap.groupfun(@sum, ...
        all_avg_stim_grouped_binned_spikes{stim_idx}(:, use_days), [], group_clusters_indices);
end

% group num trials
grouped_num_trials = cell(numel(unique_stims), 1);
for stim_idx=1:numel(unique_stims)
    grouped_num_trials{stim_idx} = accumarray(group_clusters_indices, for_psth_num_trials{stim_idx}(use_days)', ...
        [length(unique(group_clusters_indices)), 1], @unique);
end


% test = group_indices_unique_clusters(:,[1,3]);
% size(unique(test, 'rows'))
% 
% sum(cellfun(@(x) ~isempty(x), ephys.trial_stim_values) & ~isnan(bhv.days_from_learning))


%% - get baseline
time_vector = -0.5:0.001:1;
baseline_idx = time_vector >= -0.2 & time_vector <= 0;

per_stim_baseline_avg = cellfun(@(x) mean(x(baseline_idx, :), 1), all_sum_avg_stim_grouped_spikes, 'UniformOutput', false);
total_num_trials = cell2mat(grouped_num_trials')';
weighted_sum_stim_baseline = cell2mat(cellfun(@(per_stim_baseline, num_trials) per_stim_baseline .* num_trials', ...
    per_stim_baseline_avg, grouped_num_trials, 'UniformOutput', false));
weighted_avg_stim_baseline = sum(weighted_sum_stim_baseline, 1) ./ sum(total_num_trials, 1);

%% - norm smooth psths
gauss_win = gausswin(51, 3)';
normalize_and_smooth = @(psth, rep) ...
    filter(gauss_win, sum(gauss_win), (psth - weighted_avg_stim_baseline(rep)) ./ (weighted_avg_stim_baseline(rep)), [], 1);

all_smooth_stim_psths = cell(numel(unique_stims), 1);
for stim_idx=1:numel(unique_stims)
    all_smooth_stim_psths{stim_idx} = cell2mat(arrayfun(@(rep) ...
        normalize_and_smooth(squeeze(all_sum_avg_stim_grouped_spikes{stim_idx}(:, rep)), rep), ...
        1:size(all_sum_avg_stim_grouped_spikes{stim_idx}, 2), ...
        'UniformOutput', false));
end


%% - average psths across mice

% mouse idx is first
avg_animal_group_indices = group_indices_unique_clusters(:,2:3);
[unique_avg_animal_group_indices, ~, animal_group_clusters_indices] = unique(avg_animal_group_indices, 'rows');

avg_smooth_stim_psths = cell(numel(unique_stims), 1);
for stim_idx=1:numel(unique_stims)
    avg_smooth_stim_psths{stim_idx} = ap.groupfun(@mean, ...
        all_smooth_stim_psths{stim_idx}, [], animal_group_clusters_indices);
end

num_animals_stim_psths = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    num_animals_stim_psths{stim_idx} = accumarray(animal_group_clusters_indices, 1);
end

% do sem for errorbar
std_smooth_stim_psths = cell(numel(unique_stims), 1);
sem_smooth_stim_psths = cell(numel(unique_stims), 1);
for stim_idx=1:numel(unique_stims)
    std_smooth_stim_psths{stim_idx} = ap.groupfun(@std, ...
        all_smooth_stim_psths{stim_idx}, [], animal_group_clusters_indices);
    sem_smooth_stim_psths{stim_idx} = std_smooth_stim_psths{stim_idx} ./ sqrt(num_animals_stim_psths{stim_idx})';
end
