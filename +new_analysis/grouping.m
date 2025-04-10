%% grouping

%% load
save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Users\Andrada-Maria_Marica\long_str_ctx_data';
bhv_data_path = fullfile(save_path, "swr_bhv.mat");
wf_data_path = fullfile(save_path, "ctx_wf.mat");
ephys_data_path = fullfile(save_path, "ephys.mat");
ctx_str_maps_data_path = fullfile(save_path, 'ctx_maps_to_str.mat');

load(bhv_data_path)
% load(wf_data_path)
load(ephys_data_path)
load(ctx_str_maps_data_path)

master_U_fn = fullfile(save_path,'U_master.mat');
load(master_U_fn, 'U_master');
%% kmeans ids

% get number of things from total MUA - after re-run save script
all_cortex_kernel_px = cat(3, all_ctx_maps_to_str.cortex_kernel_px{:});

% run kmeans
rng('default');
rng(0);
num_clusters = 4;
kmeans_starting = mean(cell2mat(permute(cellfun(@(x) x(:,:,round(linspace(1,size(x,3),4))), ...
    all_ctx_maps_to_str.cortex_kernel_px(~cellfun(@isempty, ...
    all_ctx_maps_to_str.cortex_kernel_px)),'uni',false),[2,3,4,1])),4);
[cluster_ids, centroids, sumd] = kmeans(...
    reshape(all_cortex_kernel_px,prod(size(U_master,[1,2])),[])',num_clusters, ...
    'Distance','Correlation','start',reshape(kmeans_starting,[],num_clusters)');

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

%% stuff
num_recordings = height(ephys);
unique_stims_nan = unique(vertcat(ephys.trial_stim_values{:})); 
unique_stims = unique_stims_nan(~isnan(unique_stims_nan));


%% PSTH
%% - make big vectors of days from learning and mouse id

n_depths = arrayfun(@(rec_idx) ...
    size(all_ctx_maps_to_str.cortex_kernel_px{rec_idx}, 3) * ~isempty(all_ctx_maps_to_str.cortex_kernel_px{rec_idx}), ...
    1:height(all_ctx_maps_to_str));

for_psth_days_from_learning = repelem(bhv.days_from_learning, n_depths);
[~, ~, for_psth_animal_ids] = unique(repelem(bhv.animal, n_depths));

%% - get num trials
num_trials = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    num_trials{stim_idx} = arrayfun(@(rec_idx) ...
        sum(ephys.trial_stim_values{rec_idx} == unique_stims(stim_idx)), ...
        1:num_recordings);
end
for_psth_num_trials = cellfun(@(x) repelem(x, n_depths), num_trials, 'UniformOutput', false);

%% - get psths

all_avg_stim_grouped_binned_spikes = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    stim_grouped_binned_spikes = arrayfun(@(rec_idx) ...
        ephys.binned_spikes_stim_align{rec_idx}(ephys.trial_stim_values{rec_idx} == unique_stims(stim_idx), :, :), ...
        1:num_recordings, 'UniformOutput', false);
    avg_stim_grouped_binned_spikes = cellfun(@(rec_spikes) ...
        squeeze(mean(rec_spikes, 1)), ...
        stim_grouped_binned_spikes, 'UniformOutput', false);
    all_avg_stim_grouped_binned_spikes {stim_idx} = horzcat(avg_stim_grouped_binned_spikes{:});
end

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
softnorm = 10;
normalize_and_smooth = @(psth, rep) ...
    filter(gauss_win, sum(gauss_win), (psth - weighted_avg_stim_baseline(rep)) ./ (weighted_avg_stim_baseline(rep)+softnorm), [], 1);

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


%% - psths plot all days

stim_idx = 3;

% figure;
% tiledlayout(length(unique_cluster_ids), length(unique_days_from_learning));

days_for_plot = -3:2;
all_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

for cluster_id =1:num_clusters

    this_cluster_psth_idx = unique_avg_animal_group_indices(:,1) == cluster_id;
    these_days_from_learning = unique_avg_animal_group_indices(this_cluster_psth_idx, 2);

    plot_day_idx = ismember(these_days_from_learning, days_for_plot);
    for_plot_psth = avg_smooth_stim_psths{stim_idx}(:, this_cluster_psth_idx);


    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
    my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % time
    x = time_vector;

    % get num animals for legen
    num_animals_all = num_animals_stim_psths{stim_idx}(this_cluster_psth_idx);
    num_animals_plotted = num_animals_all(plot_day_idx);

    % make legend
    legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
        plotted_days, num_animals_plotted, 'UniformOutput', false);

    fig_test = figure;
    plot(x, for_plot_psth(:, plot_day_idx));
    colororder(gca, my_colormap);
    legend(legend_for_plot);
    xline(0, 'LineWidth', 2);

    ylim([-0.5 2.5])

    title(['Cluster' num2str(cluster_id)])
    %
    %     x_fill = [x, fliplr(x)];  % X values: time points and reversed time points
    %     y_fill = [avg + sem, fliplr(avg - sem)];  % Y values: upper and lower bounds of SEM
    %
    %     % Plot the shaded SEM region
    %     fill(x_fill, y_fill, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');


end

%% - psths sem plot per day

stim_idx = 3;

% figure;
% tiledlayout(length(unique_cluster_ids), length(unique_days_from_learning));

days_for_plot = -3:2;
all_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

for cluster_id =1:num_clusters

    this_cluster_psth_idx = unique_avg_animal_group_indices(:,1) == cluster_id;
    these_days_from_learning = unique_avg_animal_group_indices(this_cluster_psth_idx, 2);

    plot_day_idx = ismember(these_days_from_learning, days_for_plot);
    for_plot_psth = avg_smooth_stim_psths{stim_idx}(:, this_cluster_psth_idx);
    avg_all = for_plot_psth(:, plot_day_idx);
    for_plot_sem = sem_smooth_stim_psths{stim_idx}(:, this_cluster_psth_idx);
    sem_all = for_plot_sem(:, plot_day_idx);

    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
    my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % time
    x = time_vector;

    % get num animals for legen
    num_animals_all = num_animals_stim_psths{stim_idx}(this_cluster_psth_idx);
    num_animals_plotted = num_animals_all(plot_day_idx);

    figure;
    for day_idx=1:length(plotted_days)

        nexttile;
        this_day = plotted_days(day_idx);
        title_for_plot = ['Day ' num2str(this_day) ' (n = ' num2str(num_animals_plotted(day_idx)) ')'];

        avg = avg_all(:, day_idx)';
        sem = sem_all(:, day_idx)';

        x_fill = [x, fliplr(x)];  % X values: time points and reversed time points
        y_fill = [avg + sem, fliplr(avg - sem)];  % Y values: upper and lower bounds of SEM

        % Plot the shaded SEM region
        fill(x_fill, y_fill, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        hold on;

        plot(x, avg, 'Color', my_colormap(day_idx, :));
        xline(0, 'LineWidth', 2);
        hold on;
        ylim([-0.5 2.5])

        title(title_for_plot)
    end
    sgtitle(['Cluster ' num2str(cluster_id)])
end



%% - get max amplitude across days

max_window = time_vector>0 & time_vector<0.3;

max_smooth_stim_psths = arrayfun(@(stim_idx) ...
    max(all_smooth_stim_psths{stim_idx}(max_window, :), [], 1), ...
    1:numel(unique_stims), 'UniformOutput', false);

% use avg_animal_group_indices to find across animal avg
mean_max_smooth_stim_psths = cell(numel(unique_stims), 1);
for stim_idx=1:numel(unique_stims)
    mean_max_smooth_stim_psths{stim_idx} = ap.groupfun(@mean, ...
        max_smooth_stim_psths{stim_idx}, [], animal_group_clusters_indices);
end

% median
median_max_smooth_stim_psths = cell(numel(unique_stims), 1);
for stim_idx=1:numel(unique_stims)
    median_max_smooth_stim_psths{stim_idx} = ap.groupfun(@median, ...
        max_smooth_stim_psths{stim_idx}, [], animal_group_clusters_indices);
end

% find sem for max
std_max_smooth_stim_psths = cell(numel(unique_stims), 1);
sem_max_smooth_stim_psths = cell(numel(unique_stims), 1);
for stim_idx=1:numel(unique_stims)
    std_max_smooth_stim_psths{stim_idx} = ap.groupfun(@std, ...
        max_smooth_stim_psths{stim_idx}, [], animal_group_clusters_indices);
    sem_max_smooth_stim_psths{stim_idx} = std_max_smooth_stim_psths{stim_idx} ./ sqrt(num_animals_stim_psths{stim_idx})';
end

%% plot max ampl

stim_idx = 3;

days_for_plot = -3:2;
curr_color = 'r';
% max_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
% colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

figure;
tiledlayout(num_clusters,2);
for cluster_id = 1:num_clusters

    nexttile;
    imagesc(squeeze(centroid_images(cluster_id,:,:)))
    axis image;
    clim(max(abs(clim)).*[-1,1]*0.7);
    ap.wf_draw('ccf','k');
    colormap(ap.colormap('PWG'));
    ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight','bold');


    this_cluster_psth_idx = unique_avg_animal_group_indices(:,1) == cluster_id;
    these_days_from_learning = unique_avg_animal_group_indices(this_cluster_psth_idx, 2);

    plot_day_idx = ismember(these_days_from_learning, days_for_plot);
    for_plot_max_mean = mean_max_smooth_stim_psths{stim_idx}(this_cluster_psth_idx);
    for_plot_max_sem = sem_max_smooth_stim_psths{stim_idx}(this_cluster_psth_idx);

    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
%     my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % get num animals for legen
    num_animals_all = num_animals_stim_psths{stim_idx}(this_cluster_psth_idx);
    num_animals_plotted = num_animals_all(plot_day_idx);

    % make legend
    legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
        plotted_days, num_animals_plotted, 'UniformOutput', false);

    nexttile;
    errorbar(plotted_days, for_plot_max_mean(:, plot_day_idx), for_plot_max_sem(:, plot_day_idx), '-o', 'CapSize', 0, ...
        'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);
    ylim([0 8])
    title(['Cluster ' num2str(cluster_id)])

end

sgtitle('Max psth amplitude')

% ADD line per mouse in grey 


%% FRAC resp

p_thresh = 0.01;
all_unit_resp_p_value = vertcat(ephys.unit_resp_p_value);

% ger resp units
all_resp_units = cell(numel(unique_stims), 1);
cat_all_resp_units = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    all_resp_units{stim_idx} = arrayfun(@(rep) [all_unit_resp_p_value{rep}{stim_idx}{:}] < p_thresh, ...
        1:length(all_unit_resp_p_value), 'uni', false, 'ErrorHandler', @(x, varargin) []);
    cat_all_resp_units{stim_idx} = horzcat(all_resp_units{stim_idx}{:})';
end

% get single units
cat_all_single_units = vertcat(ephys.single_unit_idx{:});

% get unit depth group
all_unit_depth_group = vertcat(ephys.unit_depth_group);

% get unit cluster id based on unit depth group
depth_group_edges = vertcat(ephys.depth_group_edges);
depth_indexes = cellfun(@(x) 1:numel(x)-1, depth_group_edges, ...
    'UniformOutput', false,'ErrorHandler', @(x, varargin) []);
all_depth_indexes = horzcat(depth_indexes{:});

% get per rec cluster ids 
per_rec_cluster_ids = mat2cell(cluster_ids, n_depths);

% get units cluster ids
unit_cluster_ids = cell(length(per_rec_cluster_ids), 1);
for rec_idx=1:length(unit_cluster_ids)
    unit_depth = all_unit_depth_group{rec_idx};
    cluster = per_rec_cluster_ids{rec_idx};

    this_unit_cluster_ids = zeros(size(unit_depth));
    this_unit_cluster_ids(~isnan(unit_depth)) = cluster(unit_depth(~isnan(unit_depth)));
    unit_cluster_ids{rec_idx} = this_unit_cluster_ids;
end
all_unit_cluster_ids = vertcat(unit_cluster_ids{:});

% make days from learning and animal grouping idx
n_units = arrayfun(@(rec_idx) ...
    size(ephys.unit_depth_group{rec_idx}, 1) * ~isempty(ephys.unit_depth_group{rec_idx}), ...
    1:height(ephys));
for_units_days_from_learning = repelem(bhv.days_from_learning, n_units);
[~, ~, for_units_animal_ids] = unique(repelem(bhv.animal, n_units));

% - group per cluster

use_days_units = ~isnan(for_units_days_from_learning);

unit_group_indices = [for_units_animal_ids(use_days_units), all_unit_cluster_ids(use_days_units), for_units_days_from_learning(use_days_units)];
[unit_group_indices_unique_clusters, ~, unit_group_clusters_indices] = unique(unit_group_indices, 'rows');

% get frac resp units
num_resp_units = cell(numel(unique_stims), 1);
num_all_units = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    num_resp_units{stim_idx} = ap.groupfun(@sum, ...
        cat_all_resp_units{stim_idx}(use_days_units),  unit_group_clusters_indices, []);
    num_all_units{stim_idx} = ap.groupfun(@size, ...
        cat_all_resp_units{stim_idx}(use_days_units),  unit_group_clusters_indices, []);
end
frac_resp_units = arrayfun(@(stim_idx) num_resp_units{stim_idx}./num_all_units{stim_idx}, ...
    1:length(unique_stims), 'UniformOutput', false);

% - group per mouse and get avg
unit_avg_animal_unit_group_indices = unit_group_indices_unique_clusters(:,2:3);
[unit_unique_avg_animal_group_indices, ~, unit_animal_group_clusters_indices] = unique(unit_avg_animal_unit_group_indices, 'rows');

num_animals_stim_unit = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    num_animals_stim_unit{stim_idx} = accumarray(unit_animal_group_clusters_indices, 1);
end

avg_frac_resp_units = cell(numel(unique_stims), 1);
std_frac_resp_units = cell(numel(unique_stims), 1);
sem_frac_resp_units = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    avg_frac_resp_units{stim_idx} = ap.groupfun(@mean, ...
        frac_resp_units{stim_idx}, unit_animal_group_clusters_indices, []);
    std_frac_resp_units{stim_idx} = ap.groupfun(@std, ...
        frac_resp_units{stim_idx}, unit_animal_group_clusters_indices, []);
    sem_frac_resp_units{stim_idx} = std_frac_resp_units{stim_idx} ./ sqrt(num_animals_stim_unit{stim_idx});
end

%% - plot frac resp

days_for_plot = -3:2;
curr_color = 'k';
% max_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
% colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

for stim_idx=1:length(unique_stims)

    figure;
    tiledlayout(num_clusters,2);
    for cluster_id = 1:num_clusters

        nexttile;
        imagesc(squeeze(centroid_images(cluster_id,:,:)))
        axis image;
        clim(max(abs(clim)).*[-1,1]*0.7);
        ap.wf_draw('ccf','k');
        colormap(ap.colormap('PWG'));
        ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight','bold');


        this_cluster_unit_idx = unit_unique_avg_animal_group_indices(:,1) == cluster_id;
        these_days_from_learning = unit_unique_avg_animal_group_indices(this_cluster_unit_idx, 2);

        plot_day_idx = ismember(these_days_from_learning, days_for_plot);
        for_plot_avg_frac = avg_frac_resp_units{stim_idx}(this_cluster_unit_idx);
        for_plot_sem_frac = sem_frac_resp_units{stim_idx}(this_cluster_unit_idx);

        % get right colours
        plotted_days = these_days_from_learning(plot_day_idx);
        %     my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

        % get num animals for legen
        num_animals_all = num_animals_stim_unit{stim_idx}(this_cluster_unit_idx);
        num_animals_plotted = num_animals_all(plot_day_idx);

        % make legend
        legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
            plotted_days, num_animals_plotted, 'UniformOutput', false);

        nexttile;
        errorbar(plotted_days, for_plot_avg_frac(plot_day_idx), for_plot_sem_frac(plot_day_idx), '-o', 'CapSize', 0, ...
            'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);

        ylim([0,0.25])
        title(['Cluster ' num2str(cluster_id)])

    end
    frac_resp_title = ['Frac resp units for stim ' num2str(unique_stims(stim_idx))];
    sgtitle(frac_resp_title, 'FontSize', 20);
 end

%% STD and MEAN avg

% all_unit_mean_post_stim = vertcat(ephys.mean_post_stim);
% all_unit_mean_pre_stim = vertcat(ephys.mean_pre_stim);
% all_unit_std_post_stim = vertcat(ephys.std_post_stim);
% 
% unit_mean_post_stim = cell(length(unique_stims), 1);
% unit_mean_pre_stim = cell(length(unique_stims), 1);
% unit_std_post_stim = cell(length(unique_stims), 1);
% cat_unit_mean_post_stim = cell(length(unique_stims), 1);
% cat_unit_mean_pre_stim = cell(length(unique_stims), 1);
% cat_unit_norm_post_stim = cell(length(unique_stims), 1);
% cat_unit_std_post_stim = cell(length(unique_stims), 1);
% for stim_idx=1:length(unique_stims)
%     unit_mean_post_stim{stim_idx} = cellfun(@(x) x{stim_idx}, all_unit_mean_post_stim, 'UniformOutput', false, 'ErrorHandler', @(x, varargin) []);
%     unit_mean_pre_stim{stim_idx} = cellfun(@(x) x{stim_idx}, all_unit_mean_pre_stim, 'UniformOutput', false, 'ErrorHandler', @(x, varargin) []);
%     cat_unit_mean_post_stim{stim_idx} = vertcat(unit_mean_post_stim{stim_idx}{:});
%     cat_unit_mean_pre_stim{stim_idx} = vertcat(unit_mean_pre_stim{stim_idx}{:});
%     cat_unit_norm_post_stim{stim_idx} = (cat_unit_mean_post_stim{stim_idx} - cat_unit_mean_pre_stim{stim_idx}) ...
%         ./ (cat_unit_mean_pre_stim{stim_idx} + 0.1);
%    
%     unit_std_post_stim{stim_idx} = cellfun(@(x) x{stim_idx}, all_unit_std_post_stim, 'UniformOutput', false, 'ErrorHandler', @(x, varargin) []);
%     cat_unit_std_post_stim{stim_idx} = vertcat(unit_std_post_stim{stim_idx}{:});
% end


%%% HERE
%%% AP DEMO
cat_unit_mean_baseline = mean(cell2mat(vertcat(ephys.mean_pre_stim{:})),2);
cat_unit_mean_post_stim = mat2cell(cell2mat(vertcat(ephys.mean_post_stim{:})),sum(n_units),ones(1,3));

softnorm = 0.1;
cat_unit_norm_post_stim = cellfun(@(x) (x-cat_unit_mean_baseline)./ ...
    (cat_unit_mean_baseline+softnorm),cat_unit_mean_post_stim,'uni',false);

cat_unit_std_post_stim = mat2cell(cell2mat(vertcat(ephys.std_post_stim{:})),sum(n_units),ones(1,3));


%%% 
%%%

%% group avg mean
group_unit_mean_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    group_unit_mean_post_stim{stim_idx} = ap.groupfun(@mean, ...
        cat_unit_mean_post_stim{stim_idx}(use_days_units),  unit_group_clusters_indices, []);
end

% avg across animals
avg_unit_mean_post_stim = cell(numel(unique_stims), 1);
std_unit_mean_post_stim = cell(numel(unique_stims), 1);
sem_unit_mean_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    avg_unit_mean_post_stim{stim_idx} = ap.groupfun(@mean, ...
        group_unit_mean_post_stim{stim_idx}, unit_animal_group_clusters_indices, []);
    std_unit_mean_post_stim{stim_idx} = ap.groupfun(@std, ...
        group_unit_mean_post_stim{stim_idx}, unit_animal_group_clusters_indices, []);
    sem_unit_mean_post_stim{stim_idx} = std_unit_mean_post_stim{stim_idx} ./ sqrt(num_animals_stim_unit{stim_idx});
end

%% group avg norm mean
group_unit_norm_mean_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    group_unit_norm_mean_post_stim{stim_idx} = ap.groupfun(@mean, ...
        cat_unit_norm_post_stim{stim_idx}(use_days_units),  unit_group_clusters_indices, []);
end

% avg across animals
avg_unit_norm_mean_post_stim = cell(numel(unique_stims), 1);
std_unit_norm_mean_post_stim = cell(numel(unique_stims), 1);
sem_unit_norm_mean_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    avg_unit_norm_mean_post_stim{stim_idx} = ap.groupfun(@mean, ...
        group_unit_norm_mean_post_stim{stim_idx}, unit_animal_group_clusters_indices, []);
    std_unit_norm_mean_post_stim{stim_idx} = ap.groupfun(@std, ...
        group_unit_norm_mean_post_stim{stim_idx}, unit_animal_group_clusters_indices, []);
    sem_unit_norm_mean_post_stim{stim_idx} = std_unit_norm_mean_post_stim{stim_idx} ./ sqrt(num_animals_stim_unit{stim_idx});
end

%% group avg std
group_unit_std_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    group_unit_std_post_stim{stim_idx} = ap.groupfun(@mean, ...
        cat_unit_std_post_stim{stim_idx}(use_days_units),  unit_group_clusters_indices, []);
end

% avg across animals
avg_unit_std_post_stim = cell(numel(unique_stims), 1);
std_unit_std_post_stim = cell(numel(unique_stims), 1);
sem_unit_std_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    avg_unit_std_post_stim{stim_idx} = ap.groupfun(@mean, ...
        group_unit_std_post_stim{stim_idx}, unit_animal_group_clusters_indices, []);
    std_unit_std_post_stim{stim_idx} = ap.groupfun(@std, ...
        group_unit_std_post_stim{stim_idx}, unit_animal_group_clusters_indices, []);
    sem_unit_std_post_stim{stim_idx} = std_unit_std_post_stim{stim_idx} ./ sqrt(num_animals_stim_unit{stim_idx});
end

%% plot post stim std 

stim_idx = 3;

days_for_plot = -3:2;
curr_color = 'k';
% max_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
% colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

figure;
tiledlayout(num_clusters,2);
for cluster_id = 1:num_clusters

    nexttile;
    imagesc(squeeze(centroid_images(cluster_id,:,:)))
    axis image;
    clim(max(abs(clim)).*[-1,1]*0.7);
    ap.wf_draw('ccf','k');
    colormap(ap.colormap('PWG'));
    ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight','bold');


    this_cluster_unit_idx = unit_unique_avg_animal_group_indices(:,1) == cluster_id;
    these_days_from_learning = unit_unique_avg_animal_group_indices(this_cluster_unit_idx, 2);

    plot_day_idx = ismember(these_days_from_learning, days_for_plot);
    for_plot_avg_std_post_stim = avg_unit_std_post_stim{stim_idx}(this_cluster_unit_idx);
    for_plot_sem_std_post_stim = sem_unit_std_post_stim{stim_idx}(this_cluster_unit_idx);

    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
%     my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % get num animals for legen
    num_animals_all = num_animals_stim_unit{stim_idx}(this_cluster_unit_idx);
    num_animals_plotted = num_animals_all(plot_day_idx);

    % make legend
    legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
        plotted_days, num_animals_plotted, 'UniformOutput', false);

    nexttile;
    errorbar(plotted_days, for_plot_avg_std_post_stim(plot_day_idx), for_plot_sem_std_post_stim(plot_day_idx), '-o', 'CapSize', 0, ...
        'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);

    title(['Cluster ' num2str(cluster_id)])

end

sgtitle('STD post stim', 'FontSize', 20);


%% plot post stim mean

stim_idx = 3;

days_for_plot = -3:2;
curr_color = 'k';
% max_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
% colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

figure;
tiledlayout(num_clusters,2);
for cluster_id = 1:num_clusters

    nexttile;
    imagesc(squeeze(centroid_images(cluster_id,:,:)))
    axis image;
    clim(max(abs(clim)).*[-1,1]*0.7);
    ap.wf_draw('ccf','k');
    colormap(ap.colormap('PWG'));
    ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight','bold');


    this_cluster_unit_idx = unit_unique_avg_animal_group_indices(:,1) == cluster_id;
    these_days_from_learning = unit_unique_avg_animal_group_indices(this_cluster_unit_idx, 2);

    plot_day_idx = ismember(these_days_from_learning, days_for_plot);
    for_plot_avg_mean_post_stim = avg_unit_mean_post_stim{stim_idx}(this_cluster_unit_idx);
    for_plot_sem_mean_post_stim = sem_unit_mean_post_stim{stim_idx}(this_cluster_unit_idx);

    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
%     my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % get num animals for legen
    num_animals_all = num_animals_stim_unit{stim_idx}(this_cluster_unit_idx);
    num_animals_plotted = num_animals_all(plot_day_idx);

    % make legend
    legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
        plotted_days, num_animals_plotted, 'UniformOutput', false);

    nexttile;
    errorbar(plotted_days, for_plot_avg_mean_post_stim(plot_day_idx), for_plot_sem_mean_post_stim(plot_day_idx), '-o', 'CapSize', 0, ...
        'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);

    title(['Cluster ' num2str(cluster_id)])

end

sgtitle('MEAN post stim', 'FontSize', 20);

%% plot post norm stim mean

stim_idx = 3;

days_for_plot = -3:2;
curr_color = 'k';
% max_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
% colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

figure;
tiledlayout(num_clusters,2);
for cluster_id = 1:num_clusters

    nexttile;
    imagesc(squeeze(centroid_images(cluster_id,:,:)))
    axis image;
    clim(max(abs(clim)).*[-1,1]*0.7);
    ap.wf_draw('ccf','k');
    colormap(ap.colormap('PWG'));
    ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight','bold');


    this_cluster_unit_idx = unit_unique_avg_animal_group_indices(:,1) == cluster_id;
    these_days_from_learning = unit_unique_avg_animal_group_indices(this_cluster_unit_idx, 2);

    plot_day_idx = ismember(these_days_from_learning, days_for_plot);
    for_plot_avg_norm_mean_post_stim = avg_unit_norm_mean_post_stim{stim_idx}(this_cluster_unit_idx);
    for_plot_sem_norm_mean_post_stim = sem_unit_norm_mean_post_stim{stim_idx}(this_cluster_unit_idx);

    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
%     my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % get num animals for legen
    num_animals_all = num_animals_stim_unit{stim_idx}(this_cluster_unit_idx);
    num_animals_plotted = num_animals_all(plot_day_idx);

    % make legend
    legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
        plotted_days, num_animals_plotted, 'UniformOutput', false);

    nexttile;
    errorbar(plotted_days, for_plot_avg_norm_mean_post_stim(plot_day_idx), for_plot_sem_norm_mean_post_stim(plot_day_idx), '-o', 'CapSize', 0, ...
        'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);

    title(['Cluster ' num2str(cluster_id)])

end

sgtitle('NORM MEAN post stim', 'FontSize', 20);


%% - for resp units mean and std

% - group per cluster for each stim
resp_unit_group_indices = cell(numel(unique_stims), 1);
resp_unit_group_indices_unique_clusters = cell(numel(unique_stims), 1);
resp_unit_group_clusters_indices = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    resp_unit_group_indices{stim_idx} = unit_group_indices(logical(cat_all_resp_units{stim_idx}(use_days_units)), :);
    [resp_unit_group_indices_unique_clusters{stim_idx}, ~, ...
        resp_unit_group_clusters_indices{stim_idx}] = ...
        unique(resp_unit_group_indices{stim_idx}, 'rows');
end

% - group per mouse and get avg
resp_unit_avg_animal_unit_group_indices = cell(numel(unique_stims), 1);
resp_unit_unique_avg_animal_group_indices = cell(numel(unique_stims), 1);
resp_unit_animal_group_clusters_indices = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    resp_unit_avg_animal_unit_group_indices{stim_idx} = resp_unit_group_indices_unique_clusters{stim_idx}(:,2:3);
    [resp_unit_unique_avg_animal_group_indices{stim_idx}, ~, ...
        resp_unit_animal_group_clusters_indices{stim_idx}] = unique(resp_unit_avg_animal_unit_group_indices{stim_idx}, 'rows');
end

% get num animals
num_animals_stim_resp_unit = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    num_animals_stim_resp_unit{stim_idx} = accumarray(resp_unit_animal_group_clusters_indices{stim_idx}, 1);
end


% group per cluster
group_resp_unit_norm_mean_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    group_resp_unit_norm_mean_post_stim{stim_idx} = ap.groupfun(@mean, ...
        cat_unit_norm_post_stim{stim_idx}(logical(cat_all_resp_units{stim_idx}(use_days_units))), ...
        resp_unit_group_clusters_indices{stim_idx}, []);
end

% avg across animals
avg_resp_unit_norm_mean_post_stim = cell(numel(unique_stims), 1);
std_resp_unit_norm_mean_post_stim = cell(numel(unique_stims), 1);
sem_resp_unit_norm_mean_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    avg_resp_unit_norm_mean_post_stim{stim_idx} = ap.groupfun(@mean, ...
        group_resp_unit_norm_mean_post_stim{stim_idx}, resp_unit_animal_group_clusters_indices{stim_idx}, []);
    std_resp_unit_norm_mean_post_stim{stim_idx} = ap.groupfun(@std, ...
        group_resp_unit_norm_mean_post_stim{stim_idx}, resp_unit_animal_group_clusters_indices{stim_idx}, []);
    sem_resp_unit_norm_mean_post_stim{stim_idx} = std_resp_unit_norm_mean_post_stim{stim_idx} ./ sqrt(num_animals_stim_resp_unit{stim_idx});
end

% std 
group_resp_unit_std_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    group_resp_unit_std_post_stim{stim_idx} = ap.groupfun(@mean, ...
        cat_unit_std_post_stim{stim_idx}(logical(cat_all_resp_units{stim_idx}(use_days_units))), ...
        resp_unit_group_clusters_indices{stim_idx}, []);
end

% avg across animals
avg_resp_unit_std_post_stim = cell(numel(unique_stims), 1);
std_resp_unit_std_post_stim = cell(numel(unique_stims), 1);
sem_resp_unit_std_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    avg_resp_unit_std_post_stim{stim_idx} = ap.groupfun(@mean, ...
        group_resp_unit_std_post_stim{stim_idx}, resp_unit_animal_group_clusters_indices{stim_idx}, []);
    std_resp_unit_std_post_stim{stim_idx} = ap.groupfun(@std, ...
        group_resp_unit_std_post_stim{stim_idx}, resp_unit_animal_group_clusters_indices{stim_idx}, []);
    sem_resp_unit_std_post_stim{stim_idx} = std_resp_unit_std_post_stim{stim_idx} ./ sqrt(num_animals_stim_resp_unit{stim_idx});
end


%% plot resp units post norm stim mean

stim_idx = 3;

days_for_plot = -3:2;
curr_color = 'k';
% max_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
% colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

figure;
tiledlayout(num_clusters,2);
for cluster_id = 1:num_clusters

    nexttile;
    imagesc(squeeze(centroid_images(cluster_id,:,:)))
    axis image;
    clim(max(abs(clim)).*[-1,1]*0.7);
    ap.wf_draw('ccf','k');
    colormap(ap.colormap('PWG'));
    ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight','bold');


    this_cluster_unit_idx = resp_unit_unique_avg_animal_group_indices{stim_idx}(:,1) == cluster_id;
    these_days_from_learning = resp_unit_unique_avg_animal_group_indices{stim_idx}(this_cluster_unit_idx, 2);

    plot_day_idx = ismember(these_days_from_learning, days_for_plot);
    for_plot_avg_resp_norm_mean_post_stim = avg_resp_unit_norm_mean_post_stim{stim_idx}(this_cluster_unit_idx);
    for_plot_sem_resp_norm_mean_post_stim = sem_resp_unit_norm_mean_post_stim{stim_idx}(this_cluster_unit_idx);

    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
%     my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % get num animals for legen
    num_animals_all = num_animals_stim_resp_unit{stim_idx}(this_cluster_unit_idx);
    num_animals_plotted = num_animals_all(plot_day_idx);

    % make legend
    legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
        plotted_days, num_animals_plotted, 'UniformOutput', false);

    nexttile;
    errorbar(plotted_days, for_plot_avg_resp_norm_mean_post_stim(plot_day_idx), for_plot_sem_resp_norm_mean_post_stim(plot_day_idx), '-o', 'CapSize', 0, ...
        'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);

    title(['Cluster ' num2str(cluster_id)])

end

sgtitle('RESP units NORM MEAN post stim', 'FontSize', 20);

%% plot resp units post stim std 

stim_idx = 3;

days_for_plot = -3:2;
curr_color = 'k';
% max_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
% colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

figure;
tiledlayout(num_clusters,2);
for cluster_id = 1:num_clusters

    nexttile;
    imagesc(squeeze(centroid_images(cluster_id,:,:)))
    axis image;
    clim(max(abs(clim)).*[-1,1]*0.7);
    ap.wf_draw('ccf','k');
    colormap(ap.colormap('PWG'));
    ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight','bold');

    this_cluster_unit_idx = resp_unit_unique_avg_animal_group_indices{stim_idx}(:,1) == cluster_id;
    these_days_from_learning = resp_unit_unique_avg_animal_group_indices{stim_idx}(this_cluster_unit_idx, 2);

    plot_day_idx = ismember(these_days_from_learning, days_for_plot);
    for_plot_avg_resp_std_post_stim = avg_resp_unit_std_post_stim{stim_idx}(this_cluster_unit_idx);
    for_plot_sem_resp_std_post_stim = sem_resp_unit_std_post_stim{stim_idx}(this_cluster_unit_idx);

    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
%     my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % get num animals for legen
    num_animals_all = num_animals_stim_resp_unit{stim_idx}(this_cluster_unit_idx);
    num_animals_plotted = num_animals_all(plot_day_idx);

    % make legend
    legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
        plotted_days, num_animals_plotted, 'UniformOutput', false);

    nexttile;
    errorbar(plotted_days, for_plot_avg_resp_std_post_stim(plot_day_idx), for_plot_sem_resp_std_post_stim(plot_day_idx), '-o', 'CapSize', 0, ...
        'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);

    title(['Cluster ' num2str(cluster_id)])

end

sgtitle('RESP units STD post stim', 'FontSize', 20);

%% RESP units PSTHS

% get all psths per stim
all_unit_smooth_event_psths = vertcat(ephys.unit_smooth_event_psths{:});
cat_all_unit_smooth_event_psths = arrayfun(@(stim_idx) cell2mat(all_unit_smooth_event_psths(:, stim_idx)), ...
    1:numel(unique_stims), 'UniformOutput', false);

% resp units
% group per cluster
group_resp_unit_smooth_event_psths = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    group_resp_unit_smooth_event_psths{stim_idx} = ap.groupfun(@mean, ...
        cat_all_unit_smooth_event_psths{stim_idx}(cat_all_resp_units{stim_idx} & use_days_units, :), ...
        resp_unit_group_clusters_indices{stim_idx}, []);
end

% avg across animals
avg_resp_unit_smooth_event_psths = cell(numel(unique_stims), 1);
std_resp_unit_smooth_event_psths = cell(numel(unique_stims), 1);
sem_resp_unit_smooth_event_psths = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    avg_resp_unit_smooth_event_psths{stim_idx} = ap.groupfun(@mean, ...
        group_resp_unit_smooth_event_psths{stim_idx}, resp_unit_animal_group_clusters_indices{stim_idx}, []);
    std_resp_unit_smooth_event_psths{stim_idx} = ap.groupfun(@std, ...
        group_resp_unit_smooth_event_psths{stim_idx}, resp_unit_animal_group_clusters_indices{stim_idx}, []);
    sem_resp_unit_smooth_event_psths{stim_idx} = std_resp_unit_norm_mean_post_stim{stim_idx} ./ sqrt(num_animals_stim_resp_unit{stim_idx});
end

% % non resp
% % group per cluster
% group_resp_unit_smooth_event_psths = cell(numel(unique_stims), 1);
% for stim_idx=1:length(unique_stims)
%     group_resp_unit_smooth_event_psths{stim_idx} = ap.groupfun(@mean, ...
%         cat_all_unit_smooth_event_psths{stim_idx}(logical(cat_all_resp_units{stim_idx}), :), ...
%         resp_unit_group_clusters_indices{stim_idx}, []);
% end
% 
% % avg across animals
% avg_resp_unit_smooth_event_psths = cell(numel(unique_stims), 1);
% std_resp_unit_smooth_event_psths = cell(numel(unique_stims), 1);
% sem_resp_unit_smooth_event_psths = cell(numel(unique_stims), 1);
% for stim_idx=1:length(unique_stims)
%     avg_resp_unit_smooth_event_psths{stim_idx} = ap.groupfun(@mean, ...
%         group_resp_unit_smooth_event_psths{stim_idx}, resp_unit_animal_group_clusters_indices{stim_idx}, []);
%     std_resp_unit_smooth_event_psths{stim_idx} = ap.groupfun(@std, ...
%         group_resp_unit_smooth_event_psths{stim_idx}, resp_unit_animal_group_clusters_indices{stim_idx}, []);
%     sem_resp_unit_smooth_event_psths{stim_idx} = std_resp_unit_norm_mean_post_stim{stim_idx} ./ sqrt(num_animals_stim_resp_unit{stim_idx});
% end

%% - plot RESP psths all days

stim_idx = 3;

% figure;
% tiledlayout(length(unique_cluster_ids), length(unique_days_from_learning));

days_for_plot = -3:2;
all_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

for cluster_id =1:length(unique_cluster_ids)

    this_cluster_unit_idx = resp_unit_unique_avg_animal_group_indices{stim_idx}(:,1) == cluster_id;
    these_days_from_learning = resp_unit_unique_avg_animal_group_indices{stim_idx}(this_cluster_unit_idx, 2);

    plot_day_idx = ismember(these_days_from_learning, days_for_plot);
    for_plot_psth = avg_resp_unit_smooth_event_psths{stim_idx}(this_cluster_unit_idx, :);


    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
    my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % time
    x = time_vector;

    % get num animals for legen
    num_animals_all = num_animals_stim_resp_unit{stim_idx}(this_cluster_unit_idx);
    num_animals_plotted = num_animals_all(plot_day_idx);

    % make legend
    legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
        plotted_days, num_animals_plotted, 'UniformOutput', false);

    fig_test = figure;
    plot(x, for_plot_psth(plot_day_idx, :)');
    colororder(gca, my_colormap);
    legend(legend_for_plot);
    xline(0, 'LineWidth', 2);

    title(['RESP PSTH for Cluster ' num2str(cluster_id)])

%     ylim([-0.5 2.5])
    %
    %     x_fill = [x, fliplr(x)];  % X values: time points and reversed time points
    %     y_fill = [avg + sem, fliplr(avg - sem)];  % Y values: upper and lower bounds of SEM
    %
    %     % Plot the shaded SEM region
    %     fill(x_fill, y_fill, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');


end

%% non resp units
cat_all_non_resp_units = cellfun(@(x) ~x, cat_all_resp_units, 'UniformOutput', false);

% - group per cluster for each stim
non_resp_unit_group_indices = cell(numel(unique_stims), 1);
non_resp_unit_group_indices_unique_clusters = cell(numel(unique_stims), 1);
non_resp_unit_group_clusters_indices = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    non_resp_unit_group_indices{stim_idx} = unit_group_indices(logical(cat_all_non_resp_units{stim_idx}(use_days_units)), :);
    [non_resp_unit_group_indices_unique_clusters{stim_idx}, ~, ...
        non_resp_unit_group_clusters_indices{stim_idx}] = ...
        unique(non_resp_unit_group_indices{stim_idx}, 'rows');
end

% - group per mouse and get avg
non_resp_unit_avg_animal_unit_group_indices = cell(numel(unique_stims), 1);
non_resp_unit_unique_avg_animal_group_indices = cell(numel(unique_stims), 1);
non_resp_unit_animal_group_clusters_indices = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    non_resp_unit_avg_animal_unit_group_indices{stim_idx} = non_resp_unit_group_indices_unique_clusters{stim_idx}(:,2:3);
    [non_resp_unit_unique_avg_animal_group_indices{stim_idx}, ~, ...
        non_resp_unit_animal_group_clusters_indices{stim_idx}] = unique(non_resp_unit_avg_animal_unit_group_indices{stim_idx}, 'rows');
end

% get num animals
num_animals_stim_non_resp_unit = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    num_animals_stim_non_resp_unit{stim_idx} = accumarray(non_resp_unit_animal_group_clusters_indices{stim_idx}, 1);
end


% group per cluster
group_non_resp_unit_norm_mean_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    group_non_resp_unit_norm_mean_post_stim{stim_idx} = ap.groupfun(@mean, ...
        cat_unit_norm_post_stim{stim_idx}(logical(cat_all_non_resp_units{stim_idx}(use_days_units))), ...
        non_resp_unit_group_clusters_indices{stim_idx}, []);
end

% avg across animals
avg_non_resp_unit_norm_mean_post_stim = cell(numel(unique_stims), 1);
std_non_resp_unit_norm_mean_post_stim = cell(numel(unique_stims), 1);
sem_non_resp_unit_norm_mean_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    avg_non_resp_unit_norm_mean_post_stim{stim_idx} = ap.groupfun(@mean, ...
        group_non_resp_unit_norm_mean_post_stim{stim_idx}, non_resp_unit_animal_group_clusters_indices{stim_idx}, []);
    std_non_resp_unit_norm_mean_post_stim{stim_idx} = ap.groupfun(@std, ...
        group_non_resp_unit_norm_mean_post_stim{stim_idx}, non_resp_unit_animal_group_clusters_indices{stim_idx}, []);
    sem_non_resp_unit_norm_mean_post_stim{stim_idx} = std_non_resp_unit_norm_mean_post_stim{stim_idx} ./ sqrt(num_animals_stim_non_resp_unit{stim_idx});
end

% std 
group_non_resp_unit_std_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    group_non_resp_unit_std_post_stim{stim_idx} = ap.groupfun(@mean, ...
        cat_unit_std_post_stim{stim_idx}(logical(cat_all_non_resp_units{stim_idx}(use_days_units))), ...
        non_resp_unit_group_clusters_indices{stim_idx}, []);
end

% avg across animals
avg_non_resp_unit_std_post_stim = cell(numel(unique_stims), 1);
std_non_resp_unit_std_post_stim = cell(numel(unique_stims), 1);
sem_non_resp_unit_std_post_stim = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    avg_non_resp_unit_std_post_stim{stim_idx} = ap.groupfun(@mean, ...
        group_non_resp_unit_std_post_stim{stim_idx}, non_resp_unit_animal_group_clusters_indices{stim_idx}, []);
    std_non_resp_unit_std_post_stim{stim_idx} = ap.groupfun(@std, ...
        group_non_resp_unit_std_post_stim{stim_idx}, non_resp_unit_animal_group_clusters_indices{stim_idx}, []);
    sem_non_resp_unit_std_post_stim{stim_idx} = std_non_resp_unit_std_post_stim{stim_idx} ./ sqrt(num_animals_stim_non_resp_unit{stim_idx});
end

% psths
% group per cluster
group_non_resp_unit_smooth_event_psths = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    group_non_resp_unit_smooth_event_psths{stim_idx} = ap.groupfun(@mean, ...
        cat_all_unit_smooth_event_psths{stim_idx}(cat_all_non_resp_units{stim_idx} & use_days_units, :), ...
        non_resp_unit_group_clusters_indices{stim_idx}, []);
end

% avg across animals
avg_non_resp_unit_smooth_event_psths = cell(numel(unique_stims), 1);
std_non_resp_unit_smooth_event_psths = cell(numel(unique_stims), 1);
sem_non_resp_unit_smooth_event_psths = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    avg_non_resp_unit_smooth_event_psths{stim_idx} = ap.groupfun(@mean, ...
        group_non_resp_unit_smooth_event_psths{stim_idx}, non_resp_unit_animal_group_clusters_indices{stim_idx}, []);
    std_non_resp_unit_smooth_event_psths{stim_idx} = ap.groupfun(@std, ...
        group_non_resp_unit_smooth_event_psths{stim_idx}, non_resp_unit_animal_group_clusters_indices{stim_idx}, []);
    sem_non_resp_unit_smooth_event_psths{stim_idx} = std_non_resp_unit_norm_mean_post_stim{stim_idx} ./ sqrt(num_animals_stim_non_resp_unit{stim_idx});
end

%% - non resp plots
%% - PSTHS
stim_idx = 3;

% figure;
% tiledlayout(length(unique_cluster_ids), length(unique_days_from_learning));

days_for_plot = -3:2;
all_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

for cluster_id =1:length(unique_cluster_ids)

    this_cluster_unit_idx = non_resp_unit_unique_avg_animal_group_indices{stim_idx}(:,1) == cluster_id;
    these_days_from_learning = non_resp_unit_unique_avg_animal_group_indices{stim_idx}(this_cluster_unit_idx, 2);

    plot_day_idx = ismember(these_days_from_learning, days_for_plot);
    for_plot_psth = avg_non_resp_unit_smooth_event_psths{stim_idx}(this_cluster_unit_idx, :);


    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
    my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % time
    x = time_vector;

    % get num animals for legen
    num_animals_all = num_animals_stim_non_resp_unit{stim_idx}(this_cluster_unit_idx);
    num_animals_plotted = num_animals_all(plot_day_idx);

    % make legend
    legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
        plotted_days, num_animals_plotted, 'UniformOutput', false);

    fig_test = figure;
    plot(x, for_plot_psth(plot_day_idx, :)');
    colororder(gca, my_colormap);
    legend(legend_for_plot);
    xline(0, 'LineWidth', 2);

    title(['NON RESP PSTH for Cluster ' num2str(cluster_id)])

%     ylim([-0.5 2.5])
    %
    %     x_fill = [x, fliplr(x)];  % X values: time points and reversed time points
    %     y_fill = [avg + sem, fliplr(avg - sem)];  % Y values: upper and lower bounds of SEM
    %
    %     % Plot the shaded SEM region
    %     fill(x_fill, y_fill, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');


end

%% - plot non resp units post norm stim mean

stim_idx = 3;

days_for_plot = -3:2;
curr_color = 'k';
% max_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
% colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

figure;
tiledlayout(num_clusters,2);
for cluster_id = 1:num_clusters

    nexttile;
    imagesc(squeeze(centroid_images(cluster_id,:,:)))
    axis image;
    clim(max(abs(clim)).*[-1,1]*0.7);
    ap.wf_draw('ccf','k');
    colormap(ap.colormap('PWG'));
    ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight','bold');


    this_cluster_unit_idx = non_resp_unit_unique_avg_animal_group_indices{stim_idx}(:,1) == cluster_id;
    these_days_from_learning = non_resp_unit_unique_avg_animal_group_indices{stim_idx}(this_cluster_unit_idx, 2);

    plot_day_idx = ismember(these_days_from_learning, days_for_plot);
    for_plot_avg_non_resp_norm_mean_post_stim = avg_non_resp_unit_norm_mean_post_stim{stim_idx}(this_cluster_unit_idx);
    for_plot_sem_non_resp_norm_mean_post_stim = sem_non_resp_unit_norm_mean_post_stim{stim_idx}(this_cluster_unit_idx);

    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
%     my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % get num animals for legen
    num_animals_all = num_animals_stim_non_resp_unit{stim_idx}(this_cluster_unit_idx);
    num_animals_plotted = num_animals_all(plot_day_idx);

    % make legend
    legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
        plotted_days, num_animals_plotted, 'UniformOutput', false);

    nexttile;
    errorbar(plotted_days, for_plot_avg_non_resp_norm_mean_post_stim(plot_day_idx), for_plot_sem_non_resp_norm_mean_post_stim(plot_day_idx), '-o', 'CapSize', 0, ...
        'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);

    title(['Cluster ' num2str(cluster_id)])

end

sgtitle('NON RESP units NORM MEAN post stim', 'FontSize', 20);

%% - plot non resp units post stim std 

stim_idx = 3;

days_for_plot = -3:2;
curr_color = 'k';
% max_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
% colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

figure;
tiledlayout(num_clusters,2);
for cluster_id = 1:num_clusters

    nexttile;
    imagesc(squeeze(centroid_images(cluster_id,:,:)))
    axis image;
    clim(max(abs(clim)).*[-1,1]*0.7);
    ap.wf_draw('ccf','k');
    colormap(ap.colormap('PWG'));
    ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight','bold');

    this_cluster_unit_idx = non_resp_unit_unique_avg_animal_group_indices{stim_idx}(:,1) == cluster_id;
    these_days_from_learning = non_resp_unit_unique_avg_animal_group_indices{stim_idx}(this_cluster_unit_idx, 2);

    plot_day_idx = ismember(these_days_from_learning, days_for_plot);
    for_plot_avg_non_resp_std_post_stim = avg_non_resp_unit_std_post_stim{stim_idx}(this_cluster_unit_idx);
    for_plot_sem_non_resp_std_post_stim = sem_non_resp_unit_std_post_stim{stim_idx}(this_cluster_unit_idx);

    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
%     my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % get num animals for legen
    num_animals_all = num_animals_stim_non_resp_unit{stim_idx}(this_cluster_unit_idx);
    num_animals_plotted = num_animals_all(plot_day_idx);

    % make legend
    legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
        plotted_days, num_animals_plotted, 'UniformOutput', false);

    nexttile;
    errorbar(plotted_days, for_plot_avg_non_resp_std_post_stim(plot_day_idx), for_plot_sem_non_resp_std_post_stim(plot_day_idx), '-o', 'CapSize', 0, ...
        'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);

    title(['Cluster ' num2str(cluster_id)])

end

sgtitle('NON RESP units STD post stim', 'FontSize', 20);

%% temp check correlation
% %
% % % Initialize arrays to store the results
% % num_maps = size(all_flattened_cortex_kernel_px, 1);
% % cosine_similarities = zeros(num_maps, 1);
% % correlations = zeros(num_maps, 1);
% %
% % % Loop through each map and compute similarities with corresponding centroid
% % for i = 1:num_maps
% %     % Get the corresponding centroid for the assigned cluster
% %     assigned_cluster = cluster_ids(i);
% %     centroid_vector = centroids(assigned_cluster, :);
% %
% %     % Get the current map (flattened version)
% %     map_vector = all_flattened_cortex_kernel_px(i, :);
% %
% %     % --- Cosine Similarity ---
% %     cosine_similarities(i) = dot(map_vector, centroid_vector) / ...
% %                              (norm(map_vector) * norm(centroid_vector));
% %
% %     % --- Pearson Correlation Coefficient ---
% %     R = corrcoef(map_vector, centroid_vector);
% %     correlations(i) = R(1, 2);  % Extract the correlation coefficient
% % end
% %
% % % Display some example results
% % disp('Cosine Similarities:');
% % disp(cosine_similarities(1:10));  % Display first 10 cosine similarities
% % disp('Correlations:');
% % disp(correlations(1:10));  % Display first 10 correlation coefficients
% %
% % % Plot histograms of cosine similarities and correlations
% % figure;
% % subplot(1, 2, 1);
% % histogram(cosine_similarities, 20);
% % title('Cosine Similarities');
% % xlabel('Cosine Similarity');
% % ylabel('Frequency');
% %
% % subplot(1, 2, 2);
% % histogram(correlations, 20);
% % title('Correlation Coefficients');
% % xlabel('Correlation');
% % ylabel('Frequency');
% %
% % % visualize maps around 0.5
% %
% % % Define a threshold window around 0.5 (e.g., between 0.45 and 0.55)
% % threshold_low = 0;
% % threshold_high = 0.4;
% %
% % % Find maps near the threshold for cosine similarity
% % cosine_threshold_indices = find(cosine_similarities >= threshold_low & cosine_similarities <= threshold_high);
% %
% % % Find maps near the threshold for correlation
% % correlation_threshold_indices = find(correlations >= threshold_low & correlations <= threshold_high);
% %
% % % Combine the indices (you can also visualize them separately if needed)
% % maps_to_visualize = union(cosine_threshold_indices, correlation_threshold_indices);
% %
% % % Visualize the maps near the threshold
% % figure;
% % tiledlayout('flow');
% % for i = 90:length(maps_to_visualize)
% %     map_idx = maps_to_visualize(i);
% %
% %     % Reshape the flattened map back to its original size
% %     map = reshape(all_flattened_cortex_kernel_px(map_idx, :), [size(cortex_kernel_px, 1), size(cortex_kernel_px, 2)]);
% %
% %     % Plot the map
% %     nexttile;
% %     imagesc(map);
% %     axis image;
% %     colormap(ap.colormap('PWG'));  % Change this to your preferred colormap
% %     colorbar;
% %     title(['Map Index: ', num2str(map_idx), ', CosSim: ', num2str(cosine_similarities(map_idx)), ', Corr: ', num2str(correlations(map_idx))]);
% % end
