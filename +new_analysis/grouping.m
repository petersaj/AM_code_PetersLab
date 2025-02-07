%% grouping

%% load
save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Users\Andrada-Maria_Marica\long_str_ctx_data';
bhv_data_path = fullfile(save_path, "swr_bhv.mat");
wf_data_path = fullfile(save_path, "ctx_wf.mat");
ephys_data_path = fullfile(save_path, "ephys.mat");
ctx_str_maps_data_path = fullfile(save_path, 'ctx_maps_to_str.mat');

load(bhv_data_path)
load(wf_data_path)
load(ephys_data_path)
load(ctx_str_maps_data_path)

% load("D:\matlab_save\ctx_maps_to_str.mat")

%% replace [] with nans 
% for colIdx = 1:width(ephys)
%     colName = ephys.Properties.VariableNames{colIdx};
%     if iscell(ephys.(colName))
%         ephys.(colName)(cellfun(@isempty, ephys.(colName))) = {NaN};
%     end
% end
% 
% for colIdx = 1:width(all_ctx_maps_to_str)
%     colName = all_ctx_maps_to_str.Properties.VariableNames{colIdx};
%     if iscell(all_ctx_maps_to_str.(colName))
%         all_ctx_maps_to_str.(colName)(cellfun(@isempty, all_ctx_maps_to_str.(colName))) = {NaN};
%     end
% end

%% group
%%%%%% copy of main_ctx_str_maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get number of things from total MUA - after re-run save script
all_cortex_kernel_px = cat(3, all_ctx_maps_to_str.cortex_kernel_px{:});
all_flattened_cortex_kernel_px = reshape(all_cortex_kernel_px, [], ...
    size(all_cortex_kernel_px, 3))';
% 
% use_cortex_kernels = ~any(isnan(all_flattened_cortex_kernel_px), 2);
% nonan_all_flattened_cortex_kernel_px = all_flattened_cortex_kernel_px(use_cortex_kernels);

% reshape(cortex_kernel_px, [], ...
%         size(cortex_kernel_px, 3));
%
% total_depths = 0;
% for animal_idx=1:height(all_ctx_maps_to_str)
%     ctx_maps_to_str = all_ctx_maps_to_str.recording_data{animal_idx};
%     total_depths = total_depths + sum([ctx_maps_to_str.no_depths{:}]);
% end
% % make big cell of all the maps
% all_flattened_cortex_kernel_px = [];
% for rec_idx=1:height(all_ctx_maps_to_str)
%     ctx_maps_to_str = all_ctx_maps_to_str(rec_idx).cortex_kernel;
%     cortex_kernel_px = cat(3, ctx_maps_to_str.cortex_kernel_px{:});
%     flattened_cortex_kernel_px = reshape(cortex_kernel_px, [], ...
%         size(cortex_kernel_px, 3));
%     all_flattened_cortex_kernel_px = cat(1, all_flattened_cortex_kernel_px, flattened_cortex_kernel_px');
% end

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



%% make big vectors of days from learning and mouse id

% n_depths = arrayfun(@(rec_idx) ...
%     length(all_ctx_maps_to_str.depth_group_edges{rec_idx})-1, ...
%     1:height(all_ctx_maps_to_str));
% n_depths(n_depths<0) = 0;

n_depths = arrayfun(@(rec_idx) ...
    size(all_ctx_maps_to_str.cortex_kernel_px{rec_idx}, 3) * ~isempty(all_ctx_maps_to_str.cortex_kernel_px{rec_idx}), ...
    1:height(all_ctx_maps_to_str));

for_psth_days_from_learning = repelem(bhv.days_from_learning, n_depths);
[~, ~, for_psth_animal_ids] = unique(repelem(bhv.animal, n_depths));


%% get psths
num_recordings = height(ephys);
unique_stims_nan = unique(vertcat(ephys.trial_stim_values{:})); 
unique_stims = unique_stims_nan(~isnan(unique_stims_nan));

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


%% group according to animal, cluster and day from learning
%%%%%%% !!!!!!!! fix NaNs !!!!!!!!!! %%%%%%%%%%%%%%%%%%

% only use nonnan learning days
use_days = ~isnan(for_psth_days_from_learning);

% replace nan cluster ids with 0
% nonan_cluster_ids = cluster_ids;
% nonan_cluster_ids(isnan(cluster_ids)) = 0;

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

%% new norm smooth
% for plot
time_vector = -0.5:0.001:1;
baseline_idx = time_vector >= -0.2 & time_vector <= 0;

gauss_win = gausswin(51, 3)';
normalize_and_smooth = @(psth) ...
    filter(gauss_win, sum(gauss_win), ((psth - mean(psth(baseline_idx), 1)) ./ (mean(psth(baseline_idx), 1))), [], 1);

all_smooth_stim_psths = cell(numel(unique_stims), 1);
for stim_idx=1:numel(unique_stims)
    all_smooth_stim_psths{stim_idx} = cell2mat(arrayfun(@(rep) ...
        normalize_and_smooth(squeeze(all_sum_avg_stim_grouped_spikes{stim_idx}(:, rep))), ...
        1:size(all_sum_avg_stim_grouped_spikes{stim_idx}, 2), ...
        'UniformOutput', false));
end

%% average across mice

% mouse idx is first
avg_animal_group_indices = group_indices_unique_clusters(:,2:3);
[unique_avg_animal_group_indices, ~, animal_group_clusters_indices] = unique(avg_animal_group_indices, 'rows');

avg_smooth_stim_psths = cell(numel(unique_stims), 1);
for stim_idx=1:numel(unique_stims)
    avg_smooth_stim_psths{stim_idx} = ap.groupfun(@mean, ...
        all_smooth_stim_psths{stim_idx}, [], animal_group_clusters_indices);
end

% TO DO do numel for number
% num_animals_stim_psths = cell(numel(unique_stims), 1);
% for stim_idx=1:numel(unique_stims)
%     num_animals_stim_psths{stim_idx} = ap.groupfun(@(x) numel(x), ...
%         all_smooth_stim_psths{stim_idx}, [], animal_group_clusters_indices);
% end

num_animals_stim_psths = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    num_animals_stim_psths{stim_idx} = accumarray(animal_group_clusters_indices, 1);
end

% TO DO do sem for errorbar
std_smooth_stim_psths = cell(numel(unique_stims), 1);
sem_smooth_stim_psths = cell(numel(unique_stims), 1);
for stim_idx=1:numel(unique_stims)
    std_smooth_stim_psths{stim_idx} = ap.groupfun(@std, ...
        all_smooth_stim_psths{stim_idx}, [], animal_group_clusters_indices);
    sem_smooth_stim_psths{stim_idx} = std_smooth_stim_psths{stim_idx} ./ sqrt(num_animals_stim_psths{stim_idx})';
end


%% plot all days

stim_idx = 3;

% figure;
% tiledlayout(length(unique_cluster_ids), length(unique_days_from_learning));

days_for_plot = -3:2;
all_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

for cluster_id =1:length(unique_cluster_ids)

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
    %
    %     x_fill = [x, fliplr(x)];  % X values: time points and reversed time points
    %     y_fill = [avg + sem, fliplr(avg - sem)];  % Y values: upper and lower bounds of SEM
    %
    %     % Plot the shaded SEM region
    %     fill(x_fill, y_fill, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');


end

%% sem plot per day

stim_idx = 3;

% figure;
% tiledlayout(length(unique_cluster_ids), length(unique_days_from_learning));

days_for_plot = -3:2;
all_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

for cluster_id =1:length(unique_cluster_ids)

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



%% get max amplitude across days

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


    this_cluster_psth_idx = unique_avg_animal_group_indices(:,1) == cluster_id;
    these_days_from_learning = unique_avg_animal_group_indices(this_cluster_psth_idx, 2);

    plot_day_idx = ismember(these_days_from_learning, days_for_plot);
    for_plot_max_median = median_max_smooth_stim_psths{stim_idx}(this_cluster_psth_idx);
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
    errorbar(plotted_days, for_plot_max_median(:, plot_day_idx), for_plot_max_sem(:, plot_day_idx), '-o', 'CapSize', 0, ...
        'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);

    title(['Cluster ' num2str(cluster_id)])

end

%% frac resp

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


% get unit depth group
all_unit_depth_group = vertcat(ephys.unit_depth_group);
% if size(cat_all_resp_units{stim_idx}) == size(all_unit_depth_group)
%     disp('Sizes match: Resp units and unit depth')
% else
%     disp('Resp units diff from unit depth')
% end

%%%% HERE %%%%%%%%%%%%
% get unit cluster id based on unit depth group

depth_group_edges = vertcat(ephys.depth_group_edges);
depth_indexes = cellfun(@(x) 1:numel(x)-1, depth_group_edges, ...
    'UniformOutput', false,'ErrorHandler', @(x, varargin) []);


all_depth_indexes = horzcat(depth_indexes{:});


per_rec_cluster_ids = mat2cell(cluster_ids, n_depths);


unit_cluster_ids = cellfun(@(cluster, unit_depth) cluster(unit_depth(~isnan(unit_depth))), ...
    per_rec_cluster_ids, all_unit_depth_group, 'uni', false);


nan_idx = @(x,y) nan*~isnan(y)*x(y); 



%%%%%

%%% now some cells are missing - get idx for good cells per recording and
%%% use those

%%% probably better: put nans in unit_cluster_ids instead by doing a for
%%% loop






% make days from learning and animal grouping idx
n_units = arrayfun(@(rec_idx) ...
    size(ephys.unit_depth_group{rec_idx}, 1) * ~isempty(ephys.unit_depth_group{rec_idx}), ...
    1:height(ephys));
for_units_days_from_learning = repelem(bhv.days_from_learning, n_units);
[~, ~, for_units_animal_ids] = unique(repelem(bhv.animal, n_units));




frac_resp_units = cell(numel(unique_stims), 1);

frac_resp_units{stim_idx} = arrayfun(@(rep) sum(all_resp_units{stim_idx}{rep}) / numel(all_resp_units{stim_idx}{rep}), ...
        1:length(all_unit_resp_p_value), 'uni', false, 'ErrorHandler', @(x, varargin) nan);
    all_frac_resp_units{stim_idx} = vertcat(frac_resp_units{stim_idx}{:});

all_frac_resp_units = cell(numel(unique_stims), 1);

% group per mouse and get avg
avg_frac_resp_units = cell(numel(unique_stims), 1);
std_frac_resp_units = cell(numel(unique_stims), 1);
sem_frac_resp_units = cell(numel(unique_stims), 1);
for stim_idx=1:length(unique_stims)
    avg_frac_resp_units = ap.groupfun(@mean, ...
        frac_resp_units{stim_idx}, [], animal_group_clusters_indices);
    std_frac_resp_units{stim_idx} = ap.groupfun(@std, ...
        frac_resp_units{stim_idx}, [], animal_group_clusters_indices);
    sem_frac_resp_units{stim_idx} = std_frac_resp_units{stim_idx} ./ sqrt(num_animals_stim_psths{stim_idx})';
end

%%%%%%%% END HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OLD
%% FOR TRANSFER

%% - psths
% use for all
use_unique_days_from_learning = unique_days_from_learning(4:8);
legend_labels_days = cellfun(@(x) sprintf('Day %d', x), ...
    num2cell(use_unique_days_from_learning), 'UniformOutput', false);

n_days = length(unique_days_from_learning);
contra_good_days = ~any(isnan(contra_max_ampl_grouped_psths), 2);
figure;
tiledlayout(num_clusters,5)
for cluster_id=1:num_clusters
    for day_idx=1:n_days
        if ~ismember(unique_days_from_learning(day_idx), use_unique_days_from_learning)
            continue
        end
        ax = nexttile;
        hold on;
        plot(bin_centres, contra_avg_smooth_norm_grouped_psths{cluster_id, day_idx}, 'k', 'LineWidth', 2);
        xlim([0 0.3])
        ylim([-1 2.5])
        box off;
        ax.XColor = 'none';
        ax.YColor = 'none';
        ax.Color = 'white';
        hold on;
        AP_scalebar(0.1, 0.5)
    end
end
% sgtitle('Contra Stim', 'FontSize', 20, 'FontWeight','bold')

%     nexttile;
%     imagesc(squeeze(centroid_images(cluster_id,:,:)))
%     axis image;
%     clim(max(abs(clim)).*[-1,1]*0.7);
%     ap.wf_draw('ccf','k');
%     colormap(ap.colormap('PWG'));
%     box off;
%     axis off;
%     ylabel(['Map ', num2str(cluster_id)], 'FontSize', 20, 'FontWeight','bold');

%% - max ampl
figure;
contra_days_on_plot = unique_days_from_learning(contra_good_days);
tiledlayout(num_clusters,2);
for cluster_id=1:num_clusters
    nexttile;
    imagesc(squeeze(centroid_images(cluster_id,:,:)))
    axis image;
    clim(max(abs(clim)).*[-1,1]*0.7);
    ap.wf_draw('ccf','k');
    colormap(ap.colormap('PWG'));
    ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight','bold');

    nexttile;
    for mouse_id = 1:length(unique_animal_ids)
        plot(unique_days_from_learning, contra_per_mouse_max_ampl{mouse_id, cluster_id}, ...
            '-o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5], 'Color', [0.5 0.5 0.5]);
        hold on;
    end
    plot(unique_days_from_learning, contra_max_ampl_grouped_psths(:,cluster_id), ...
        '-o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k');    xlabel('Days from association day')
    ylabel('Increase in firing rate')
    xlim([use_unique_days_from_learning(1), use_unique_days_from_learning(end)])
    ylim([-0.5 1])
end
sgtitle('Contra Stim', 'FontSize', 20, 'FontWeight','bold')

% NEW

figure;
contra_days_on_plot = unique_days_from_learning(contra_good_days);
tiledlayout(num_clusters, 2);

for cluster_id = 1:num_clusters
    nexttile;
    imagesc(squeeze(centroid_images(cluster_id, :, :)))
    axis image;
    clim(max(abs(clim)) .* [-1, 1] * 0.7);
    ap.wf_draw('ccf', 'k');
    colormap(ap.colormap('PWG'));
    ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight', 'bold');

    nexttile;
    % Initialize a matrix to hold the data for SEM calculation
    data_for_sem = zeros(length(unique_animal_ids), length(unique_days_from_learning));

    for mouse_id = 1:length(unique_animal_ids)
        % Collect data for each mouse
        data_for_sem(mouse_id, :) = contra_per_mouse_max_ampl{mouse_id, cluster_id};
    end

    % Calculate mean and SEM
    mean_ampl = mean(data_for_sem, 1, 'omitnan');  % Mean along the mouse_id dimension
    sem_ampl = std(data_for_sem, 0, 1, 'omitnan') ./ sqrt(size(data_for_sem, 1));  % SEM

    % Plot mean with error bars
    errorbar(unique_days_from_learning, mean_ampl, sem_ampl, '-o', ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k');
    xlabel('Days from association day');
    ylabel('Increase in firing rate');
    xlim([use_unique_days_from_learning(1), use_unique_days_from_learning(end)]);
    ylim([-0.5 1]);
end

sgtitle('Contra Stim', 'FontSize', 20, 'FontWeight', 'bold');

%% - NEW MAX ALL ON ONE
use_unique_days_from_learning = unique_days_from_learning(3:9);

my_colors = ap.colormap('KR', num_clusters);
new_my_colors = my_colors([3 1 4 2], :);
figure('Position', [680 50 320 950]);
for cluster_id = 1:num_clusters
    % Initialize a matrix to hold the data for SEM calculation
    data_for_sem = zeros(length(unique_animal_ids), length(unique_days_from_learning));

    for mouse_id = 1:length(unique_animal_ids)
        % Collect data for each mouse
        data_for_sem(mouse_id, :) = contra_per_mouse_max_ampl{mouse_id, cluster_id};
    end

    % Calculate mean and SEM
    mean_ampl = mean(data_for_sem, 1, 'omitnan');  % Mean along the mouse_id dimension
    sem_ampl = std(data_for_sem, 0, 1, 'omitnan') ./ sqrt(size(data_for_sem, 1));  % SEM

    median_ampl = median(data_for_sem, 1, 'omitnan');  % Mean along the mouse_id dimension

    % get colour
    curr_color = new_my_colors(cluster_id,:);
    % Plot mean with error bars
    hold on;
    errorbar(unique_days_from_learning, median_ampl, sem_ampl, '-o', 'CapSize', 0, ...
        'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);
end

xlabel('Days from association day', 'FontSize', 20);
ylabel('{\Delta}R/R', 'FontSize', 20);
xlim([use_unique_days_from_learning(1), use_unique_days_from_learning(end)]);

ax = get(gca);
% Customize tick labels and spacing
ax.XAxis.FontSize = 16;  % Set X-axis tick label font size
ax.YAxis.FontSize = 16;  % Set Y-axis tick label font size

% Set specific intervals for ticks to make them more sparse
xticks(-3:2:3); % 5 ticks on X-axis
yticks(0:0.5:2.5); % Customize Y-axis tick intervals


%% - same but different shape for poster
% - NEW MAX ALL ON ONE
use_unique_days_from_learning = unique_days_from_learning(3:9);

my_colors = ap.colormap('KR', num_clusters);
new_my_colors = my_colors([3 1 4 2], :);
figure('Position', [680 50 880 750]);
for cluster_id = 1:num_clusters
    % Initialize a matrix to hold the data for SEM calculation
    data_for_sem = zeros(length(unique_animal_ids), length(unique_days_from_learning));

    for mouse_id = 1:length(unique_animal_ids)
        % Collect data for each mouse
        data_for_sem(mouse_id, :) = contra_per_mouse_max_ampl{mouse_id, cluster_id};
    end

    % Calculate mean and SEM
    mean_ampl = mean(data_for_sem, 1, 'omitnan');  % Mean along the mouse_id dimension
    sem_ampl = std(data_for_sem, 0, 1, 'omitnan') ./ sqrt(size(data_for_sem, 1));  % SEM

    median_ampl = median(data_for_sem, 1, 'omitnan');  % Mean along the mouse_id dimension

    % get colour
    curr_color = new_my_colors(cluster_id,:);
    % Plot mean with error bars
    hold on;
    errorbar(unique_days_from_learning, median_ampl, sem_ampl, '-o', 'CapSize', 0, ...
        'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);
end

xlabel('Days from association day', 'FontSize', 40);
ylabel('{\Delta}R/R', 'FontSize', 40);
xlim([use_unique_days_from_learning(1), use_unique_days_from_learning(end)]);

ax = get(gca);
% Customize tick labels and spacing
ax.XAxis.FontSize = 30;  % Set X-axis tick label font size
ax.YAxis.FontSize = 30;  % Set Y-axis tick label font size

% Set specific intervals for ticks to make them more sparse
xticks(-3:2:3); % 5 ticks on X-axis
yticks(0:0.5:2.5); % Customize Y-axis tick intervals

%%%%%% copy of main_ctx_str_maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% temp check correlation
%
% % Initialize arrays to store the results
% num_maps = size(all_flattened_cortex_kernel_px, 1);
% cosine_similarities = zeros(num_maps, 1);
% correlations = zeros(num_maps, 1);
%
% % Loop through each map and compute similarities with corresponding centroid
% for i = 1:num_maps
%     % Get the corresponding centroid for the assigned cluster
%     assigned_cluster = cluster_ids(i);
%     centroid_vector = centroids(assigned_cluster, :);
%
%     % Get the current map (flattened version)
%     map_vector = all_flattened_cortex_kernel_px(i, :);
%
%     % --- Cosine Similarity ---
%     cosine_similarities(i) = dot(map_vector, centroid_vector) / ...
%                              (norm(map_vector) * norm(centroid_vector));
%
%     % --- Pearson Correlation Coefficient ---
%     R = corrcoef(map_vector, centroid_vector);
%     correlations(i) = R(1, 2);  % Extract the correlation coefficient
% end
%
% % Display some example results
% disp('Cosine Similarities:');
% disp(cosine_similarities(1:10));  % Display first 10 cosine similarities
% disp('Correlations:');
% disp(correlations(1:10));  % Display first 10 correlation coefficients
%
% % Plot histograms of cosine similarities and correlations
% figure;
% subplot(1, 2, 1);
% histogram(cosine_similarities, 20);
% title('Cosine Similarities');
% xlabel('Cosine Similarity');
% ylabel('Frequency');
%
% subplot(1, 2, 2);
% histogram(correlations, 20);
% title('Correlation Coefficients');
% xlabel('Correlation');
% ylabel('Frequency');
%
% % visualize maps around 0.5
%
% % Define a threshold window around 0.5 (e.g., between 0.45 and 0.55)
% threshold_low = 0;
% threshold_high = 0.4;
%
% % Find maps near the threshold for cosine similarity
% cosine_threshold_indices = find(cosine_similarities >= threshold_low & cosine_similarities <= threshold_high);
%
% % Find maps near the threshold for correlation
% correlation_threshold_indices = find(correlations >= threshold_low & correlations <= threshold_high);
%
% % Combine the indices (you can also visualize them separately if needed)
% maps_to_visualize = union(cosine_threshold_indices, correlation_threshold_indices);
%
% % Visualize the maps near the threshold
% figure;
% tiledlayout('flow');
% for i = 90:length(maps_to_visualize)
%     map_idx = maps_to_visualize(i);
%
%     % Reshape the flattened map back to its original size
%     map = reshape(all_flattened_cortex_kernel_px(map_idx, :), [size(cortex_kernel_px, 1), size(cortex_kernel_px, 2)]);
%
%     % Plot the map
%     nexttile;
%     imagesc(map);
%     axis image;
%     colormap(ap.colormap('PWG'));  % Change this to your preferred colormap
%     colorbar;
%     title(['Map Index: ', num2str(map_idx), ', CosSim: ', num2str(cosine_similarities(map_idx)), ', Corr: ', num2str(correlations(map_idx))]);
% end
