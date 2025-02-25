%% SAVE
%% - Load dataset
%% load
save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Users\Andrada-Maria_Marica\long_str_ctx_data';
bhv_data_path = fullfile(save_path, "swr_bhv.mat");
wf_data_path = fullfile(save_path, "ctx_wf.mat");

load(bhv_data_path)
load(wf_data_path)

% temp
ephys_data_path = fullfile(save_path, "ephys.mat");
load(ephys_data_path)

%% - make roi for vis
% animal = 'AM021';
% workflow = {'lcr_passive'};
% recordings = plab.find_recordings(animal, [], workflow);
% use_rec = length(recordings)-1;
% rec_day = recordings(use_rec).day;
% rec_time = recordings(use_rec).recording{end};
% verbose = true;
% load_parts.behavior = true;
% load_parts.widefield = true;
% ap.load_recording
% 
% % vis
% figure; colormap(gray);
% imagesc(wf_avg);
% axis image off
% clim([0, 20000])
% ap.wf_draw('ccf', 'y')
% vis_roi_poly = drawpolygon;
% vis_roi_mask = createMask(vis_roi_poly);
% 
% % pfc
% figure; colormap(gray);
% imagesc(wf_avg);
% axis image off
% clim([0, 20000])
% ap.wf_draw('ccf', 'y')
% pfc_roi_poly = drawpolygon;
% pfc_roi_mask = createMask(pfc_roi_poly);
% % save('D:\matlab_save\new_pfc_ROI', "new_pfc_roi_mask", "-v7.3");
% 
% figure; imagesc(vis_roi_mask); axis image;
% title('vis ROI mask');
% 
% figure; imagesc(pfc_roi_mask); axis image;
% title('pfc ROI mask');
% 
% save('D:\matlab_save\ROIs', "vis_roi_mask", "pfc_roi_mask", "-v7.3");

%% load ROIs
load('D:\matlab_save\OLD\ROIs.mat');

load("D:\matlab_save\OLD\new_pfc_ROI.mat");
figure; imagesc(new_pfc_roi_mask); axis image;
title('new pfc ROI mask');

master_U_fn = fullfile(save_path,'U_master.mat');

load(master_U_fn, 'U_master');

%% random
num_recordings = height(wf);
unique_stims_nan = unique(vertcat(wf.trial_stim_values{:})); 
unique_stims = unique_stims_nan(~isnan(unique_stims_nan));

wf_stim_time = wf.wf_stim_time{1};

%% Get avg Vs and baseline subtract

% group Vs per stim
all_avg_stim_Vs = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    stim_grouped_Vs = arrayfun(@(rec_idx) ...
        wf.V_no_move_stim_align{rec_idx}(wf.trial_stim_values{rec_idx} == unique_stims(stim_idx), :, :), ...
        1:num_recordings, 'UniformOutput', false);
    temp_avg_stim_grouped_Vs = cellfun(@(rec_V) ...
        squeeze(mean(rec_V, 1)), ...
        stim_grouped_Vs, 'UniformOutput', false);
    all_avg_stim_Vs {stim_idx} = cat(3, temp_avg_stim_grouped_Vs{:});
end

% get baseline
baseline_idx = wf_stim_time > -0.2 & wf_stim_time < 0;
all_stim_baseline = cellfun(@(V) mean(V(baseline_idx, :, :), 1), all_avg_stim_Vs, 'UniformOutput', false);
avg_stim_baseline = squeeze(mean(cell2mat(all_stim_baseline), 1));

all_norm_avg_stim_Vs = cellfun(@(V) V - permute(repmat(avg_stim_baseline, [1, 1, length(wf_stim_time)]), [3, 1, 2]), all_avg_stim_Vs, 'UniformOutput', false);


%% group Vs
% get grouping by learning day 
for_wf_days_from_learning = bhv.days_from_learning;
[~, ~, for_wf_animal_ids] = unique(bhv.animal);

use_days = ~isnan(for_wf_days_from_learning);
[wf_unique_avg_animal_group_indices, ~, wf_animal_group_clusters_indices] = unique(for_wf_days_from_learning(use_days));

% avg across animals
avg_grouped_norm_stim_Vs = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    avg_grouped_norm_stim_Vs{stim_idx} = ap.groupfun(@nanmean, ...
        all_norm_avg_stim_Vs{stim_idx}(:,:,use_days), [], [], wf_animal_group_clusters_indices);
end

for_plot_Vs = cellfun(@(x) permute(x, [2, 1, 3]), avg_grouped_norm_stim_Vs, 'UniformOutput', false);

%% count
num_animals_stim_wf = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    num_animals_stim_wf{stim_idx} = accumarray(wf_animal_group_clusters_indices, 1);
end

%% get pfc roi
all_norm_stim_pfc_roi = cellfun(@(V) cell2mat(arrayfun(@(rec_idx) ap.wf_roi(U_master,V(:,:, rec_idx)',[],[],new_pfc_roi_mask), ...
    1:num_recordings, 'UniformOutput', false)'), ...
    all_norm_avg_stim_Vs, 'UniformOutput', false);

avg_grouped_norm_stim_pfc_roi = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    avg_grouped_norm_stim_pfc_roi{stim_idx} = ap.groupfun(@nanmean, ...
        all_norm_stim_pfc_roi{stim_idx}(use_days, :)', [], wf_animal_group_clusters_indices);
end

% do sem for errorbar
std_grouped_norm_stim_pfc_roi = cell(numel(unique_stims), 1);
sem_grouped_norm_stim_pfc_roi = cell(numel(unique_stims), 1);
for stim_idx=1:numel(unique_stims)
    std_grouped_norm_stim_pfc_roi{stim_idx} = ap.groupfun(@nanstd, ...
        all_norm_stim_pfc_roi{stim_idx}(use_days, :)', [], wf_animal_group_clusters_indices);
    sem_grouped_norm_stim_pfc_roi{stim_idx} = std_grouped_norm_stim_pfc_roi{stim_idx} ./ sqrt(num_animals_stim_wf{stim_idx})';
end

%% get vis ROI 
all_norm_stim_vis_roi = cellfun(@(V) cell2mat(arrayfun(@(rec_idx) ap.wf_roi(U_master,V(:,:, rec_idx)',[],[],vis_roi_mask), ...
    1:num_recordings, 'UniformOutput', false)'), ...
    all_norm_avg_stim_Vs, 'UniformOutput', false);

avg_grouped_norm_stim_vis_roi = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    avg_grouped_norm_stim_vis_roi{stim_idx} = ap.groupfun(@nanmean, ...
        all_norm_stim_vis_roi{stim_idx}(use_days, :)', [], wf_animal_group_clusters_indices);
end

% do sem for errorbar
std_grouped_norm_stim_vis_roi = cell(numel(unique_stims), 1);
sem_grouped_norm_stim_vis_roi = cell(numel(unique_stims), 1);
for stim_idx=1:numel(unique_stims)
    std_grouped_norm_stim_vis_roi{stim_idx} = ap.groupfun(@nanstd, ...
        all_norm_stim_vis_roi{stim_idx}(use_days, :)', [], wf_animal_group_clusters_indices);
    sem_grouped_norm_stim_vis_roi{stim_idx} = std_grouped_norm_stim_vis_roi{stim_idx} ./ sqrt(num_animals_stim_wf{stim_idx})';
end

%% get avg image
avg_image_time_idx = wf_stim_time > 0 & wf_stim_time < 0.3;
avg_image_time_px = cellfun(@(x) plab.wf.svd2px(U_master, x(:, avg_image_time_idx, :)), for_plot_Vs, 'uni', false);

mean_avg_image_px = cellfun(@(x) squeeze(nanmean(x, 3)), avg_image_time_px, 'UniformOutput', false);
max_avg_image_px = cellfun(@(x) squeeze(max(abs(x), [], 3)), avg_image_time_px, 'UniformOutput', false);

%% get max amplitude pfc ROI in window 
max_ampl_window = wf_stim_time > 0 & wf_stim_time < 0.3;

% get max per rec
all_max_ampl_pfc_roi = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    all_max_ampl_pfc_roi{stim_idx} = max(all_norm_stim_pfc_roi{stim_idx}(:, max_ampl_window)', [], 1);
end

% group and get mean across mice
mean_max_ampl_pfc_roi = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    mean_max_ampl_pfc_roi{stim_idx} = ap.groupfun(@nanmean, ...
        all_max_ampl_pfc_roi{stim_idx}(use_days)', wf_animal_group_clusters_indices);
end

median_max_ampl_pfc_roi = cell(numel(unique_stims), 1);
for stim_idx = 1:numel(unique_stims)
    median_max_ampl_pfc_roi{stim_idx} = ap.groupfun(@median, ...
        all_max_ampl_pfc_roi{stim_idx}(use_days)', wf_animal_group_clusters_indices);
end

% do sem for errorbar
std_max_ampl_pfc_roi = cell(numel(unique_stims), 1);
sem_max_ampl_pfc_roi = cell(numel(unique_stims), 1);
for stim_idx=1:numel(unique_stims)
    std_max_ampl_pfc_roi{stim_idx} = ap.groupfun(@nanstd, ...
        all_max_ampl_pfc_roi{stim_idx}(use_days)', wf_animal_group_clusters_indices);
    sem_max_ampl_pfc_roi{stim_idx} = std_max_ampl_pfc_roi{stim_idx} ./ sqrt(num_animals_stim_wf{stim_idx});
end

%% pfc max ampl plots

days_for_plot = -3:2;
curr_color = 'k';
these_days_from_learning = wf_unique_avg_animal_group_indices;
plot_day_idx = ismember(these_days_from_learning, days_for_plot);

for stim_idx=1:length(unique_stims)
        
    for_plot_mean_max_pfc_roi = mean_max_ampl_pfc_roi{stim_idx};
    for_plot_sem_max_pfc_roi = sem_max_ampl_pfc_roi{stim_idx};

    plotted_days = these_days_from_learning(plot_day_idx);

%     % get num animals for legend
%     num_animals_plotted = num_animals_stim_wf{stim_idx}(plot_day_idx);
% 
%     % make legend
%     legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
%         plotted_days, num_animals_plotted, 'UniformOutput', false);

    figure;
    errorbar(plotted_days, for_plot_mean_max_pfc_roi(plot_day_idx), for_plot_sem_max_pfc_roi(plot_day_idx), '-o', 'CapSize', 0, ...
        'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);
    title(['Mean max ampl pfc ROI for stim ' num2str(unique_stims(stim_idx))])
end

days_for_plot = -3:2;
curr_color = 'k';
these_days_from_learning = wf_unique_avg_animal_group_indices;
plot_day_idx = ismember(these_days_from_learning, days_for_plot);

for stim_idx=1:length(unique_stims)
        
    for_plot_median_max_pfc_roi = median_max_ampl_pfc_roi{stim_idx};
    for_plot_sem_max_pfc_roi = sem_max_ampl_pfc_roi{stim_idx};

    plotted_days = these_days_from_learning(plot_day_idx);

%     % get num animals for legend
%     num_animals_plotted = num_animals_stim_wf{stim_idx}(plot_day_idx);
% 
%     % make legend
%     legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
%         plotted_days, num_animals_plotted, 'UniformOutput', false);

    figure;
    errorbar(plotted_days, for_plot_median_max_pfc_roi(plot_day_idx), for_plot_sem_max_pfc_roi(plot_day_idx), '-o', 'CapSize', 0, ...
        'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'Color', curr_color);
    title(['Median max ampl pfc ROI for stim ' num2str(unique_stims(stim_idx))])
end


%% pfc ROI plots

days_for_plot = -3:2;
all_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));
these_days_from_learning = wf_unique_avg_animal_group_indices;
plot_day_idx = ismember(these_days_from_learning, days_for_plot);

for stim_idx=1:numel(unique_stims)

    for_plot_wf_roi = avg_grouped_norm_stim_pfc_roi{stim_idx}(:, plot_day_idx);

    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
    my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % get num animals for legend
    num_animals_plotted = num_animals_stim_wf{stim_idx}(plot_day_idx);

    % make legend
    legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
        plotted_days, num_animals_plotted, 'UniformOutput', false);

    fig_test = figure;
    plot(wf_stim_time, for_plot_wf_roi);
    colororder(gca, my_colormap);
    legend(legend_for_plot);
    xline(0, 'LineWidth', 2);

%     ylim([-0.5 2.5])

    title(['pfc ROI for stim ' num2str(unique_stims(stim_idx))])
end

%% vis ROI plots

days_for_plot = -3:2;
all_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));
these_days_from_learning = wf_unique_avg_animal_group_indices;
plot_day_idx = ismember(these_days_from_learning, days_for_plot);

for stim_idx=1:numel(unique_stims)

    for_plot_wf_roi = avg_grouped_norm_stim_vis_roi{stim_idx}(:, plot_day_idx);

    % get right colours
    plotted_days = these_days_from_learning(plot_day_idx);
    my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

    % get num animals for legend
    num_animals_plotted = num_animals_stim_wf{stim_idx}(plot_day_idx);

    % make legend
    legend_for_plot = arrayfun(@(day, num) ['Day ' num2str(day) ' (n = ' num2str(num) ')'], ...
        plotted_days, num_animals_plotted, 'UniformOutput', false);

    fig_test = figure;
    plot(wf_stim_time, for_plot_wf_roi);
    colororder(gca, my_colormap);
    legend(legend_for_plot);
    xline(0, 'LineWidth', 2);

%     ylim([-0.5 2.5])

    title(['right vis ROI for stim ' num2str(unique_stims(stim_idx))])
end
%% pfc ROI sem plots per day 

days_for_plot = -3:2;
all_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

these_days_from_learning = wf_unique_avg_animal_group_indices;
plot_day_idx = ismember(these_days_from_learning, days_for_plot);

% get right colours
plotted_days = these_days_from_learning(plot_day_idx);
my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

for stim_idx=1:numel(unique_stims)

    avg_all = avg_grouped_norm_stim_pfc_roi{stim_idx}(:, plot_day_idx);
    sem_all = sem_grouped_norm_stim_pfc_roi{stim_idx}(:, plot_day_idx);

    % time
    x = wf_stim_time;

    % get num animals for legen
    num_animals_plotted = num_animals_stim_wf{stim_idx}(plot_day_idx);

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
        ylim([-5 15]*10^(-4))

        title(title_for_plot)
    end
    sgtitle(['pfc ROI Stim ' num2str(unique_stims(stim_idx))])
end

%% vis ROI sem plots per day 

days_for_plot = -3:2;
all_colormap = ap.colormap('BKR', 2*max(abs(days_for_plot))+1);
colormap_days = -max(abs(days_for_plot)):max(abs(days_for_plot));

these_days_from_learning = wf_unique_avg_animal_group_indices;
plot_day_idx = ismember(these_days_from_learning, days_for_plot);

% get right colours
plotted_days = these_days_from_learning(plot_day_idx);
my_colormap = all_colormap(ismember(colormap_days, plotted_days), :);

for stim_idx=1:numel(unique_stims)

    avg_all = avg_grouped_norm_stim_vis_roi{stim_idx}(:, plot_day_idx);
    sem_all = sem_grouped_norm_stim_vis_roi{stim_idx}(:, plot_day_idx);

    % time
    x = wf_stim_time;

    % get num animals for legen
    num_animals_plotted = num_animals_stim_wf{stim_idx}(plot_day_idx);

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
        ylim([-5*10^(-4) 3.5*10^(-3)])

        title(title_for_plot)
    end
    sgtitle(['right vis ROI Stim ' num2str(unique_stims(stim_idx))])
end

%% movies

days_for_plot = -3:2;
these_days_from_learning = wf_unique_avg_animal_group_indices;
plot_day_idx = ismember(these_days_from_learning, days_for_plot);

for stim_idx=1:length(unique_stims)
    movie = plab.wf.svd2px(U_master, for_plot_Vs{stim_idx}(:, :, plot_day_idx));
    ap.imscroll(movie, wf_stim_time)
end

% test_idx = find(plot_day_idx);
% weird_day_idx = test_idx(1);
% ap.imscroll(avg_image_time_px{3}(:,:,:,weird_day_idx), wf_stim_time(avg_image_time_idx))
% 
% imagesc(nanmean(avg_image_time_px{3}(:,:,:,weird_day_idx), 3))

%% test avg vs max for one day
% figure; imagesc(mean_avg_image_px{3}(:,:,6))
% axis image;
% axis off;
% clim_val = max(clim);
% clim([-clim_val, clim_val]);
% ap.wf_draw('ccf','k');
% colormap(ap.colormap('PWG', [], 1.2));
% colorbar;
% title('Mean')
% 
% 
% figure; imagesc(max_avg_image_px{3}(:,:,6))
% axis image;
% axis off;
% clim_val = max(abs(test2(:))); 
% clim([-clim_val, clim_val]);
% ap.wf_draw('ccf','k');
% colormap(ap.colormap('PWG', [], 1.2));
% colorbar;
% title('Max')

%% mean avg image

days_for_plot = -3:2;
these_days_from_learning = wf_unique_avg_animal_group_indices;
plot_day_idx = ismember(these_days_from_learning, days_for_plot);
plotted_days = these_days_from_learning(plot_day_idx);

for stim_idx=1:length(unique_stims)
    avg_image_px = mean_avg_image_px{stim_idx}(:,:,plot_day_idx);
    
    figure;
    tiledlayout("flow")
    for day_idx=1:length(plotted_days)
        nexttile;
        imagesc(avg_image_px(:, :, day_idx))
        axis image
        ap.wf_draw('ccf','k');
        colormap(ap.colormap('PWG', [], 1.2));
        clim_val = max(abs(clim)) * 0.7;
        clim([-clim_val, clim_val]); % Set the color limits
        clim(clim/2)
%         colorbar_handle = colorbar;
%         colorbar_handle.Ticks = [0, clim_val];
%         colorbar_handle.TickLabels = {'0', 'max'};
        %         colorbar_handle.FontSize = 40;
        axis off;
        title(['Day ' num2str(plotted_days(day_idx))])
    end
    sgtitle(['Avg Stim ' num2str(unique_stims(stim_idx))])
end

%% max avg image

days_for_plot = -3:2;
these_days_from_learning = wf_unique_avg_animal_group_indices;
plot_day_idx = ismember(these_days_from_learning, days_for_plot);
plotted_days = these_days_from_learning(plot_day_idx);

for stim_idx=1:length(unique_stims)
    avg_image_px = max_avg_image_px{stim_idx}(:,:,plot_day_idx);
    
    figure;
    tiledlayout("flow")
    for day_idx=1:length(plotted_days)
        nexttile;
        imagesc(avg_image_px(:, :, day_idx))
        axis image
        ap.wf_draw('ccf','k');
        colormap(ap.colormap('PWG', [], 1.2));
        clim_val = max(abs(clim)) * 0.7;
        clim([-clim_val, clim_val]); % Set the color limits
        clim(clim/2)
%         colorbar_handle = colorbar;
%         colorbar_handle.Ticks = [0, clim_val];
%         colorbar_handle.TickLabels = {'0', 'max'};
        %         colorbar_handle.FontSize = 40;
        axis off;
        title(['Day ' num2str(plotted_days(day_idx))])
    end
    sgtitle(['Max avg Stim ' num2str(unique_stims(stim_idx))])
end

%% OLD
% 
% 
% % Initialize a structure to hold grouped data
% all_new_pfc_avg_grouped = struct();
% 
% % Iterate over all animals
% for animal_idx = 1:numel(ctx_wf.animals)
%     % Get recording data for the current animal
%     recording_data = ctx_wf.recording_data{animal_idx};
%     
%     % Extract the 'days from learning' column
%     days_from_learning = recording_data.days_from_learning;
%     
%     % get U
%     U_master = ctx_wf.U_master{animal_idx};
% 
%     % Iterate over each day and its corresponding vector
%     for day_idx = 1:numel(days_from_learning)
%         this_day = days_from_learning{day_idx};
% 
%         % get the Vs
%         V_avg_contra_stim_align = recording_data.V_avg_contra_stim_align{day_idx};
% 
%         % get new ROI
%         day_new_pfc_avg_contra_roi_stim_align = ap.wf_roi(U_master,V_avg_contra_stim_align,[],[],new_pfc_roi_mask);
%         
%         % Create a valid field name for the day
%         day_field = ['day_', strrep(num2str(this_day), '-', 'neg')];
%         
%         % Check if the day already exists in the grouped data
%         if ~isfield(all_new_pfc_avg_grouped, day_field)
%             all_new_pfc_avg_grouped.(day_field) = {}; % Initialize as a cell array
%         end
%         
%         % Append the vector to the corresponding day
%         all_new_pfc_avg_grouped.(day_field){end+1} = day_new_pfc_avg_contra_roi_stim_align';
%     end
% end
% 
% % !!!!!!!!! get wf stim time
% wf_stim_time = recording_data.wf_stim_time{1};
% 
% % combine and avg for each day
% day_fields = fieldnames(all_new_pfc_avg_grouped);
% all_new_pfc_avg_grouped_combined = struct();
% for day_idx = 1:numel(day_fields)
%     day_field = day_fields{day_idx};
%     all_new_pfc_avg_grouped_combined.(day_field) = horzcat(all_new_pfc_avg_grouped.(day_field){:}); % Concatenate vectors for this day
% end
% 
% 
% % Get all day fields from grouped_combined
% day_fields = fieldnames(all_new_pfc_avg_grouped_combined);
% 
% % Initialize storage for days, averages, and SEMs
% sorted_days = zeros(numel(day_fields), 1);
% mean_new_pfc_avg_grouped = {};
% sems_new_pfc_avg_grouped = {};
% 
% % Process each day
% for i = 1:numel(day_fields)
%     day_field = day_fields{i};
%     
%     % Extract the numerical day value (handle 'neg' in field names)
%     sorted_days(i) = str2double(strrep(strrep(day_field, 'day_neg', '-'), 'day_', ''));
%     
%     % Extract data for the current day (already a matrix)
%     day_matrix = all_new_pfc_avg_grouped_combined.(day_field); % Each column is a variable
%     
%     % Compute average and SEM across columns
%     avg = mean(day_matrix, 2);             
%     sem = std(day_matrix, 0, 2) / sqrt(size(day_matrix, 2)); 
%     
%     % Store results
%     mean_new_pfc_avg_grouped{i} = avg;
%     sems_new_pfc_avg_grouped{i} = sem;
% end
% 
% % Sort results by day
% [sorted_days, sort_idx] = sort(sorted_days);
% mean_new_pfc_avg_grouped = mean_new_pfc_avg_grouped(sort_idx);
% sems_new_pfc_avg_grouped = sems_new_pfc_avg_grouped(sort_idx);
% 
% % pfc 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Number of days to plot
% num_days = numel(sorted_days);
% 
% figure;
% tiledlayout(4,3);  
% 
% % Loop through each day and plot in a separate subplot
% for i = 1:num_days
% 
%     if sorted_days(i) < 0
%         day_field = sprintf('day_neg%d', abs(sorted_days(i)));  % For negative days, use day_neg3, day_neg2, etc.
%     else
%         day_field = sprintf('day_%d', sorted_days(i));  % For non-negative days, use day_3, day_2, etc.
%     end
% 
%     % Data for the current day
%     x = wf_stim_time;
%     avg = mean_new_pfc_avg_grouped{i}';        % Average (1x151)
%     sem = sems_new_pfc_avg_grouped{i}';            % SEM (1x151)
%     
%     % Get the size of the data for the current day (number of repetitions)
%     num_reps = size(all_new_pfc_avg_grouped_combined.(day_field), 2);  % Number of rows (repetitions)
% 
%     % Create a new tile for each day
%     nexttile;
%     hold on;
%     
%     % Create the x and y values for the SEM region
%     x_fill = [x, fliplr(x)];  % X values: time points and reversed time points
%     y_fill = [avg + sem, fliplr(avg - sem)];  % Y values: upper and lower bounds of SEM
%     
%     % Plot the shaded SEM region
%     fill(x_fill, y_fill, [0 0.5 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 
%     % Plot average as a line
%     plot(x, avg, 'LineWidth', 1.5, 'Color', [0 0.5 0]);
% 
%     % Title and labels for each subplot
%     title_str = sprintf('Day %d (Repetitions: %d)', sorted_days(i), num_reps);  % Add repetitions to the title
%     title(title_str);
%     xlabel('Time from stim onset');
%     ylabel('\Delta F/F', 'FontSize', 12);
%     xline(0, 'LineWidth', 2)
%     ylim([-10^(-3) 2*10^(-3)])
%     hold off;
%     xlim([-0.1 0.3])
%     AP_scalebar(0.1, 10^(-3))
%     ap.prettyfig
% end
% 
% %% SPARE
% % figure; imagesc(vis_roi_mask); axis image;
% % title('ROI mask');
% 
% vis_roi_trace = ap.wf_roi(wf_U,wf_V,wf_avg,[],vis_roi_mask);
% 
% figure('Position', [680 460 860 520]);
% plot(wf_t, vis_roi_trace, 'color', [0 0.7 0]);
% % title('Visual ROI fluorescence');
% xline(stimOn_times(contra_good_trials), 'k', 'LineWidth', 2)
% ylim([-0.02, 0.02])
% xlim([60 120])
% xlabel('Time (s)', 'FontSize', 34)
% ylabel('{\Delta}F/F', 'FontSize', 34)
% box off
% AP_scalebar(10, 0.01)
% 
% %% - wf avg trials
% roi_stim_align = interp1(wf_t,vis_roi_trace,time_stimulus)';
% avg_contra_roi_stim_align = mean(roi_stim_align(:, contra_good_trials), 2);
% figure('Position', [680 460 860 520]);
% plot(timevec, avg_contra_roi_stim_align, 'color', [0 0.7 0])
% xline(0, 'k', 'LineWidth', 2)
% xline(0.5, 'k', 'LineWidth', 2)
% xlabel('Time from stim onset (s)', 'FontSize', 34)
% ylabel('{\Delta}F/F', 'FontSize', 34)
% % title('Visual ROI fluorescence');
% box off
% AP_scalebar(0.2, 1*10^-3)
