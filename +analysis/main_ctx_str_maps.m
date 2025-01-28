load('D:\matlab_save\new_all_ctx_maps_and_passive_str_ephys_data');

%%% test
this_test = all_ctx_maps_to_str.recording_data{1,1};
this_test.cortex_kernel_px(:)

data = struct;
for i = 1:5
    data(i).wf = rand(100,30);
end

%%%%%%%%%%

% get number of things from total MUA - after re-run save script
total_depths = 0;
for animal_idx=1:height(all_ctx_maps_to_str)
    ctx_maps_to_str = all_ctx_maps_to_str.recording_data{animal_idx};
    total_depths = total_depths + sum([ctx_maps_to_str.no_depths{:}]);
end
% make big cell of all the maps
all_flattened_cortex_kernel_px = [];
for animal_idx=1:height(all_ctx_maps_to_str)
    ctx_maps_to_str = all_ctx_maps_to_str.recording_data{animal_idx};
    cortex_kernel_px = cat(3, ctx_maps_to_str.cortex_kernel_px{:});
    flattened_cortex_kernel_px = reshape(cortex_kernel_px, [], ...
        size(cortex_kernel_px, 3));
    all_flattened_cortex_kernel_px = cat(1, all_flattened_cortex_kernel_px, flattened_cortex_kernel_px');
end

% run kmeans
rng('default');
rng(0);
num_clusters = 4;
[cluster_ids, centroids, sumd] = kmeans(double(all_flattened_cortex_kernel_px), num_clusters, 'Distance', 'correlation',  'Replicates',5);


% check who the masters are
centroid_images = reshape(centroids, [num_clusters, [size(cortex_kernel_px, 1) size(cortex_kernel_px, 2)]]);
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

%% make big vectors of days from learning and mouse id
for_psth_days_from_learning_cell={};
for_psth_mouse_ids_cell={};
for animal_idx=1:height(all_ctx_maps_to_str)
    days_from_learning = [all_ctx_maps_to_str.recording_data{animal_idx}.days_from_learning{:}];
    no_depths = [all_ctx_maps_to_str.recording_data{animal_idx}.no_depths{:}];
    for_psth_days_from_learning_cell{animal_idx} = repelem(days_from_learning, no_depths);
    for_psth_mouse_ids_cell{animal_idx} = repelem(animal_idx, length(for_psth_days_from_learning_cell{animal_idx}));
end
for_psth_days_from_learning = cat(2, for_psth_days_from_learning_cell{:});
for_psth_mouse_ids = cat(2, for_psth_mouse_ids_cell{:});


%% get psths
all_contra_no_move_psth_stim_align = [];
for animal_idx=1:height(all_ctx_maps_to_str)
    %     sum([all_ctx_maps_to_str.recording_data{animal_idx}.no_depths{:}]);
    ctx_maps_to_str = all_ctx_maps_to_str.recording_data{animal_idx};
    contra_no_move_psth_stim_align = cat(1, ctx_maps_to_str.contra_no_move_psth_stim_align{:});
    all_contra_no_move_psth_stim_align = vertcat(all_contra_no_move_psth_stim_align, contra_no_move_psth_stim_align);
end

all_centre_no_move_psth_stim_align = [];
for animal_idx=1:height(all_ctx_maps_to_str)
    %     sum([all_ctx_maps_to_str.recording_data{animal_idx}.no_depths{:}]);
    ctx_maps_to_str = all_ctx_maps_to_str.recording_data{animal_idx};
    centre_no_move_psth_stim_align = cat(1, ctx_maps_to_str.centre_no_move_psth_stim_align{:});
    all_centre_no_move_psth_stim_align = vertcat(all_centre_no_move_psth_stim_align, centre_no_move_psth_stim_align);
end

all_ipsi_no_move_psth_stim_align = [];
for animal_idx=1:height(all_ctx_maps_to_str)
    %     sum([all_ctx_maps_to_str.recording_data{animal_idx}.no_depths{:}]);
    ctx_maps_to_str = all_ctx_maps_to_str.recording_data{animal_idx};
    ipsi_no_move_psth_stim_align = cat(1, ctx_maps_to_str.ipsi_no_move_psth_stim_align{:});
    all_ipsi_no_move_psth_stim_align = vertcat(all_ipsi_no_move_psth_stim_align, ipsi_no_move_psth_stim_align);
end
% % get norm psth 
% all_norm_no_move_psth_stim_align = [];
% for animal_idx=1:height(all_ctx_maps_to_str)
%     %     sum([all_ctx_maps_to_str.recording_data{animal_idx}.no_depths{:}]);
%     ctx_maps_to_str = all_ctx_maps_to_str.recording_data{animal_idx};
%     norm_no_move_psth_stim_align = cat(1, ctx_maps_to_str.norm_no_move_psth_stim_align{:});
%     all_norm_no_move_psth_stim_align = vertcat(all_norm_no_move_psth_stim_align, norm_no_move_psth_stim_align);
% end

%% group according to animal, cluster and day from learning
% cluster and learning subscripts
[~, ~, cluster_subs] = unique(cluster_ids);
[unique_days_from_learning, ~, learning_day_subs] = unique(for_psth_days_from_learning);
[unique_mouse_ids, ~, mouse_subs] = unique(for_psth_mouse_ids);

% Combine cluster_ids and days into a single matrix for indexing
group_indices = [mouse_subs, cluster_ids, learning_day_subs];

% contra
% Define a custom function to handle the grouping for 3D matrices
contra_group_function = @(indices) {all_contra_no_move_psth_stim_align(indices, :, :)};

% Use accumarray to group the values
contra_grouped_psths = accumarray(group_indices, (1:length(cluster_ids))', [], contra_group_function, {});

% centre
% Define a custom function to handle the grouping for 3D matrices
centre_group_function = @(indices) {all_centre_no_move_psth_stim_align(indices, :, :)};

% Use accumarray to group the values
centre_grouped_psths = accumarray(group_indices, (1:length(cluster_ids))', [], centre_group_function, {});

% ipsi
% Define a custom function to handle the grouping for 3D matrices
ipsi_group_function = @(indices) {all_ipsi_no_move_psth_stim_align(indices, :, :)};

% Use accumarray to group the values
ipsi_grouped_psths = accumarray(group_indices, (1:length(cluster_ids))', [], ipsi_group_function, {});

% test = {contra_grouped_psths(:,1,:)};

%% for each mouse average and normalize and smooth in each group
% smoothing window
gauss_win = gausswin(51, 3)';

% contra
bin_centres = all_ctx_maps_to_str.recording_data{1,1}.bin_centres{1,1};  
contra_smooth_norm_grouped_psths = cell(size(contra_grouped_psths));
for mouse_id = 1:length(unique_mouse_ids)
    this_grouped_psth = squeeze(contra_grouped_psths(mouse_id, :, :));
    this_avg_grouped_psth = cellfun(@(x) sum(x, 1), this_grouped_psth, 'UniformOutput', false);
    % get non-empty elements
    non_empty_mask = cellfun(@(x) ~isempty(x), this_avg_grouped_psth);
    non_empty_grouped_psths = cell2mat(this_avg_grouped_psth(non_empty_mask));
    % normalize and smooth non-empty elements
    psth_baseline = mean(non_empty_grouped_psths(:, bin_centres > -0.2 & bin_centres < 0), 2);
%     psth_std = std(psth_baseline);
    normalized_psths = (non_empty_grouped_psths - psth_baseline) ./ (psth_baseline);
    smooth_norm_psths = filter(gauss_win,sum(gauss_win),normalized_psths, [], 2);
    % update the cell array with normalized values
    this_norm_avg_grouped_psth = this_avg_grouped_psth;
    this_norm_avg_grouped_psth(non_empty_mask) = mat2cell(smooth_norm_psths, ...
        ones(1, size(smooth_norm_psths, 1)), size(smooth_norm_psths, 2));
    % add to big normalized matrix
    contra_smooth_norm_grouped_psths(mouse_id, :, :) = this_norm_avg_grouped_psth;
end

% centre
bin_centres = all_ctx_maps_to_str.recording_data{1,1}.bin_centres{1,1};
centre_smooth_norm_grouped_psths = cell(size(centre_grouped_psths));
for mouse_id = 1:length(unique_mouse_ids)
    this_grouped_psth = squeeze(centre_grouped_psths(mouse_id, :, :));
    this_avg_grouped_psth = cellfun(@(x) sum(x, 1), this_grouped_psth, 'UniformOutput', false);
    % get non-empty elements
    non_empty_mask = cellfun(@(x) ~isempty(x), this_avg_grouped_psth);
    non_empty_grouped_psths = cell2mat(this_avg_grouped_psth(non_empty_mask));
    % normalize non-empty elements
    psth_baseline = mean(non_empty_grouped_psths(:, bin_centres > -0.2 & bin_centres < 0), 2);
%     psth_std = std(psth_baseline);
    normalized_psths = (non_empty_grouped_psths - psth_baseline) ./ (psth_baseline);
    smooth_norm_psths = filter(gauss_win,sum(gauss_win),normalized_psths, [], 2);

    % update the cell array with normalized values
    this_norm_avg_grouped_psth = this_avg_grouped_psth;
    this_norm_avg_grouped_psth(non_empty_mask) = mat2cell(smooth_norm_psths, ...
        ones(1, size(smooth_norm_psths, 1)), size(smooth_norm_psths, 2));
    % add to big normalized matrix
    centre_smooth_norm_grouped_psths(mouse_id, :, :) = this_norm_avg_grouped_psth;
end

% ipsi
bin_centres = all_ctx_maps_to_str.recording_data{1,1}.bin_centres{1,1};
ipsi_smooth_norm_grouped_psths = cell(size(ipsi_grouped_psths));
for mouse_id = 1:length(unique_mouse_ids)
    this_grouped_psth = squeeze(ipsi_grouped_psths(mouse_id, :, :));
    this_avg_grouped_psth = cellfun(@(x) sum(x, 1), this_grouped_psth, 'UniformOutput', false);
    % get non-empty elements
    non_empty_mask = cellfun(@(x) ~isempty(x), this_avg_grouped_psth);
    non_empty_grouped_psths = cell2mat(this_avg_grouped_psth(non_empty_mask));
    % normalize non-empty elements
    psth_baseline = mean(non_empty_grouped_psths(:, bin_centres > -0.2 & bin_centres < 0), 2);
%     psth_std = std(psth_baseline);
    normalized_psths = (non_empty_grouped_psths - psth_baseline) ./ (psth_baseline);
    smooth_norm_psths = filter(gauss_win,sum(gauss_win),normalized_psths, [], 2);

    % update the cell array with normalized values
    this_norm_avg_grouped_psth = this_avg_grouped_psth;
    this_norm_avg_grouped_psth(non_empty_mask) = mat2cell(smooth_norm_psths, ...
        ones(1, size(smooth_norm_psths, 1)), size(smooth_norm_psths, 2));
    % add to big normalized matrix
    ipsi_smooth_norm_grouped_psths(mouse_id, :, :) = this_norm_avg_grouped_psth;
end

%% take these psths and average across mice
contra_avg_smooth_norm_grouped_psths = cell(num_clusters, length(unique_days_from_learning));
for cluster_id = 1:num_clusters
    for day_idx = 1:length(unique_days_from_learning)
        norm_grouped_psth = contra_smooth_norm_grouped_psths(:, cluster_id, day_idx);
        non_empty_mask = ~cellfun(@isempty, norm_grouped_psth);
        non_empty_data = cat(1, norm_grouped_psth{non_empty_mask});
        if size(non_empty_data, 1)>=3
            contra_avg_smooth_norm_grouped_psths{cluster_id, day_idx} = mean(non_empty_data, 1);
        else 
            contra_avg_smooth_norm_grouped_psths{cluster_id, day_idx} = nan;
        end
    end
end

centre_avg_smooth_norm_grouped_psths = cell(num_clusters, length(unique_days_from_learning));
for cluster_id = 1:num_clusters
    for day_idx = 1:length(unique_days_from_learning)
        norm_grouped_psth = centre_smooth_norm_grouped_psths(:, cluster_id, day_idx);
        non_empty_mask = ~cellfun(@isempty, norm_grouped_psth);
        non_empty_data = cat(1, norm_grouped_psth{non_empty_mask});
        if size(non_empty_data, 1)>=3
            centre_avg_smooth_norm_grouped_psths{cluster_id, day_idx} = mean(non_empty_data, 1);
        else
            centre_avg_smooth_norm_grouped_psths{cluster_id, day_idx} = nan;
        end
    end
end

ipsi_avg_smooth_norm_grouped_psths = cell(num_clusters, length(unique_days_from_learning));
for cluster_id = 1:num_clusters
    for day_idx = 1:length(unique_days_from_learning)
        norm_grouped_psth = ipsi_smooth_norm_grouped_psths(:, cluster_id, day_idx);
        non_empty_mask = ~cellfun(@isempty, norm_grouped_psth);
        non_empty_data = cat(1, norm_grouped_psth{non_empty_mask});
        if size(non_empty_data, 1)>=3
            ipsi_avg_smooth_norm_grouped_psths{cluster_id, day_idx} = mean(non_empty_data, 1);
        else 
            ipsi_avg_smooth_norm_grouped_psths{cluster_id, day_idx} = nan;
        end
    end
end



%% get max amplitude across days

max_window = bin_centres>0&bin_centres<0.3;
% % CONTRA
% do max for each mouse
contra_per_mouse_max_ampl = {};
for cluster_id = 1:num_clusters
    for mouse_id = 1:length(unique_mouse_ids)
        this_norm_grouped_psth = squeeze(contra_smooth_norm_grouped_psths(mouse_id, cluster_id, :));
        this_mouse_max_ampl_grouped_psths_cell = cellfun(@(x) ...
            max(x(max_window), [], 'all'), ...
            this_norm_grouped_psth, ...
            'UniformOutput', false, ...
            'ErrorHandler', @(x, varargin) nan);  % Handle errors (e.g., empty cells) by returning NaN
        this_mouse_empty_groups_idx = cellfun(@isempty, this_mouse_max_ampl_grouped_psths_cell);
        this_mouse_max_ampl_grouped_psths_cell(this_mouse_empty_groups_idx) = {nan};
        this_mouse_max_ampl_grouped_psths = cell2mat(this_mouse_max_ampl_grouped_psths_cell);
        contra_per_mouse_max_ampl(mouse_id, cluster_id) = {this_mouse_max_ampl_grouped_psths};
    end
end
contra_max_ampl_grouped_psths = zeros(11, 4);
for cluster_id = 1:size(contra_per_mouse_max_ampl,2)
    concatenated_contra_per_mouse_max_ampl = cat(2, contra_per_mouse_max_ampl{:, cluster_id}); % count how many animals are in each day with non-nan
    contra_max_ampl_grouped_psths(:, cluster_id) = nanmean(concatenated_contra_per_mouse_max_ampl, 2);
end

% CENTRE
% max for each mouse
centre_per_mouse_max_ampl = {};
for cluster_id = 1:num_clusters
    for mouse_id = 1:length(unique_mouse_ids)
        this_norm_grouped_psth = squeeze(centre_smooth_norm_grouped_psths(mouse_id, cluster_id, :));
           this_mouse_max_ampl_grouped_psths_cell = cellfun(@(x) ...
            max(x(max_window), [], 'all'), ...
            this_norm_grouped_psth, ...
            'UniformOutput', false, ...
            'ErrorHandler', @(x, varargin) nan);  % Handle errors (e.g., empty cells) by returning NaN
        this_mouse_empty_groups_idx = cellfun(@isempty, this_mouse_max_ampl_grouped_psths_cell);
        this_mouse_max_ampl_grouped_psths_cell(this_mouse_empty_groups_idx) = {nan};
        this_mouse_max_ampl_grouped_psths = cell2mat(this_mouse_max_ampl_grouped_psths_cell);
        centre_per_mouse_max_ampl(mouse_id, cluster_id) = {this_mouse_max_ampl_grouped_psths};
    end
end
centre_max_ampl_grouped_psths = zeros(11, 4);
for cluster_id = 1:size(centre_per_mouse_max_ampl,2)
    concatenated_centre_per_mouse_max_ampl = cat(2, centre_per_mouse_max_ampl{:, cluster_id});
    centre_max_ampl_grouped_psths(:, cluster_id) = nanmean(concatenated_centre_per_mouse_max_ampl, 2);
end

% IPSI
% max for each mouse
ipsi_per_mouse_max_ampl = {};
for cluster_id = 1:num_clusters
    for mouse_id = 1:length(unique_mouse_ids)
        this_norm_grouped_psth = squeeze(ipsi_smooth_norm_grouped_psths(mouse_id, cluster_id, :));
           this_mouse_max_ampl_grouped_psths_cell = cellfun(@(x) ...
            max(x(max_window), [], 'all'), ...
            this_norm_grouped_psth, ...
            'UniformOutput', false, ...
            'ErrorHandler', @(x, varargin) nan);  % Handle errors (e.g., empty cells) by returning NaN
        this_mouse_empty_groups_idx = cellfun(@isempty, this_mouse_max_ampl_grouped_psths_cell);
        this_mouse_max_ampl_grouped_psths_cell(this_mouse_empty_groups_idx) = {nan};
        this_mouse_max_ampl_grouped_psths = cell2mat(this_mouse_max_ampl_grouped_psths_cell);
        ipsi_per_mouse_max_ampl(mouse_id, cluster_id) = {this_mouse_max_ampl_grouped_psths};
    end
end
ipsi_max_ampl_grouped_psths = zeros(11, 4);
for cluster_id = 1:size(ipsi_per_mouse_max_ampl,2)
    concatenated_ipsi_per_mouse_max_ampl = cat(2, ipsi_per_mouse_max_ampl{:, cluster_id});
    ipsi_max_ampl_grouped_psths(:, cluster_id) = nanmean(concatenated_ipsi_per_mouse_max_ampl, 2);
end

%% plots 

%% ADD HERE PSTH PER MOUSE
% % % contra
% bin_window = 0.001;
% bin_centres = -0.5:bin_window:2;
% for cluster_id=1:num_clusters
%     figure;
%     tiledlayout('flow');
%     for day_idx=1:max(learning_day_subs)
%         if isempty(contra_avg_norm_grouped_psths{cluster_id, day_idx})
%             continue
%         end
%         nexttile;
%         plot(bin_centres, contra_avg_norm_grouped_psths{cluster_id, day_idx})
%         title(['Day ' num2str(unique_days_from_learning(day_idx))]);
%     end
%     sgtitle(['Contra Cluster ' num2str(cluster_id)])
% end
% 
% % centre
% bin_window = 0.001;
% bin_centres = -0.5:bin_window:2;
% for cluster_id=1:num_clusters
%     figure;
%     tiledlayout('flow');
%     for day_idx=1:max(learning_day_subs)
%         if isempty(centre_avg_norm_grouped_psths{cluster_id, day_idx})
%             continue
%         end
%         nexttile;
%         plot(bin_centres, centre_avg_norm_grouped_psths{cluster_id, day_idx})
%         title(['Day ' num2str(unique_days_from_learning(day_idx))]);
%     end
%     sgtitle(['Centre Cluster ' num2str(cluster_id)])
% end


%% psths across days

% use for all
use_unique_days_from_learning = unique_days_from_learning(3:9);
legend_labels_days = cellfun(@(x) sprintf('Day %d', x), ...
         num2cell(use_unique_days_from_learning), 'UniformOutput', false);

% CONTRA
n_days = length(unique_days_from_learning);
contra_good_days = ~any(isnan(contra_max_ampl_grouped_psths), 2);
figure;
tiledlayout(num_clusters,2)
for cluster_id=1:num_clusters
    nexttile;
    imagesc(squeeze(centroid_images(cluster_id,:,:)))
    axis image;
    clim(max(abs(clim)).*[-1,1]*0.7);
    ap.wf_draw('ccf','k');
    colormap(ap.colormap('PWG'));
    ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight','bold');

    all_psths = nexttile;
    legend_contra_good_days = contra_good_days;
    for day_idx=1:n_days
        if contra_good_days(day_idx) == 1
            if isnan(contra_avg_smooth_norm_grouped_psths{cluster_id, day_idx})
                legend_contra_good_days(day_idx) = 0;
                continue
            end
    
            if unique_days_from_learning(day_idx)<0
                plot(bin_centres, contra_avg_smooth_norm_grouped_psths{cluster_id, day_idx}, '--');
            else
                plot(bin_centres, contra_avg_smooth_norm_grouped_psths{cluster_id, day_idx});
            end
            hold on;
        end
    end
    all_psths.ColorOrder = brewermap(n_days,"PuRd");
    xlabel('Time from stim onset (s)')
    ylabel('Firing rate')
%     legend_labels_days = cellfun(@(x) sprintf('Day %d', x), ...
%         num2cell(unique_days_from_learning(logical(legend_contra_good_days))), 'UniformOutput', false);
    legend(legend_labels_days)
    ylim([-1 2.5])
end
sgtitle('Contra Stim', 'FontSize', 20, 'FontWeight','bold')

% CENTRE
n_days = length(unique_days_from_learning);
centre_good_days = ~any(isnan(centre_max_ampl_grouped_psths), 2);
figure;
tiledlayout(num_clusters,2)
for cluster_id=1:num_clusters
    nexttile;
    imagesc(squeeze(centroid_images(cluster_id,:,:)))
    axis image;
    clim(max(abs(clim)).*[-1,1]*0.7);
    ap.wf_draw('ccf','k');
    colormap(ap.colormap('PWG'));
    ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight','bold');

    all_psths = nexttile;
    legend_centre_good_days = centre_good_days;
    for day_idx=1:n_days
        if centre_good_days(day_idx) == 1
            if isnan(centre_avg_smooth_norm_grouped_psths{cluster_id, day_idx})
                legend_centre_good_days(day_idx) = 0;
                continue
            end
    
            if unique_days_from_learning(day_idx)<0
                plot(bin_centres, centre_avg_smooth_norm_grouped_psths{cluster_id, day_idx}, '--');
            else
                plot(bin_centres, centre_avg_smooth_norm_grouped_psths{cluster_id, day_idx});
            end
            hold on;
        end
    end
    all_psths.ColorOrder = brewermap(n_days,"PuRd");
    xlabel('Time from stim onset (s)')
    ylabel('Firing rate')
%     legend_labels_days = cellfun(@(x) sprintf('Day %d', x), ...
%         num2cell(unique_days_from_learning(logical(legend_centre_good_days))), 'UniformOutput', false);
    legend(legend_labels_days)
    ylim([-1 2.5])
end
sgtitle('Centre Stim', 'FontSize', 20, 'FontWeight','bold')

% ipsi
n_days = length(unique_days_from_learning);
ipsi_good_days = ~any(isnan(ipsi_max_ampl_grouped_psths), 2);
figure;
tiledlayout(num_clusters,2)
for cluster_id=1:num_clusters
    nexttile;
    imagesc(squeeze(centroid_images(cluster_id,:,:)))
    axis image;
    clim(max(abs(clim)).*[-1,1]*0.7);
    ap.wf_draw('ccf','k');
    colormap(ap.colormap('PWG'));
    ylabel(['Cluster ', num2str(cluster_id)], 'FontSize', 14, 'FontWeight','bold');

    all_psths = nexttile;
    legend_ipsi_good_days = ipsi_good_days;
    for day_idx=1:n_days
        if ipsi_good_days(day_idx) == 1
            if isnan(ipsi_avg_smooth_norm_grouped_psths{cluster_id, day_idx})
                legend_ipsi_good_days(day_idx) = 0;
                continue
            end
    
            if unique_days_from_learning(day_idx)<0
                plot(bin_centres, ipsi_avg_smooth_norm_grouped_psths{cluster_id, day_idx}, '--');
            else
                plot(bin_centres, ipsi_avg_smooth_norm_grouped_psths{cluster_id, day_idx});
            end
            hold on;
        end
    end
    all_psths.ColorOrder = brewermap(n_days,"PuRd");
    xlabel('Time from stim onset (s)')
    ylabel('Firing rate')
%     legend_labels_days = cellfun(@(x) sprintf('Day %d', x), ...
%         num2cell(unique_days_from_learning(logical(legend_ipsi_good_days))), 'UniformOutput', false);
    legend(legend_labels_days)
    ylim([-1 2.5])
end
sgtitle('Ipsi Stim', 'FontSize', 20, 'FontWeight','bold')
%% make plot of max amplitude next to the master cortical maps
% CONTRA
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
    for mouse_id = 1:length(unique_mouse_ids)
        plot(unique_days_from_learning, contra_per_mouse_max_ampl{mouse_id, cluster_id}, ...
            '-o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5], 'Color', [0.5 0.5 0.5]);
        hold on;
    end
    plot(unique_days_from_learning, contra_max_ampl_grouped_psths(:,cluster_id), ...
        '-o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k');    xlabel('Days from association day')
    ylabel('Increase in firing rate')
    xlim([contra_days_on_plot(1), contra_days_on_plot(end)])
end
sgtitle('Contra Stim', 'FontSize', 20, 'FontWeight','bold')

% CENTRE
figure;
centre_days_on_plot = unique_days_from_learning(centre_good_days);
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
    for mouse_id = 1:length(unique_mouse_ids)
        plot(unique_days_from_learning, centre_per_mouse_max_ampl{mouse_id, cluster_id}, ...
            '-o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5], 'Color', [0.5 0.5 0.5]);
        hold on;
    end
    plot(unique_days_from_learning, centre_max_ampl_grouped_psths(:,cluster_id), ...
        '-o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k');    xlabel('Days from association day')
    ylabel('Increase in firing rate')
    xlim([centre_days_on_plot(1), centre_days_on_plot(end)])
end
sgtitle('Centre Stim', 'FontSize', 20, 'FontWeight','bold')

% IPSI
figure;
ipsi_days_on_plot = unique_days_from_learning(ipsi_good_days);
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
    for mouse_id = 1:length(unique_mouse_ids)
        plot(unique_days_from_learning, ipsi_per_mouse_max_ampl{mouse_id, cluster_id}, ...
            '-o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5], 'Color', [0.5 0.5 0.5]);
        hold on;
    end
    plot(unique_days_from_learning, ipsi_max_ampl_grouped_psths(:,cluster_id), ...
        '-o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k');    xlabel('Days from association day')
    ylabel('Increase in firing rate')
    xlim([ipsi_days_on_plot(1), ipsi_days_on_plot(end)])
end
sgtitle('Ipsi Stim', 'FontSize', 20, 'FontWeight','bold')

% figure;
% tiledlayout('flow');
% for cluster_id=1:num_clusters
%     nexttile;
%     for mouse_id = 1:length(unique_mouse_ids)
%         plot(unique_days_from_learning, per_mouse_max_ampl{mouse_id, cluster_id}, ...
%             '-o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5], 'Color', [0.5 0.5 0.5]);
%         hold on;
%     end
%     plot(unique_days_from_learning, max_ampl_grouped_psths(cluster_id,:), ...
%         '-o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k');
%     title(['Cluster ' num2str(cluster_id)])
% end


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
    for mouse_id = 1:length(unique_mouse_ids)
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
    data_for_sem = zeros(length(unique_mouse_ids), length(unique_days_from_learning));

    for mouse_id = 1:length(unique_mouse_ids)
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
    data_for_sem = zeros(length(unique_mouse_ids), length(unique_days_from_learning));

    for mouse_id = 1:length(unique_mouse_ids)
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
    data_for_sem = zeros(length(unique_mouse_ids), length(unique_days_from_learning));

    for mouse_id = 1:length(unique_mouse_ids)
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
