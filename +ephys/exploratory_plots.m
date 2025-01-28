function exploratory_plots(animal, days_from_learning, ephys_data, save_fig_path, cell_group)

save_fig_path_cell = fullfile(save_fig_path, cell_group);

if ~isfolder(save_fig_path_cell)
    mkdir(save_fig_path_cell)
end

for use_rec=1:height(ephys_data)

    rec_day = ephys_data.rec_day{use_rec};
    bin_centres = ephys_data.bin_centres{use_rec, :};

    good_templates = ephys_data.good_templates{use_rec};
    good_templates_idx = find(good_templates);

    if strcmpi(cell_group,'all')
        units_cell_group = ones(length(good_templates_idx), 1);
    else
        units_cell_group = ephys_data.(lower(cell_group)){use_rec};
    end

    contra_good_trials = ephys_data.contra_good_trials{use_rec};
    contra_all_spikes_in_stim_time = ephys_data.contra_all_spikes_in_stim_time{use_rec};
    contra_sharp_p_units = ephys_data.contra_sharp_p_units{use_rec};
    contra_wide_p_units = ephys_data.contra_wide_p_units{use_rec};

    centre_good_trials = ephys_data.centre_good_trials{use_rec};
    centre_all_spikes_in_stim_time = ephys_data.centre_all_spikes_in_stim_time{use_rec};
    centre_sharp_p_units = ephys_data.centre_sharp_p_units{use_rec};
    centre_wide_p_units = ephys_data.centre_wide_p_units{use_rec};

    %% get resp and unresp cells in this group

    % contra
    contra_sharp_responsive_units = find(units_cell_group & contra_sharp_p_units < 0.05);
    contra_sharp_unresponsive_units = find(units_cell_group & contra_sharp_p_units >= 0.05);

    contra_wide_responsive_units = find(units_cell_group & contra_wide_p_units < 0.05);
    contra_wide_unresponsive_units = find(units_cell_group & contra_wide_p_units >= 0.05);

    % centre
    centre_sharp_responsive_units = find(units_cell_group & centre_sharp_p_units < 0.05);
    centre_sharp_unresponsive_units = find(units_cell_group & centre_sharp_p_units >= 0.05);

    centre_wide_responsive_units = find(units_cell_group & centre_wide_p_units < 0.05);
    centre_wide_unresponsive_units = find(units_cell_group & centre_wide_p_units >= 0.05);

    %% - odd vs even

    % contra
    odd_contra_trials = 1:2:sum(contra_good_trials);
    even_contra_trials = 2:2:sum(contra_good_trials);

    odd_contra_mean_all_spikes_in_stim_time = squeeze(mean(contra_all_spikes_in_stim_time(:,odd_contra_trials,:), 2));
    even_contra_mean_all_spikes_in_stim_time = squeeze(mean(contra_all_spikes_in_stim_time(:,even_contra_trials,:), 2));

    % centre
    odd_centre_trials = 1:2:sum(centre_good_trials);
    even_centre_trials = 2:2:sum(centre_good_trials);

    odd_centre_mean_all_spikes_in_stim_time = squeeze(mean(centre_all_spikes_in_stim_time(:,odd_centre_trials,:), 2));
    even_centre_mean_all_spikes_in_stim_time = squeeze(mean(centre_all_spikes_in_stim_time(:,even_centre_trials,:), 2));

    %% - smooth
    % define gaussian window
    gauss_win = gausswin(51, 3)';

    % contra all units
    even_contra_smooth_all_units = filter(gauss_win,sum(gauss_win),even_contra_mean_all_spikes_in_stim_time, [], 2);
    odd_contra_smooth_all_units = filter(gauss_win,sum(gauss_win),odd_contra_mean_all_spikes_in_stim_time, [], 2);

    % centre all units
    even_centre_smooth_all_units = filter(gauss_win,sum(gauss_win),even_centre_mean_all_spikes_in_stim_time, [], 2);
    odd_centre_smooth_all_units = filter(gauss_win,sum(gauss_win),odd_centre_mean_all_spikes_in_stim_time, [], 2);


    %% - get baseline

    %     contra_baseline = mean(mean(contra_all_spikes_in_stim_time(:,:,bin_centres>-0.2 & bin_centres<0), 3), 2);
    %     centre_baseline = mean(mean(centre_all_spikes_in_stim_time(:,:,bin_centres>-0.2 & bin_centres<0), 3), 2);

    even_contra_smooth_baseline = mean(even_contra_smooth_all_units(:,bin_centres>-0.2 & bin_centres<0), 2);
    odd_contra_smooth_baseline = mean(odd_contra_smooth_all_units(:,bin_centres>-0.2 & bin_centres<0), 2);

    even_centre_smooth_baseline = mean(even_centre_smooth_all_units(:,bin_centres>-0.2 & bin_centres<0), 2);
    odd_centre_smooth_baseline = mean(odd_centre_smooth_all_units(:,bin_centres>-0.2 & bin_centres<0), 2);

    %% - normalize

    % contra
    even_contra_norm_smooth_all_units = (even_contra_smooth_all_units - even_contra_smooth_baseline) ...
        ./ (even_contra_smooth_baseline + std(even_contra_smooth_baseline));
    odd_contra_norm_smooth_all_units = (odd_contra_smooth_all_units - odd_contra_smooth_baseline) ...
        ./ (odd_contra_smooth_baseline + std(odd_contra_smooth_baseline));

    % centre
    even_centre_norm_smooth_all_units = (even_centre_smooth_all_units - even_centre_smooth_baseline) ...
        ./ (even_centre_smooth_baseline + std(even_centre_smooth_baseline));
    odd_centre_norm_smooth_all_units = (odd_centre_smooth_all_units - odd_centre_smooth_baseline) ...
        ./ (odd_centre_smooth_baseline + std(odd_centre_smooth_baseline));

    %% all trials 

    % contra
    all_contra_mean_all_spikes_in_stim_time = squeeze(mean(contra_all_spikes_in_stim_time, 2));

    % - smooth
    % define gaussian window
    gauss_win = gausswin(51, 3)';

    % contra all units
    all_contra_smooth_all_units = filter(gauss_win,sum(gauss_win),all_contra_mean_all_spikes_in_stim_time, [], 2);

    % - get baseline

    all_contra_smooth_baseline = mean(all_contra_smooth_all_units(:,bin_centres>-0.2 & bin_centres<0), 2);

    % - normalize
    all_contra_norm_smooth_all_units = (all_contra_smooth_all_units - all_contra_smooth_baseline) ...
        ./ (all_contra_smooth_baseline + std(all_contra_smooth_baseline));

    
    % centre
    all_centre_mean_all_spikes_in_stim_time = squeeze(mean(centre_all_spikes_in_stim_time, 2));

    % - smooth
    % define gaussian window
    gauss_win = gausswin(51, 3)';

    % centre all units
    all_centre_smooth_all_units = filter(gauss_win,sum(gauss_win),all_centre_mean_all_spikes_in_stim_time, [], 2);

    % - get baseline

    all_centre_smooth_baseline = mean(all_centre_smooth_all_units(:,bin_centres>-0.2 & bin_centres<0), 2);

    % - normalize
    all_centre_norm_smooth_all_units = (all_centre_smooth_all_units - all_centre_smooth_baseline) ...
        ./ (all_centre_smooth_baseline + std(all_centre_smooth_baseline));

    %% get number of resp cells
    % contra
    num_cells_contra_sharp_responsive_units(use_rec) = length(contra_sharp_responsive_units);
    num_cells_contra_wide_responsive_units(use_rec) = length(contra_wide_responsive_units);

    % centre
    num_cells_centre_sharp_responsive_units(use_rec) = length(centre_sharp_responsive_units);
    num_cells_centre_wide_responsive_units(use_rec) = length(centre_wide_responsive_units);

    %% get fraction of resp cells
    % contra
    frac_cells_contra_sharp_responsive_units(use_rec) = length(contra_sharp_responsive_units) / length(units_cell_group);
    frac_cells_contra_wide_responsive_units(use_rec) = length(contra_wide_responsive_units) / length(units_cell_group);

    % centre
    frac_cells_centre_sharp_responsive_units(use_rec) = length(centre_sharp_responsive_units) / length(units_cell_group);
    frac_cells_centre_wide_responsive_units(use_rec) = length(centre_wide_responsive_units) / length(units_cell_group);

    %% max amplitude
    % contra
    contra_sharp_mean_max_amplitude = max(max(abs(mean(even_contra_norm_smooth_all_units(contra_sharp_responsive_units, :), 1))), ...
        max(abs(mean(odd_contra_norm_smooth_all_units(contra_sharp_responsive_units, :), 1))));
    contra_wide_mean_max_amplitude = max(max(abs(mean(even_contra_norm_smooth_all_units(contra_wide_responsive_units, :), 1))), ...
        max(abs(mean(odd_contra_norm_smooth_all_units(contra_wide_responsive_units, :), 1))));
    all_contra_sharp_mean_max_amplitude(use_rec) = contra_sharp_mean_max_amplitude;
    all_contra_wide_mean_max_amplitude(use_rec) = contra_wide_mean_max_amplitude;

    % centre
    centre_sharp_mean_max_amplitude = max(max(abs(mean(even_centre_norm_smooth_all_units(centre_sharp_responsive_units, :), 1))), ...
        max(abs(mean(odd_centre_norm_smooth_all_units(centre_sharp_responsive_units, :), 1))));
    centre_wide_mean_max_amplitude = max(max(abs(mean(even_centre_norm_smooth_all_units(centre_wide_responsive_units, :), 1))), ...
        max(abs(mean(odd_centre_norm_smooth_all_units(centre_wide_responsive_units, :), 1))));
    all_centre_sharp_mean_max_amplitude(use_rec) = centre_sharp_mean_max_amplitude;
    all_centre_wide_mean_max_amplitude(use_rec) = centre_wide_mean_max_amplitude;


    %% - sort units on even trials based on normalized values in odd trials

    % get spikes 50-150ms post stim
    post_stim_time = [0.05 0.15];
    post_bin_window = diff(post_stim_time);
    this_post_stim_time = bin_centres>post_stim_time(1) & bin_centres<post_stim_time(2);

    % contra
    odd_contra_post_stim = mean(odd_contra_norm_smooth_all_units(:, this_post_stim_time) / post_bin_window, 2);

    % sort for plotting
    [contra_sorted_post_stim, contra_sorted_units] = sort(odd_contra_post_stim, 'descend');
    even_contra_sorted_norm_smooth_all_units = even_contra_norm_smooth_all_units(contra_sorted_units, :);

    % get sharp responsive sorted units for plotting
    [~, contra_sorted_sharp_responsive_units_idx] = sort(odd_contra_post_stim(contra_sharp_responsive_units), 'descend');
    contra_sorted_sharp_responsive_units = contra_sharp_responsive_units(contra_sorted_sharp_responsive_units_idx);
    even_contra_sorted_norm_smooth_sharp_responsive_units = even_contra_norm_smooth_all_units(contra_sorted_sharp_responsive_units, :);

    % get sharp unresponsive sorted units for plotting
    [~, contra_sharp_unresponsive_units_idx] = sort(odd_contra_post_stim(contra_sharp_unresponsive_units), 'descend');
    contra_sorted_sharp_unresponsive_units = contra_sharp_unresponsive_units(contra_sharp_unresponsive_units_idx);
    even_contra_sorted_norm_smooth_sharp_unresponsive_units = even_contra_norm_smooth_all_units(contra_sorted_sharp_unresponsive_units, :);

    % get wide responsive sorted units for plotting
    [~, contra_sorted_wide_responsive_units_idx] = sort(odd_contra_post_stim(contra_wide_responsive_units), 'descend');
    contra_sorted_wide_responsive_units = contra_wide_responsive_units(contra_sorted_wide_responsive_units_idx);
    even_contra_sorted_norm_smooth_wide_responsive_units = even_contra_norm_smooth_all_units(contra_sorted_wide_responsive_units, :);

    % get wide unresponsive sorted units for plotting
    [~, contra_sorted_wide_unresponsive_units_idx] = sort(odd_contra_post_stim(contra_wide_unresponsive_units), 'descend');
    contra_sorted_wide_unresponsive_units = contra_wide_unresponsive_units(contra_sorted_wide_unresponsive_units_idx);
    even_contra_sorted_norm_smooth_wide_unresponsive_units = even_contra_norm_smooth_all_units(contra_sorted_wide_unresponsive_units, :);

    % get psths for plotting
    all_contra_sorted_norm_smooth_sharp_responsive_units = all_contra_norm_smooth_all_units(contra_sorted_sharp_responsive_units,:);
    all_contra_sorted_norm_smooth_wide_responsive_units = all_contra_norm_smooth_all_units(contra_sorted_wide_responsive_units,:);
    all_contra_sorted_norm_smooth_sharp_unresponsive_units = all_contra_norm_smooth_all_units(contra_sorted_sharp_unresponsive_units,:);
    all_contra_sorted_norm_smooth_wide_unresponsive_units = all_contra_norm_smooth_all_units(contra_sorted_wide_unresponsive_units,:);

    % centre
    odd_centre_post_stim = mean(odd_centre_norm_smooth_all_units(:, this_post_stim_time) / post_bin_window, 2);

    % sort for plotting
    [centre_sorted_post_stim, centre_sorted_units] = sort(odd_centre_post_stim, 'descend');
    even_centre_sorted_norm_smooth_all_units = even_centre_norm_smooth_all_units(centre_sorted_units, :);

    % get sharp responsive sorted units for plotting
    [~, centre_sorted_sharp_responsive_units_idx] = sort(odd_centre_post_stim(centre_sharp_responsive_units), 'descend');
    centre_sorted_sharp_responsive_units = centre_sharp_responsive_units(centre_sorted_sharp_responsive_units_idx);
    even_centre_sorted_norm_smooth_sharp_responsive_units = even_centre_norm_smooth_all_units(centre_sorted_sharp_responsive_units, :);

    % get sharp unresponsive sorted units for plotting
    [~, centre_sharp_unresponsive_units_idx] = sort(odd_centre_post_stim(centre_sharp_unresponsive_units), 'descend');
    centre_sorted_sharp_unresponsive_units = centre_sharp_unresponsive_units(centre_sharp_unresponsive_units_idx);
    even_centre_sorted_norm_smooth_sharp_unresponsive_units = even_centre_norm_smooth_all_units(centre_sorted_sharp_unresponsive_units, :);

    % get wide responsive sorted units for plotting
    [~, centre_sorted_wide_responsive_units_idx] = sort(odd_centre_post_stim(centre_wide_responsive_units), 'descend');
    centre_sorted_wide_responsive_units = centre_wide_responsive_units(centre_sorted_wide_responsive_units_idx);
    even_centre_sorted_norm_smooth_wide_responsive_units = even_centre_norm_smooth_all_units(centre_sorted_wide_responsive_units, :);

    % get wide unresponsive sorted units for plotting
    [~, centre_sorted_wide_unresponsive_units_idx] = sort(odd_centre_post_stim(centre_wide_unresponsive_units), 'descend');
    centre_sorted_wide_unresponsive_units = centre_wide_unresponsive_units(centre_sorted_wide_unresponsive_units_idx);
    even_centre_sorted_norm_smooth_wide_unresponsive_units = even_centre_norm_smooth_all_units(centre_sorted_wide_unresponsive_units, :);

    % get psths for plotting
    all_centre_sorted_norm_smooth_sharp_responsive_units = all_centre_norm_smooth_all_units(centre_sorted_sharp_responsive_units,:);
    all_centre_sorted_norm_smooth_wide_responsive_units = all_centre_norm_smooth_all_units(centre_sorted_wide_responsive_units,:);
    all_centre_sorted_norm_smooth_sharp_unresponsive_units = all_centre_norm_smooth_all_units(centre_sorted_sharp_unresponsive_units,:);
    all_centre_sorted_norm_smooth_wide_unresponsive_units = all_centre_norm_smooth_all_units(centre_sorted_wide_unresponsive_units,:);


    %% - plots per day
    % - contra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contra_stim_fig = figure('Position', get(0, 'Screensize'), Visible="off");
    sgtitle([animal ' ' cell_group ' Contra responses ' rec_day]);

    % -- sharp
    if ~isempty(max(abs(even_contra_sorted_norm_smooth_sharp_responsive_units), [],"all"))
        upper_caxis = max(max(abs(even_contra_sorted_norm_smooth_sharp_responsive_units), [],"all"), ...
            max(abs(even_contra_sorted_norm_smooth_sharp_unresponsive_units), [],"all"));
    else
        upper_caxis = max(abs(even_contra_sorted_norm_smooth_sharp_unresponsive_units), [],"all");
    end

    % responsive cells
    resp = subplot(2,3,1);
    imagesc(bin_centres, [], even_contra_sorted_norm_smooth_sharp_responsive_units);
    %     yticks(1:length(contra_sorted_sharp_responsive_units));
    %     yticklabels(good_templates_idx(contra_sorted_sharp_responsive_units));
    xline(0, 'LineWidth', 1);
    xline(0.5, 'LineWidth', 1);
    colormap(resp, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;
    title('Responsive cells');
    ylabel('Sharp', 'FontWeight', 'bold', 'FontSize', 14);


    % unresponsive cells
    unresp = subplot(2,3,2);
    imagesc(bin_centres, [], even_contra_sorted_norm_smooth_sharp_unresponsive_units);
    %     yticks(1:length(contra_sorted_sharp_unresponsive_units));
    %     yticklabels(good_templates_idx(contra_sorted_sharp_unresponsive_units));
    xline(0, 'LineWidth', 1);
    xline(0.5, 'LineWidth', 1);
    colormap(unresp, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;
    title('Unresponsive cells');

    % mean traces
    subplot(2,3,3);
    plot(bin_centres, mean(even_contra_sorted_norm_smooth_sharp_responsive_units, 1));
    hold on;
    plot(bin_centres, mean(even_contra_sorted_norm_smooth_sharp_unresponsive_units, 1));
    hold on;
    xline(0, 'LineWidth', 1);
    xline(0.5, 'LineWidth', 1);
    legend({'Responsive', 'Unresponsive'});


    % -- wide
    if ~isempty(max(abs(even_contra_sorted_norm_smooth_wide_responsive_units), [],"all"))
        upper_caxis = max(max(abs(even_contra_sorted_norm_smooth_wide_responsive_units), [],"all"), ...
            max(abs(even_contra_sorted_norm_smooth_wide_unresponsive_units), [],"all"));
    else
        upper_caxis = max(abs(even_contra_sorted_norm_smooth_wide_unresponsive_units), [],"all");
    end

    % responsive cells
    resp = subplot(2,3,4);
    imagesc(bin_centres, [], even_contra_sorted_norm_smooth_wide_responsive_units);
    %     yticks(1:length(contra_sorted_wide_responsive_units));
    %     yticklabels(good_templates_idx(contra_sorted_wide_responsive_units));
    xline(0, 'LineWidth', 1);
    xline(0.5, 'LineWidth', 1);
    colormap(resp, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;
    title('Responsive cells');
    ylabel('Wide', 'FontWeight', 'bold', 'FontSize', 14);


    % unresponsive cells
    unresp = subplot(2,3,5);
    imagesc(bin_centres, [], even_contra_sorted_norm_smooth_wide_unresponsive_units);
    %     yticks(1:length(contra_sorted_wide_unresponsive_units));
    %     yticklabels(good_templates_idx(contra_sorted_wide_unresponsive_units));
    xline(0, 'LineWidth', 1);
    xline(0.5, 'LineWidth', 1);
    colormap(unresp, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;
    title('Unresponsive cells');

    % mean traces
    subplot(2,3,6);
    plot(bin_centres, mean(even_contra_sorted_norm_smooth_wide_responsive_units, 1));
    hold on;
    plot(bin_centres, mean(even_contra_sorted_norm_smooth_wide_unresponsive_units, 1));
    hold on;
    xline(0, 'LineWidth', 1);
    xline(0.5, 'LineWidth', 1);
    legend({'Responsive', 'Unresponsive'});

    % save this fig
    contra_stim_fig_name = [animal '_' rec_day '_' cell_group '_contra.tif'];
    contra_stim_fig_path = fullfile(save_fig_path_cell, contra_stim_fig_name);
    saveas(contra_stim_fig, contra_stim_fig_path);

  %%%%%%%%%%%%% TEMP %%%%%%%%%%%%%%%%%%%%%%%%

%   figure;
%   subplot(2,2,1)
%   plot(bin_centres, mean(even_contra_sorted_norm_smooth_wide_responsive_units, 1));
%   hold on;
%   plot(bin_centres, mean(even_contra_sorted_norm_smooth_wide_unresponsive_units, 1));
%   hold on;
%   xline(0, 'LineWidth', 1);
%   xline(0.5, 'LineWidth', 1);
%   legend({'Responsive', 'Unresponsive'});
%   ylabel('Wide', 'FontWeight', 'bold', 'FontSize', 14);
%   title('Even trials', 'FontWeight', 'bold', 'FontSize', 14);
% 
%   subplot(2,2,2)
%   plot(bin_centres, mean(all_contra_sorted_norm_smooth_wide_responsive_units, 1));
%   hold on;
%   plot(bin_centres, mean(all_contra_sorted_norm_smooth_wide_unresponsive_units, 1));
%   hold on;
%   xline(0, 'LineWidth', 1);
%   xline(0.5, 'LineWidth', 1);
%   legend({'Responsive', 'Unresponsive'});
%   title('All trials', 'FontWeight', 'bold', 'FontSize', 14);
% 
%   subplot(2,2,3);
%   plot(bin_centres, mean(even_contra_sorted_norm_smooth_sharp_responsive_units, 1));
%   hold on;
%   plot(bin_centres, mean(even_contra_sorted_norm_smooth_sharp_unresponsive_units, 1));
%   hold on;
%   xline(0, 'LineWidth', 1);
%   xline(0.5, 'LineWidth', 1);
%   legend({'Responsive', 'Unresponsive'});
%   ylabel('Sharp', 'FontWeight', 'bold', 'FontSize', 14);
% 
%   subplot(2,2,4);
%   plot(bin_centres, mean(all_contra_sorted_norm_smooth_sharp_responsive_units, 1));
%   hold on;
%   plot(bin_centres, mean(all_contra_sorted_norm_smooth_sharp_unresponsive_units, 1));
%   hold on;
%   xline(0, 'LineWidth', 1);
%   xline(0.5, 'LineWidth', 1);
%   legend({'Responsive', 'Unresponsive'});


    % - centre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    centre_stim_fig = figure('Position', get(0, 'Screensize'), Visible="off");
    sgtitle([animal ' ' cell_group ' Centre responses ' rec_day]);

    % -- sharp
    if ~isempty(max(abs(even_centre_sorted_norm_smooth_sharp_responsive_units), [],"all"))
        upper_caxis = max(max(abs(even_centre_sorted_norm_smooth_sharp_responsive_units), [],"all"), ...
            max(abs(even_centre_sorted_norm_smooth_sharp_unresponsive_units), [],"all"));
    else
        upper_caxis = max(abs(even_centre_sorted_norm_smooth_sharp_unresponsive_units), [],"all");
    end

    % responsive cells
    resp = subplot(2,3,1);
    imagesc(bin_centres, [], even_centre_sorted_norm_smooth_sharp_responsive_units);
    %     yticks(1:length(centre_sorted_sharp_responsive_units));
    %     yticklabels(good_templates_idx(centre_sorted_sharp_responsive_units));
    xline(0, 'LineWidth', 1);
    xline(0.5, 'LineWidth', 1);
    colormap(resp, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;
    title('Responsive cells');
    ylabel('Sharp', 'FontWeight', 'bold', 'FontSize', 14);


    % unresponsive cells
    unresp = subplot(2,3,2);
    imagesc(bin_centres, [], even_centre_sorted_norm_smooth_sharp_unresponsive_units);
    %     yticks(1:length(centre_sorted_sharp_unresponsive_units));
    %     yticklabels(good_templates_idx(centre_sorted_sharp_unresponsive_units));
    xline(0, 'LineWidth', 1);
    xline(0.5, 'LineWidth', 1);
    colormap(unresp, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;
    title('Unresponsive cells');

    % mean traces
    subplot(2,3,3);
    plot(bin_centres, mean(even_centre_sorted_norm_smooth_sharp_responsive_units, 1));
    hold on;
    plot(bin_centres, mean(even_centre_sorted_norm_smooth_sharp_unresponsive_units, 1));
    hold on;
    xline(0, 'LineWidth', 1);
    xline(0.5, 'LineWidth', 1);
    legend({'Responsive', 'Unresponsive'});


    % -- wide
    if ~isempty(max(abs(even_centre_sorted_norm_smooth_wide_responsive_units), [],"all"))
        upper_caxis = max(max(abs(even_centre_sorted_norm_smooth_wide_responsive_units), [],"all"), ...
            max(abs(even_centre_sorted_norm_smooth_wide_unresponsive_units), [],"all"));
    else
        upper_caxis = max(abs(even_centre_sorted_norm_smooth_wide_unresponsive_units), [],"all");
    end

    % responsive cells
    resp = subplot(2,3,4);
    imagesc(bin_centres, [], even_centre_sorted_norm_smooth_wide_responsive_units);
    %     yticks(1:length(centre_sorted_wide_responsive_units));
    %     yticklabels(good_templates_idx(centre_sorted_wide_responsive_units));
    xline(0, 'LineWidth', 1);
    xline(0.5, 'LineWidth', 1);
    colormap(resp, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;
    title('Responsive cells');
    ylabel('Wide', 'FontWeight', 'bold', 'FontSize', 14);


    % unresponsive cells
    unresp = subplot(2,3,5);
    imagesc(bin_centres, [], even_centre_sorted_norm_smooth_wide_unresponsive_units);
    %     yticks(1:length(centre_sorted_wide_unresponsive_units));
    %     yticklabels(good_templates_idx(centre_sorted_wide_unresponsive_units));
    xline(0, 'LineWidth', 1);
    xline(0.5, 'LineWidth', 1);
    colormap(unresp, AP_colormap('BWR', [], 0.7));
    if upper_caxis
        caxis([-upper_caxis upper_caxis]);
    end
    colorbar;
    title('Unresponsive cells');

    % mean traces
    subplot(2,3,6);
    plot(bin_centres, mean(even_centre_sorted_norm_smooth_wide_responsive_units, 1));
    hold on;
    plot(bin_centres, mean(even_centre_sorted_norm_smooth_wide_unresponsive_units, 1));
    hold on;
    xline(0, 'LineWidth', 1);
    xline(0.5, 'LineWidth', 1);
    legend({'Responsive', 'Unresponsive'});

    % save this fig
    centre_stim_fig_name = [animal '_' rec_day '_' cell_group '_centre.tif'];
    centre_stim_fig_path = fullfile(save_fig_path_cell, centre_stim_fig_name);
    saveas(centre_stim_fig, centre_stim_fig_path);


    %% psth plots  

    % get resp cells in this group
    % contra templates
    contra_sharp_responsive_templates = good_templates_idx(contra_sharp_p_units < 0.05);
    contra_wide_responsive_templates = good_templates_idx(contra_wide_p_units < 0.05);

    % find template at borders
    template_depths = ephys_data.template_depths{use_rec};
    [sorted_template_depths, sorted_template_depth_idx] = sort(template_depths);

    % get str depth
    str_depth = ephys_data.str_depth{use_rec};
    str_start = str_depth(1);
    str_end = str_depth(2);

    % define str split
    str_split_step = floor((str_end-str_start)/3);
    str_split = [str_start ...
        str_start+str_split_step ...
        str_start+2*str_split_step ...
        str_end];

    % find str split and get psths
    template_str_1_idx = [find(sorted_template_depths>str_split(1) & sorted_template_depths<=str_split(2))];
    str_1_templates = sorted_template_depth_idx(template_str_1_idx);
    psth_str_1_stim = mean(all_contra_norm_smooth_all_units(str_1_templates,:), 1);

    template_str_2_idx = [find(sorted_template_depths>str_split(2) & sorted_template_depths<=str_split(3))];
    str_2_templates = sorted_template_depth_idx(template_str_2_idx);
    psth_str_2_stim = mean(all_contra_norm_smooth_all_units(str_2_templates,:), 1);

    template_str_3_idx = [find(sorted_template_depths>str_split(3) & sorted_template_depths<=str_split(4))];
    str_3_templates = sorted_template_depth_idx(template_str_3_idx);
    psth_str_3_stim = mean(all_contra_norm_smooth_all_units(str_3_templates,:), 1);

    % get psth for all str and save out
    psths_str_all_stim{use_rec} = mean([psth_str_1_stim; psth_str_2_stim; psth_str_3_stim], 1);

    % cortex
    cortex_templates = sorted_template_depth_idx(1:template_str_1_idx(1));
    psth_cortex_stim = mean(all_contra_norm_smooth_all_units(cortex_templates,:), 1);

    % under striatum
    if isempty(template_str_3_idx)
        psth_under_str_stim = nan(1, length(psth_str_1_stim));
    else
        under_str_templates = sorted_template_depth_idx(template_str_3_idx(end)+1:end, 1);
        psth_under_str_stim = mean(all_contra_norm_smooth_all_units(under_str_templates,:), 1);
    end

    % append all together
    all_psths_stim = vertcat(psth_cortex_stim, ...
        psth_str_1_stim, ...
        psth_str_2_stim, ...
        psth_str_3_stim, ...
        psth_under_str_stim);

    %         %plot to check
    %     figure;
    %     sgtitle(['Task ' rec_day]);
    %
    % %     stim = subplot(1,3,1);
    %     title('Stim')
    %     plot(psth_str_1_stim)
    %     hold on;
    %     plot(psth_str_2_stim)
    %     hold on;
    %     plot(psth_str_3_stim)
    %     legend({'1', '2', '3'})

    % figure
    depth_psth_per_day_fig = figure('Position', get(0, 'Screensize'), Visible='off');
    plot_depth = tiledlayout(3,4,'TileSpacing','Compact','Padding','Compact');

    % get spike_templates from struct
    spike_templates = ephys_data.spike_templates{use_rec};

    % plot by depth
    unit_axes = nexttile(1,[3 1]);;
    set(unit_axes,'YDir','reverse');
    hold on;

    norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
    unit_dots = scatter3(norm_spike_n,template_depths(unique(spike_templates)), ...
        unique(spike_templates),20,'k','filled');

    sharp_responsive_unit_dots = scatter3(norm_spike_n(contra_sharp_responsive_units),template_depths(contra_sharp_responsive_templates), ...
        contra_sharp_responsive_templates,20,'magenta','filled');

    wide_responsive_unit_dots = scatter3(norm_spike_n(contra_wide_responsive_units),template_depths(contra_wide_responsive_templates), ...
        contra_wide_responsive_templates,20,'blue','filled');

    both_responsive_templates = intersect(contra_wide_responsive_templates, contra_sharp_responsive_templates);
    both_responsive_units = intersect(contra_wide_responsive_units, contra_sharp_responsive_units);

    both_responsive_unit_dots = scatter3(norm_spike_n(both_responsive_units),template_depths(both_responsive_templates), ...
        both_responsive_templates,20,'green','filled');

    yline(str_depth, 'red')
    yline(str_split(2:3), 'blue')
    xlim(unit_axes,[-0.1,1]);
    ylim([-50, max(template_depths)+50]);
    ylabel('Depth (\mum)')
    xlabel('Normalized log rate')
    sgtitle(rec_day, 'FontWeight', 'b', 'FontSize', 14)

    legend([sharp_responsive_unit_dots wide_responsive_unit_dots both_responsive_unit_dots], {'Sharp resp', 'Wide resp', 'Both resp'})

    % plot to check
    nexttile(2, [3 3]);
    plot(bin_centres, psth_str_1_stim)
    hold on;
    plot(bin_centres, psth_str_2_stim)
    hold on;
    plot(bin_centres, psth_str_3_stim)
    xlabel('Time from stim onset')
    legend({'str 1', 'str 2', 'str 3'})
    title(['PSTHs']);

    % save 
    depth_psth_per_day_fig_name = [animal '_' rec_day '_' cell_group '_depth_psth_per_day.tif'];
    depth_psth_per_day_fig_path = fullfile(save_fig_path_cell, depth_psth_per_day_fig_name);
    saveas(depth_psth_per_day_fig, depth_psth_per_day_fig_path);

end

%% plots across days
%% plots

if height(ephys_data)>1
    % num cells
    contra_num_cells_fig = figure(Visible="off");
    plot(days_from_learning, num_cells_contra_sharp_responsive_units, '-o')
    hold on;
    plot(days_from_learning, num_cells_contra_wide_responsive_units, '-o')
    xlim([days_from_learning(1) days_from_learning(end)])
    xline(0, 'LineWidth', 2, 'Color', 'k')
    legend({'Sharp', 'Wide'})
    title([animal ' ' cell_group ' Contra Number of responsive cells across days'])
    contra_stim_fig_name = [animal '_' cell_group '_contra_stim_Number_responsive_cells.tif'];
    contra_stim_fig_path = fullfile(save_fig_path_cell, contra_stim_fig_name);
    saveas(contra_num_cells_fig, contra_stim_fig_path);

    centre_num_cells_fig = figure(Visible="off");
    plot(days_from_learning, num_cells_centre_sharp_responsive_units, '-o')
    hold on;
    plot(days_from_learning, num_cells_centre_wide_responsive_units, '-o')
    xlim([days_from_learning(1) days_from_learning(end)])
    xline(0, 'LineWidth', 2, 'Color', 'k')
    legend({'Sharp', 'Wide'})
    title([animal ' ' cell_group ' Centre Number of responsive cells across days'])
    centre_stim_fig_name = [animal '_' cell_group '_centre_stim_Number_responsive_cells.tif'];
    centre_stim_fig_path = fullfile(save_fig_path_cell, centre_stim_fig_name);
    saveas(centre_num_cells_fig, centre_stim_fig_path);

    % frac cells
    contra_frac_cells_fig = figure(Visible="off");
    plot(days_from_learning, frac_cells_contra_sharp_responsive_units, '-o')
    hold on;
    plot(days_from_learning, frac_cells_contra_wide_responsive_units, '-o')
    xlim([days_from_learning(1) days_from_learning(end)])
    xline(0, 'LineWidth', 2, 'Color', 'k')
    legend({'Sharp', 'Wide'})
    title([animal ' ' cell_group ' Contra Fraction of responsive cells across days'])
    contra_stim_fig_name = [animal '_' cell_group '_contra_stim_Fraction_responsive_cells.tif'];
    contra_stim_fig_path = fullfile(save_fig_path_cell, contra_stim_fig_name);
    saveas(contra_frac_cells_fig, contra_stim_fig_path);

    centre_frac_cells_fig = figure(Visible="off");
    plot(days_from_learning, frac_cells_centre_sharp_responsive_units, '-o')
    hold on;
    plot(days_from_learning, frac_cells_centre_wide_responsive_units, '-o')
    xlim([days_from_learning(1) days_from_learning(end)])
    xline(0, 'LineWidth', 2, 'Color', 'k')
    legend({'Sharp', 'Wide'})
    title([animal ' ' cell_group ' Centre Fraction of responsive cells across days'])
    centre_stim_fig_name = [animal '_' cell_group '_centre_stim_Fraction_responsive_cells.tif'];
    centre_stim_fig_path = fullfile(save_fig_path_cell, centre_stim_fig_name);
    saveas(centre_frac_cells_fig, centre_stim_fig_path);

    % max ampl
    contra_max_ampl_fig = figure(Visible="off");
    plot(days_from_learning, all_contra_sharp_mean_max_amplitude, '-o')
    hold on;
    plot(days_from_learning, all_contra_wide_mean_max_amplitude, '-o')
    xlim([days_from_learning(1) days_from_learning(end)])
    xline(0, 'LineWidth', 2, 'Color', 'k')
    legend({'Sharp', 'Wide'})
    title([animal ' ' cell_group ' Contra Max amplitudes'])
    contra_stim_fig_name = [animal '_' cell_group '_contra_stim_Max_amplitudes.tif'];
    contra_stim_fig_path = fullfile(save_fig_path_cell, contra_stim_fig_name);
    saveas(contra_max_ampl_fig, contra_stim_fig_path);

    centre_max_ampl_fig = figure(Visible="off");
    plot(days_from_learning, all_centre_sharp_mean_max_amplitude, '-o')
    hold on;
    plot(days_from_learning, all_centre_wide_mean_max_amplitude, '-o')
    xlim([days_from_learning(1) days_from_learning(end)])
    xline(0, 'LineWidth', 2, 'Color', 'k')
    legend({'Sharp', 'Wide'})
    title([animal ' ' cell_group ' Centre Max amplitudes'])
    centre_stim_fig_name = [animal '_' cell_group '_centre_stim_Max_amplitudes.tif'];
    centre_stim_fig_path = fullfile(save_fig_path_cell, centre_stim_fig_name);
    saveas(centre_max_ampl_fig, centre_stim_fig_path);

    % psth across days
    psths_fig = figure; %(Visible="off");
    originalColormap = brewermap(length(psths_str_all_stim)+1, 'RdPu');
%     darkeningFactor = 0.1; 
%     modifiedColormap = originalColormap;
    modifiedColormap = originalColormap(2:end,:);
    axes('ColorOrder', modifiedColormap,'NextPlot','replacechildren');
    plot(bin_centres, vertcat(psths_str_all_stim{days_from_learning<0}), '--', 'linewidth', 2);
    hold on
    plot(bin_centres, vertcat(psths_str_all_stim{days_from_learning>=0}), 'linewidth', 2);
    xlim([-0.2 1])
    xlabel('Time from stim onset (s)')
    ylabel('Spike rate (spks/s)')
    title([animal ' ' cell_group ' PSTHs across days'])
    days_from_learning_str = arrayfun(@num2str, days_from_learning, 'UniformOutput', false);
    legend_for_plot = cellfun(@(x) ['Day ' x], days_from_learning_str, 'UniformOutput', false);
    legend(legend_for_plot);
    psths_fig_name = [animal '_' cell_group '_PSTHs_across_days.tif'];
    psths_fig_path = fullfile(save_fig_path_cell, psths_fig_name);
    saveas(psths_fig, psths_fig_path);

end
