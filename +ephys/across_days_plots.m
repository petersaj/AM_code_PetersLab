function across_days_plots(animal, ephys_data, save_fig_path, cell_group)

save_fig_path_cell = fullfile(save_fig_path, cell_group);

if ~isfolder(save_fig_path_cell)
    mkdir(save_fig_path_cell)
end

% initialize arrays
num_cells_contra_sharp_responsive_units = nan(height(ephys_data), 1);
num_cells_contra_wide_responsive_units = nan(height(ephys_data), 1);
num_cells_centre_sharp_responsive_units = nan(height(ephys_data), 1);
num_cells_centre_wide_responsive_units = nan(height(ephys_data), 1);

all_contra_sharp_mean_max_amplitude = nan(height(ephys_data), 1);
all_contra_wide_mean_max_amplitude = nan(height(ephys_data), 1);
all_centre_sharp_mean_max_amplitude = nan(height(ephys_data), 1);
all_centre_wide_mean_max_amplitude = nan(height(ephys_data), 1);

for use_rec=1:height(ephys_data)

    good_templates = ephys_data.good_templates{use_rec};
    good_templates_idx = find(good_templates);

    if strcmpi(cell_group,'all')
        units_cell_group = ones(length(good_templates_idx), 1);
    else
        units_cell_group = ephys_data.(lower(cell_group)){use_rec};
    end

    contra_sharp_p_units = ephys_data.contra_sharp_p_units{use_rec};
    contra_wide_p_units = ephys_data.contra_wide_p_units{use_rec};

    centre_sharp_p_units = ephys_data.centre_sharp_p_units{use_rec};
    centre_wide_p_units = ephys_data.centre_wide_p_units{use_rec};

    %% get responsive units in this cell group
    % contra
    contra_sharp_responsive_units = find(units_cell_group & contra_sharp_p_units < 0.05);
    contra_wide_responsive_units = find(units_cell_group & contra_wide_p_units < 0.05);

    % centre
    centre_sharp_responsive_units = find(units_cell_group & centre_sharp_p_units < 0.05);
    centre_wide_responsive_units = find(units_cell_group & centre_wide_p_units < 0.05);

 

end

%% plots
% num cells
contra_num_cells_fig = figure;
plot(num_cells_contra_sharp_responsive_units, '-o')
hold on;
plot(num_cells_contra_wide_responsive_units, '-o')
legend({'Sharp', 'Wide'})
title('Contra stim. Number of responsive cells across days')
contra_stim_fig_name = [animal '_' cell_group '_contra_stim_Number_responsive_cells.tif'];
contra_stim_fig_path = fullfile(save_fig_path_cell, contra_stim_fig_name);
saveas(contra_num_cells_fig, contra_stim_fig_path);

centre_num_cells_fig = figure;
plot(num_cells_centre_sharp_responsive_units, '-o')
hold on;
plot(num_cells_centre_wide_responsive_units, '-o')
legend({'Sharp', 'Wide'})
title('Centre stim. Number of responsive cells across days')
centre_stim_fig_name = [animal '_' cell_group '_centre_stim_Number_responsive_cells.tif'];
centre_stim_fig_path = fullfile(save_fig_path_cell, centre_stim_fig_name);
saveas(centre_num_cells_fig, centre_stim_fig_path);

% max ampl
contra_max_ampl_fig = figure;
plot(all_contra_sharp_mean_max_amplitude, '-o')
hold on;
plot(all_contra_wide_mean_max_amplitude, '-o')
legend({'Sharp', 'Wide'})
title('Contra Max amplitudes')
contra_stim_fig_name = [animal '_' cell_group '_contra_stim_Max_amplitudes.tif'];
contra_stim_fig_path = fullfile(save_fig_path_cell, contra_stim_fig_name);
saveas(contra_max_ampl_fig, contra_stim_fig_path);

centre_max_ampl_fig = figure;
plot(all_centre_sharp_mean_max_amplitude, '-o')
hold on;
plot(all_centre_wide_mean_max_amplitude, '-o')
legend({'Sharp', 'Wide'})
title('Centre Max amplitudes')
centre_stim_fig_name = [animal '_' cell_group '_centre_stim_Max_amplitudes.tif'];
centre_stim_fig_path = fullfile(save_fig_path_cell, centre_stim_fig_name);
saveas(centre_max_ampl_fig, centre_stim_fig_path);