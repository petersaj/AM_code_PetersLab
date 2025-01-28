animal = 'AM021';

load([animal '_swr_bhv_data.mat']);
load([animal '_task_str_ephys_data']);
load([animal '_str_ephys_data.mat']);
load([animal '_task_regression_data']);
load([animal '_per_neuron_task_regression_data']);


% passive
all_str_1_passive_max_amplitude = cell(size(passive_ephys_data, 1), 1);
all_str_2_passive_max_amplitude = cell(size(passive_ephys_data, 1), 1);
all_str_3_passive_max_amplitude = cell(size(passive_ephys_data, 1), 1);

% task
all_str_1_task_max_amplitude = cell(size(task_ephys_data, 1), 1);
all_str_2_task_max_amplitude = cell(size(task_ephys_data, 1), 1);
all_str_3_task_max_amplitude = cell(size(task_ephys_data, 1), 1);

all_days_from_learning = swr_bhv_data.days_from_learning{:};
days_from_learning = all_days_from_learning(3:end);


for this_rec=1:size(task_ephys_data, 1)

    rec_day = task_ephys_data.rec_day{this_rec};

    str_1_templates = task_ephys_data.str_1_templates{this_rec};
    str_2_templates = task_ephys_data.str_2_templates{this_rec};
    str_3_templates = task_ephys_data.str_3_templates{this_rec};
    str_all_templates = regression_per_neuron.str_all_templates{this_rec};

    [~, idx_str_1_templates, ~] = intersect(str_all_templates, str_1_templates, 'stable');
    [~, idx_str_2_templates, ~] = intersect(str_all_templates, str_2_templates, 'stable');
    [~, idx_str_3_templates, ~] = intersect(str_all_templates, str_3_templates, 'stable');

    % find resp neurons in passive
    passive_sharp_p_units = passive_ephys_data.contra_sharp_p_units{this_rec};
    passive_sharp_responsive_units = find(passive_sharp_p_units < 0.05);
    passive_sharp_unresponsive_units = find(passive_sharp_p_units >= 0.05);

    passive_wide_p_units = passive_ephys_data.contra_wide_p_units{this_rec};
    passive_wide_responsive_units = find(passive_wide_p_units < 0.05);
    passive_wide_unresponsive_units = find(passive_wide_p_units >= 0.05);

    % get for each part of striatum
    [str_1_passive_sharp_responsive_units{this_rec}, idx_str_1_passive_sharp_responsive_units{this_rec}, ~] = ...
        intersect(str_1_templates, passive_sharp_responsive_units, 'stable');
    [str_2_passive_sharp_responsive_units{this_rec}, idx_str_2_passive_sharp_responsive_units{this_rec}, ~] = ...
        intersect(str_2_templates, passive_sharp_responsive_units, 'stable');
    [str_3_passive_sharp_responsive_units{this_rec}, idx_str_3_passive_sharp_responsive_units{this_rec}, ~] = ...
        intersect(str_3_templates, passive_sharp_responsive_units, 'stable');
    [str_all_passive_sharp_responsive_units{this_rec}, idx_str_all_passive_sharp_responsive_units{this_rec}, ~] = ...
        intersect(str_all_templates, passive_sharp_responsive_units, 'stable');

    [str_1_passive_wide_responsive_units{this_rec}, idx_str_1_passive_wide_responsive_units{this_rec}, ~] = ...
        intersect(str_1_templates, passive_wide_responsive_units, 'stable');
    [str_2_passive_wide_responsive_units{this_rec}, idx_str_2_passive_wide_responsive_units{this_rec}, ~] = ...
        intersect(str_2_templates, passive_wide_responsive_units, 'stable');
    [str_3_passive_wide_responsive_units{this_rec}, idx_str_3_passive_wide_responsive_units{this_rec}, ~] = ...
        intersect(str_3_templates, passive_wide_responsive_units, 'stable');
    [str_all_passive_wide_responsive_units{this_rec}, idx_str_all_passive_wide_responsive_units{this_rec}, ~] = ...
        intersect(str_all_templates, passive_wide_responsive_units, 'stable');

    % find responsive neurons in task
    task_str_all_sign_rank_test_p_val = regression_per_neuron.str_all_sign_rank_test_p_val{this_rec};
    str_all_task_responsive_units{this_rec} = str_all_templates(task_str_all_sign_rank_test_p_val < 0.05);
    [str_1_task_responsive_units{this_rec}, idx_str_1_task_responsive_units{this_rec}, ~]  = ...
        intersect(str_1_templates, str_all_task_responsive_units{this_rec}, 'stable');
    [str_2_task_responsive_units{this_rec}, idx_str_2_task_responsive_units{this_rec}, ~]  = ...
        intersect(str_2_templates, str_all_task_responsive_units{this_rec}, 'stable');
    [str_3_task_responsive_units{this_rec}, idx_str_3_task_responsive_units{this_rec}, ~]  = ...
        intersect(str_3_templates, str_all_task_responsive_units{this_rec}, 'stable');


    %% AMPLITUDE
    %% - passive
    % LIKE IN PASSIVE CODE
%     passive_bin_centres = passive_ephys_data.bin_centres{this_rec};
%     passive_all_spikes_in_stim_time = passive_ephys_data.contra_all_spikes_in_stim_time{this_rec};
%     passive_mean_all_spikes_in_stim_time = squeeze(mean(passive_all_spikes_in_stim_time, 2));
% 
%     % smooth and normalize
%     gauss_win = gausswin(51, 3)';
%     passive_smooth_all_units = filter(gauss_win,sum(gauss_win),passive_mean_all_spikes_in_stim_time, [], 2);
%     passive_smooth_baseline = mean(passive_smooth_all_units(:,passive_bin_centres>-0.2 & passive_bin_centres<0), 2);
%     passive_norm_smooth_all_units = (passive_smooth_all_units - passive_smooth_baseline) ...
%         ./ (passive_smooth_baseline + std(passive_smooth_baseline));
% 
%     % get max amplitude for each striatum bit
%     str_1_passive_sharp_mean_max_amplitude = max(abs(mean(passive_norm_smooth_all_units(str_1_passive_sharp_responsive_units, :), 1)));
%     str_1_passive_wide_mean_max_amplitude = max(abs(mean(passive_norm_smooth_all_units(str_1_passive_wide_responsive_units, :), 1)));
%     all_str_1_passive_sharp_mean_max_amplitude(this_rec) = str_1_passive_sharp_mean_max_amplitude;
%     all_str_1_passive_wide_mean_max_amplitude(this_rec) = str_1_passive_wide_mean_max_amplitude;
% 
%     str_2_passive_sharp_mean_max_amplitude = max(abs(mean(passive_norm_smooth_all_units(str_2_passive_sharp_responsive_units, :), 1)));
%     str_2_passive_wide_mean_max_amplitude = max(abs(mean(passive_norm_smooth_all_units(str_2_passive_wide_responsive_units, :), 1)));
%     all_str_2_passive_sharp_mean_max_amplitude(this_rec) = str_2_passive_sharp_mean_max_amplitude;
%     all_str_2_passive_wide_mean_max_amplitude(this_rec) = str_2_passive_wide_mean_max_amplitude;
% 
%     str_3_passive_sharp_mean_max_amplitude = max(abs(mean(passive_norm_smooth_all_units(str_3_passive_sharp_responsive_units, :), 1)));
%     str_3_passive_wide_mean_max_amplitude = max(abs(mean(passive_norm_smooth_all_units(str_3_passive_wide_responsive_units, :), 1)));
%     all_str_3_passive_sharp_mean_max_amplitude(this_rec) = str_3_passive_sharp_mean_max_amplitude;
%     all_str_3_passive_wide_mean_max_amplitude(this_rec) = str_3_passive_wide_mean_max_amplitude;

    % MORE SIMILAR TO TASK
    passive_rec_day = find(ismember(passive_ephys_data.rec_day, rec_day));
    passive_stimOn_times = passive_ephys_data.stimOn_times{passive_rec_day};
    passive_spike_times_timelite = passive_ephys_data.spike_times_timeline{passive_rec_day};
    passive_spike_templates = passive_ephys_data.spike_templates{passive_rec_day};
    passive_good_trials = passive_ephys_data.contra_good_trials{passive_rec_day};
    good_templates_idx = passive_ephys_data.good_templates{passive_rec_day};

    % create time vector around all stim onsets
    bin_window = 0.01;
    bin_edges = -0.5:bin_window:2;
    passive_around_stim_time = passive_stimOn_times(passive_good_trials) + bin_edges;

    % for plot
    bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;

    % task stim
    passive_spikes_in_stim_time_downsampled = nan(length(good_templates_idx), length(passive_stimOn_times(passive_good_trials)), size(passive_around_stim_time, 2)-1);
    for unit_idx=1:length(good_templates_idx)
        unit_spikes = passive_spike_templates == unit_idx;
        unit_spike_times = passive_spike_times_timelite(unit_spikes);
        passive_spikes_in_stim_time_downsampled(unit_idx, :, :) = cell2mat(arrayfun(@(trial_id) histcounts(unit_spike_times, passive_around_stim_time(trial_id,:))', ...
            1:size(passive_around_stim_time, 1), 'UniformOutput',false))' / bin_window;
    end

    % get avg across trials
    passive_psth_stim_str_1_all = squeeze(mean(passive_spikes_in_stim_time_downsampled(str_1_templates, :, :), 2));
    passive_psth_stim_str_2_all = squeeze(mean(passive_spikes_in_stim_time_downsampled(str_2_templates, :, :), 2));
    passive_psth_stim_str_3_all = squeeze(mean(passive_spikes_in_stim_time_downsampled(str_3_templates, :, :), 2));

    % baseline
    passive_str_1_stim_baseline = mean(passive_psth_stim_str_1_all(bin_edges<0), 2);
    passive_str_2_stim_baseline = mean(passive_psth_stim_str_2_all(bin_edges<0), 2);
    passive_str_3_stim_baseline = mean(passive_psth_stim_str_3_all(bin_edges<0), 2);

    % baseline subtract
    passive_psth_stim_str_1_all_baseline_sub = passive_psth_stim_str_1_all - repmat(passive_str_1_stim_baseline,1,size(passive_psth_stim_str_1_all,2));
    %     passive_psth_stim_str_1_mean_baseline_sub = mean(passive_psth_stim_str_1_all_baseline_sub, 1)';
    passive_psth_stim_str_2_all_baseline_sub = passive_psth_stim_str_2_all - repmat(passive_str_2_stim_baseline,1,size(passive_psth_stim_str_2_all,2));
    %     passive_psth_stim_str_2_mean_baseline_sub = mean(passive_psth_stim_str_2_all_baseline_sub, 1)';
    passive_psth_stim_str_3_all_baseline_sub = passive_psth_stim_str_3_all - repmat(passive_str_3_stim_baseline,1,size(passive_psth_stim_str_3_all,2));
    %     passive_psth_stim_str_3_mean_baseline_sub = mean(passive_psth_stim_str_3_all_baseline_sub, 1)';

    % get max ampl
    all_str_1_passive_max_amplitude{this_rec} = max(passive_psth_stim_str_1_all_baseline_sub, [], 2);
    all_str_2_passive_max_amplitude{this_rec} = max(passive_psth_stim_str_2_all_baseline_sub, [], 2);
    all_str_3_passive_max_amplitude{this_rec} = max(passive_psth_stim_str_3_all_baseline_sub, [], 2);

    %% - task
    task_mean_neuron_str_all_trace_residual_stim_align = regression_per_neuron.mean_neuron_str_all_trace_residual_stim_align{this_rec};
    all_str_all_task_max_amplitude = cellfun(@(neuron_residual_stim_align) max(neuron_residual_stim_align), ...
        task_mean_neuron_str_all_trace_residual_stim_align, 'UniformOutput', true);
    all_str_1_task_max_amplitude{this_rec} = all_str_all_task_max_amplitude(idx_str_1_templates);
    all_str_2_task_max_amplitude{this_rec} = all_str_all_task_max_amplitude(idx_str_2_templates);
    all_str_3_task_max_amplitude{this_rec} = all_str_all_task_max_amplitude(idx_str_3_templates);

    %     temp_all_str_1_task_max_amplitude = all_str_1_task_max_amplitude{this_rec};
    %     temp_all_str_2_task_max_amplitude = all_str_2_task_max_amplitude{this_rec};
    %     temp_all_str_3_task_max_amplitude = all_str_3_task_max_amplitude{this_rec};

    %     %% get responsive amplitudes
    %     responsive_str_1_task_max_amplitude = temp_all_str_1_task_max_amplitude(idx_str_1_task_responsive_units{this_rec});
    %     responsive_str_1_task_max_amplitude = temp_all_str_1_task_max_amplitude(idx_str_1_task_responsive_units{this_rec});
    %     responsive_str_1_task_max_amplitude = temp_all_str_1_task_max_amplitude(idx_str_1_task_responsive_units{this_rec});
    %
    %% frac of cells
    str_1_task_frac_cells(this_rec) = length(str_1_task_responsive_units{this_rec})/length(str_1_templates);
    str_2_task_frac_cells(this_rec) = length(str_2_task_responsive_units{this_rec})/length(str_2_templates);
    str_3_task_frac_cells(this_rec) = length(str_3_task_responsive_units{this_rec})/length(str_3_templates);

    str_1_passive_sharp_frac_cells(this_rec) = length(str_1_passive_sharp_responsive_units{this_rec})/length(str_1_templates);
    str_2_passive_sharp_frac_cells(this_rec) = length(str_2_passive_sharp_responsive_units{this_rec})/length(str_2_templates);
    str_3_passive_sharp_frac_cells(this_rec) = length(str_3_passive_sharp_responsive_units{this_rec})/length(str_3_templates);

    str_1_passive_wide_frac_cells(this_rec) = length(str_1_passive_wide_responsive_units{this_rec})/length(str_1_templates);
    str_2_passive_wide_frac_cells(this_rec) = length(str_2_passive_wide_responsive_units{this_rec})/length(str_2_templates);
    str_3_passive_wide_frac_cells(this_rec) = length(str_3_passive_wide_responsive_units{this_rec})/length(str_3_templates);

end

% max amplitudes
figure;
tiledlayout("flow")
for this_rec=1:size(task_ephys_data, 1)
    nexttile
    learning_day = days_from_learning(this_rec);
    plot(all_str_1_passive_max_amplitude{this_rec}, all_str_1_task_max_amplitude{this_rec}', '.')
    title(['Day ' num2str(learning_day)])
    xlabel('Passive')
    ylabel('Task')
end
sgtitle('Striatum 1')

figure;
tiledlayout("flow")
for this_rec=1:size(task_ephys_data, 1)
    nexttile
    learning_day = days_from_learning(this_rec);
    plot(all_str_2_passive_max_amplitude{this_rec}, all_str_2_task_max_amplitude{this_rec}', '.')
    title(['Day ' num2str(learning_day)])
    xlabel('Passive')
    ylabel('Task')
end
sgtitle('Striatum 2')

figure;
tiledlayout("flow")
for this_rec=1:size(task_ephys_data, 1)
    nexttile
    learning_day = days_from_learning(this_rec);
    plot(all_str_3_passive_max_amplitude{this_rec}, all_str_3_task_max_amplitude{this_rec}', '.')
    title(['Day ' num2str(learning_day)])
    xlabel('Passive')
    ylabel('Task')
end
sgtitle('Striatum 3')

% frac across days
figure;
tiledlayout(3,1);
nexttile;
plot(days_from_learning, str_1_task_frac_cells)
hold on;
plot(days_from_learning, str_1_passive_sharp_frac_cells)
hold on;
plot(days_from_learning, str_1_passive_wide_frac_cells)
legend({'Task', 'Passive Sharp', 'Passive Wide'})
title('Striatum 1')

nexttile;
plot(days_from_learning, str_2_task_frac_cells)
hold on;
plot(days_from_learning, str_2_passive_sharp_frac_cells)
hold on;
plot(days_from_learning, str_2_passive_wide_frac_cells)
legend({'Task', 'Passive Sharp', 'Passive Wide'})
title('Striatum 2')

nexttile;
plot(days_from_learning, str_3_task_frac_cells)
hold on;
plot(days_from_learning, str_3_passive_sharp_frac_cells)
hold on;
plot(days_from_learning, str_3_passive_wide_frac_cells)
legend({'Task', 'Passive Sharp', 'Passive Wide'})
title('Striatum 3')

% compare responsive cells
str_1_same_responsive_cells = cell(1,size(task_ephys_data, 1));
str_1_task_not_passive_cells = cell(1,size(task_ephys_data, 1));
str_1_passive_not_task_cells = cell(1,size(task_ephys_data, 1));
for this_rec=1:size(task_ephys_data, 1)
    str_1_same_responsive_cells{this_rec} = intersect(str_1_task_responsive_units{this_rec}, str_1_passive_sharp_responsive_units{this_rec});
    str_1_task_not_passive_cells{this_rec} = setdiff(str_1_task_responsive_units{this_rec}, str_1_passive_sharp_responsive_units{this_rec});
    str_1_passive_not_task_cells{this_rec} = setdiff(str_1_passive_sharp_responsive_units{this_rec}, str_1_task_responsive_units{this_rec});

    disp(['Str 1 Day ' num2str(days_from_learning(this_rec))])
    disp(['Same ' num2str(length(str_1_same_responsive_cells{this_rec})) ' units:'])
    disp(str_1_same_responsive_cells{this_rec})
    disp(['Task not passive ' num2str(length(str_1_task_not_passive_cells{this_rec})) ' units:'])
    disp(str_1_task_not_passive_cells{this_rec})
    disp(['Passive not task ' num2str(length(str_1_passive_not_task_cells{this_rec})) ' units:'])
    disp(str_1_passive_not_task_cells{this_rec})
end


figure;
tiledlayout('flow');
for this_rec=1:size(task_ephys_data, 1)
    nexttile;
    temp_max_ampl = all_str_1_passive_max_amplitude{this_rec};
    plot(temp_max_ampl, task_ephys_data.str_1_templates{this_rec}, '.')
    hold on;
    plot(temp_max_ampl(str_1_same_responsive_cells{this_rec}), str_1_same_responsive_cells{this_rec}, 'o')
    hold on;
    plot(temp_max_ampl(), str_1_task_not_passive_cells{this_rec}, '*')
end

%% RANDOM
this_rec=size(task_ephys_data, 1);
str_all_templates = regression_per_neuron.str_all_templates{this_rec};
str_all_sign_rank_test_p_val = regression_per_neuron.str_all_sign_rank_test_p_val{this_rec};
mean_neuron_str_all_trace_residual_stim_align = regression_per_neuron.mean_neuron_str_all_trace_residual_stim_align{this_rec};
mean_neuron_stim_align_str_all_trace_baseline_sub = regression_per_neuron.mean_neuron_stim_align_str_all_trace_baseline_sub{this_rec};
bin_window_for_test = regression_per_neuron.bin_window_for_test;

figure;
tiledlayout('flow');
for neuron_idx = 1:length(str_all_templates)
    if str_all_sign_rank_test_p_val(neuron_idx)<0.05
        nexttile
        plot(bin_edges, mean_neuron_str_all_trace_residual_stim_align{neuron_idx})
        xline([-bin_window_for_test bin_window_for_test])
        title(['Neuron ' num2str(neuron_idx) ' test val ' num2str(str_all_sign_rank_test_p_val(neuron_idx))])
    end
end
sgtitle(['Residual trace Responsive bin ' num2str(bin_window_for_test)])

figure;
tiledlayout('flow');
for neuron_idx = 1:length(str_all_templates)
    if str_all_sign_rank_test_p_val(neuron_idx)<0.05
        nexttile
        plot(bin_edges, mean_neuron_stim_align_str_all_trace_baseline_sub{neuron_idx})
        xline([-bin_window_for_test bin_window_for_test])
        title(['Neuron ' num2str(neuron_idx) ' test val ' num2str(str_all_sign_rank_test_p_val(neuron_idx))])
    end
end
sgtitle('Full trace Responsive')


% make e.g. plot for neuron
neuron_idx = 20;

figure;
tiledlayout('flow');
res_plot = nexttile;
plot(bin_edges, mean_neuron_str_all_trace_residual_stim_align{neuron_idx})
xline([-bin_window_for_test bin_window_for_test])
title('Residual')
ylim([-11 40])
xlim([-0.5 1])
res_plot.ColorOrder =  brewermap(5, 'PiYg');

full_plot = nexttile;
plot(bin_edges, mean_neuron_stim_align_str_all_trace_baseline_sub{neuron_idx})
xline([-bin_window_for_test bin_window_for_test])
title('Full')
ylim([-11 40])
xlim([-0.5 1])
full_plot.ColorOrder =  brewermap(5, 'PiYg');




        %     figure;
%     tiledlayout('flow');
%     for neuron_idx = 1:10%length(str_1_templates)
%         if str_1_sign_rank_test_p_val(neuron_idx)>=0.05
%             nexttile
%             plot(bin_edges, mean_neuron_stim_align_str_1_trace_baseline_sub{neuron_idx})
%             xline([-bin_window_for_test bin_window_for_test])
%             title(['Neuron ' num2str(neuron_idx) ' test val ' num2str(str_1_sign_rank_test_p_val(neuron_idx))])
%         end
%     end
%     sgtitle('Full trace Non-Responsive')

