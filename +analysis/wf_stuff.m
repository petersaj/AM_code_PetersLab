%% SAVE
%% - Load dataset
load('D:\matlab_save\all_swr_bhv_data.mat');

animals = ["AM011"; "AM012"; "AM014"; "AM015"; "AM016"; ...
    "AM017"; "AM018"; "AM019"; "AM021"; "AM022"; "AM026"]; %; "AM029"];

%% - make roi for vis
animal = 'AM021';
workflow = {'lcr_passive'};
recordings = plab.find_recordings(animal, [], workflow);
use_rec = length(recordings)-1;
rec_day = recordings(use_rec).day;
rec_time = recordings(use_rec).recording{end};
verbose = true;
load_parts.behavior = true;
load_parts.widefield = true;
ap.load_recording

% vis
figure; colormap(gray);
imagesc(wf_avg);
axis image off
clim([0, 20000])
ap.wf_draw('ccf', 'y')
vis_roi_poly = drawpolygon;
vis_roi_mask = createMask(vis_roi_poly);

% pfc
figure; colormap(gray);
imagesc(wf_avg);
axis image off
clim([0, 20000])
ap.wf_draw('ccf', 'y')
pfc_roi_poly = drawpolygon;
pfc_roi_mask = createMask(pfc_roi_poly);
% save('D:\matlab_save\new_pfc_ROI', "new_pfc_roi_mask", "-v7.3");

figure; imagesc(vis_roi_mask); axis image;
title('vis ROI mask');

figure; imagesc(pfc_roi_mask); axis image;
title('pfc ROI mask');

save('D:\matlab_save\ROIs', "vis_roi_mask", "pfc_roi_mask", "-v7.3");

%% load ROIs
load('D:\matlab_save\ROIs');

master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');

%% - go through each animal
for animal_idx=1:length(animals)

    animal = animals{animal_idx};
    disp(['Start ' animal])
    workflow_passive = {'lcr_passive'};
    recordings_passive = plab.find_recordings(animal, [], workflow_passive);
    workflow_task = {'stim_wheel_right*'};
    recordings_task = plab.find_recordings(animal, [], workflow_task);
    training_days = ismember({recordings_passive.day}, {recordings_task.day});
    train_rec_passive = recordings_passive(training_days);

    bhv_days = swr_bhv_data.bhv_days{animal_idx};
    days_from_learning = swr_bhv_data.days_from_learning{animal_idx};
    
    wf_save = table;

    for use_rec=1:length(train_rec_passive)
        rec_day = train_rec_passive(use_rec).day;
        rec_time = train_rec_passive(use_rec).recording{end};
        verbose = true;
        load_parts.behavior = true;
        load_parts.widefield = true;
        ap.load_recording

        this_day_from_learning = days_from_learning(ismember(bhv_days, rec_day));

        wf_save.rec_day(use_rec) = {rec_day};
        wf_save.days_from_learning(use_rec) = {this_day_from_learning};

        %% -- get contra stuff
        contra_stim = 90;
        trial_stim_values = vertcat(trial_events.values.TrialStimX);
        trial_stim_values = trial_stim_values(1:length(stimOn_times));

        % stim aligned wheel move
        % create matrix of times for stim onset
        timestep = 0.01;
        start_time = -0.5;
        end_time = 1;  
        timevec = start_time:timestep:end_time;
        stim_frame = (-start_time)*(1/timestep)+1;
        time_stimulus = stimOn_times+timevec;
        % stim aligned wheel move
        t = timelite.timestamps;
        stim_wheel_move = interp1(t,+wheel_move,time_stimulus);
        no_move_trials = sum(stim_wheel_move(:,stim_frame:end),2)==0;

        contra_good_trials = (trial_stim_values == contra_stim) & no_move_trials;

        %% -- align to master
        [U_master,V_master] = plab.wf.u2master(wf_U,wf_V);

        %% -- interp wf
%         V_interp = interp1(wf_t,wf_V',grab_time)';
        
        wf_stim_time = -0.2:0.01:0.7;
        wf_stim_interp_time = stimOn_times + wf_stim_time;

        V_stim_align = interp1(wf_t,V_master',wf_stim_interp_time);
        V_avg_contra_stim_align = squeeze(mean(V_stim_align(contra_good_trials, :, :), 1))';

%         px_avg_contra_stim_align = plab.wf.svd2px(U_master,V_avg_contra_stim_align);
%         ap.imscroll(px_avg_contra_stim_align);
%         axis image
%         colormap(ap.colormap('PWG'));
%         title('Example map');
%         clim_val = max(abs(clim)) * 0.7;
%         clim([-clim_val, clim_val]); % Set the color limits
% %         clim(clim/8)
        
        %% -- get ROIs
        roi_stim_time = -0.5:0.01:1;
        roi_stim_interp_time = stimOn_times + roi_stim_time;

        vis_roi_trace = ap.wf_roi(U_master,V_master,wf_avg,[],vis_roi_mask);
        vis_roi_stim_align = interp1(wf_t,vis_roi_trace,roi_stim_interp_time)';
        vis_avg_contra_roi_stim_align = mean(vis_roi_stim_align(:, contra_good_trials), 2);
%         figure;
%         plot(vis_avg_contra_roi_stim_align)
      
        pfc_roi_trace = ap.wf_roi(U_master,V_master,wf_avg,[],pfc_roi_mask);
        pfc_roi_stim_align = interp1(wf_t,pfc_roi_trace,roi_stim_interp_time)';
        pfc_avg_contra_roi_stim_align = mean(pfc_roi_stim_align(:, contra_good_trials), 2);

        %% -- save

        wf_save.wf_stim_time(use_rec) = {wf_stim_time};
        wf_save.V_avg_contra_stim_align(use_rec) = {V_avg_contra_stim_align};

        wf_save.roi_stim_time(use_rec) = {roi_stim_time};
        wf_save.vis_avg_contra_roi_stim_align(use_rec) = {vis_avg_contra_roi_stim_align};
        wf_save.pfc_avg_contra_roi_stim_align(use_rec) = {pfc_avg_contra_roi_stim_align};

        disp(['Done day ' num2str(use_rec)])
        
    end
    all_wf_save_cell{animal_idx} = wf_save;
    disp(['Done ' animal])
end

ctx_wf = table();
ctx_wf.animals = animals;
ctx_wf.recording_data = all_wf_save_cell';

ctx_wf.U_master(:) = {U_master};

% ctx_wf.vis_roi_mask(:) = {vis_roi_mask};
% ctx_wf.pfc_roi_mask(:) = {pfc_roi_mask};


save_name = 'D:\matlab_save\ctx_wf';
save(save_name, "ctx_wf", "-v7.3");


%% LOAD and AVG 
load('D:\matlab_save\ctx_wf');

animals = ctx_wf.animals;

%% - visual
% Initialize a structure to hold grouped data
all_vis_avg_grouped = struct();

% Iterate over all animals
for animal_idx = 1:numel(ctx_wf.animals)
    % Get recording data for the current animal
    recording_data = ctx_wf.recording_data{animal_idx};
    
    % Extract the 'days from learning' column
    days_from_learning = recording_data.days_from_learning;
    
    % Extract the corresponding 'vis_avg_contra_roi_stim_align' vectors
    vis_avg_contra_roi_stim_align = recording_data.vis_avg_contra_roi_stim_align;
    
    % Iterate over each day and its corresponding vector
    for day_idx = 1:numel(days_from_learning)
        this_day = days_from_learning{day_idx};
        day_vis_avg_contra_roi_stim_align = vis_avg_contra_roi_stim_align{day_idx}; % assuming vectors are stored in a cell array
        
        % Create a valid field name for the day
        day_field = ['day_', strrep(num2str(this_day), '-', 'neg')];
        
        % Check if the day already exists in the grouped data
        if ~isfield(all_vis_avg_grouped, day_field)
            all_vis_avg_grouped.(day_field) = {}; % Initialize as a cell array
        end
        
        % Append the vector to the corresponding day
        all_vis_avg_grouped.(day_field){end+1} = day_vis_avg_contra_roi_stim_align;
    end
end

% get roi stim time
roi_stim_time = recording_data.roi_stim_time{1};

% combine and avg for each day
day_fields = fieldnames(all_vis_avg_grouped);
all_vis_avg_grouped_combined = struct();
for day_idx = 1:numel(day_fields)
    day_field = day_fields{day_idx};
    all_vis_avg_grouped_combined.(day_field) = horzcat(all_vis_avg_grouped.(day_field){:}); % Concatenate vectors for this day
end


% Get all day fields from grouped_combined
day_fields = fieldnames(grouped_combined);

% Initialize storage for days, averages, and SEMs
sorted_days = zeros(numel(day_fields), 1);
mean_vis_avg_grouped = {};
sems_vis_avg_grouped = {};

% Process each day
for i = 1:numel(day_fields)
    day_field = day_fields{i};
    
    % Extract the numerical day value (handle 'neg' in field names)
    sorted_days(i) = str2double(strrep(strrep(day_field, 'day_neg', '-'), 'day_', ''));
    
    % Extract data for the current day (already a matrix)
    day_matrix = all_vis_avg_grouped_combined.(day_field); % Each column is a variable
    
    % Compute average and SEM across columns
    avg = mean(day_matrix, 2);             
    sem = std(day_matrix, 0, 2) / sqrt(size(day_matrix, 2)); 
    
    % Store results
    mean_vis_avg_grouped{i} = avg;
    sems_vis_avg_grouped{i} = sem;
end

% Sort results by day
[sorted_days, sort_idx] = sort(sorted_days);
mean_vis_avg_grouped = mean_vis_avg_grouped(sort_idx);
sems_vis_avg_grouped = sems_vis_avg_grouped(sort_idx);

%% - pfc
% Initialize a structure to hold grouped data
all_pfc_avg_grouped = struct();

% Iterate over all animals
for animal_idx = 1:numel(ctx_wf.animals)
    % Get recording data for the current animal
    recording_data = ctx_wf.recording_data{animal_idx};
    
    % Extract the 'days from learning' column
    days_from_learning = recording_data.days_from_learning;
    
    % Extract the corresponding 'pfc_avg_contra_roi_stim_align' vectors
    pfc_avg_contra_roi_stim_align = recording_data.pfc_avg_contra_roi_stim_align;
    
    % Iterate over each day and its corresponding vector
    for day_idx = 1:numel(days_from_learning)
        this_day = days_from_learning{day_idx};
        day_pfc_avg_contra_roi_stim_align = pfc_avg_contra_roi_stim_align{day_idx}; % assuming vectors are stored in a cell array
        
        % Create a valid field name for the day
        day_field = ['day_', strrep(num2str(this_day), '-', 'neg')];
        
        % Check if the day already exists in the grouped data
        if ~isfield(all_pfc_avg_grouped, day_field)
            all_pfc_avg_grouped.(day_field) = {}; % Initialize as a cell array
        end
        
        % Append the vector to the corresponding day
        all_pfc_avg_grouped.(day_field){end+1} = day_pfc_avg_contra_roi_stim_align;
    end
end

% get roi stim time
roi_stim_time = recording_data.roi_stim_time{1};

% combine and avg for each day
day_fields = fieldnames(all_pfc_avg_grouped);
all_pfc_avg_grouped_combined = struct();
for day_idx = 1:numel(day_fields)
    day_field = day_fields{day_idx};
    all_pfc_avg_grouped_combined.(day_field) = horzcat(all_pfc_avg_grouped.(day_field){:}); % Concatenate vectors for this day
end


% Get all day fields from grouped_combined
day_fields = fieldnames(all_pfc_avg_grouped_combined);

% Initialize storage for days, averages, and SEMs
sorted_days = zeros(numel(day_fields), 1);
mean_pfc_avg_grouped = {};
sems_pfc_avg_grouped = {};

% Process each day
for i = 1:numel(day_fields)
    day_field = day_fields{i};
    
    % Extract the numerical day value (handle 'neg' in field names)
    sorted_days(i) = str2double(strrep(strrep(day_field, 'day_neg', '-'), 'day_', ''));
    
    % Extract data for the current day (already a matrix)
    day_matrix = all_pfc_avg_grouped_combined.(day_field); % Each column is a variable
    
    % Compute average and SEM across columns
    avg = mean(day_matrix, 2);             
    sem = std(day_matrix, 0, 2) / sqrt(size(day_matrix, 2)); 
    
    % Store results
    mean_pfc_avg_grouped{i} = avg;
    sems_pfc_avg_grouped{i} = sem;
end

% Sort results by day
[sorted_days, sort_idx] = sort(sorted_days);
mean_pfc_avg_grouped = mean_pfc_avg_grouped(sort_idx);
sems_pfc_avg_grouped = sems_pfc_avg_grouped(sort_idx);


%% - plot
% visual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of days to plot
num_days = numel(sorted_days);

figure;
tiledlayout(4,3);  

% Loop through each day and plot in a separate subplot
for i = 1:num_days
    
    if sorted_days(i) < 0
        day_field = sprintf('day_neg%d', abs(sorted_days(i)));  % For negative days, use day_neg3, day_neg2, etc.
    else
        day_field = sprintf('day_%d', sorted_days(i));  % For non-negative days, use day_3, day_2, etc.
    end

    % Data for the current day
    x = roi_stim_time;
    avg = mean_vis_avg_grouped{i}';        % Average (1x151)
    sem = sems_vis_avg_grouped{i}';            % SEM (1x151)

    % Get the size of the data for the current day (number of repetitions)
    num_reps = size(grouped_combined.(day_field), 2);  % Number of rows (repetitions)

    % Create a new tile for each day
    nexttile;
    hold on;
    
    % Plot average as a line
    plot(x, avg, 'LineWidth', 1.5, 'Color', 'k');
    
    % Create the x and y values for the SEM region
    x_fill = [x, fliplr(x)];  % X values: time points and reversed time points
    y_fill = [avg + sem, fliplr(avg - sem)];  % Y values: upper and lower bounds of SEM
    
    % Plot the shaded SEM region
    fill(x_fill, y_fill, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % Title and labels for each subplot
    title_str = sprintf('Day %d (Repetitions: %d)', sorted_days(i), num_reps);  % Add repetitions to the title
    title(title_str);
    xlabel('Time from stim onset');
    ylabel('\Delta F/F', 'FontSize', 12);
    hold off;
end

% pfc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of days to plot
num_days = numel(sorted_days);

figure;
tiledlayout(4,3);  

% Loop through each day and plot in a separate subplot
for i = 1:num_days

    if sorted_days(i) < 0
        day_field = sprintf('day_neg%d', abs(sorted_days(i)));  % For negative days, use day_neg3, day_neg2, etc.
    else
        day_field = sprintf('day_%d', sorted_days(i));  % For non-negative days, use day_3, day_2, etc.
    end

    % Data for the current day
    x = roi_stim_time;
    avg = mean_pfc_avg_grouped{i}';        % Average (1x151)
    sem = sems_pfc_avg_grouped{i}';            % SEM (1x151)
    
    % Get the size of the data for the current day (number of repetitions)
    num_reps = size(all_pfc_avg_grouped_combined.(day_field), 2);  % Number of rows (repetitions)

    % Create a new tile for each day
    nexttile;
    hold on;
    
    % Create the x and y values for the SEM region
    x_fill = [x, fliplr(x)];  % X values: time points and reversed time points
    y_fill = [avg + sem, fliplr(avg - sem)];  % Y values: upper and lower bounds of SEM
    
    % Plot the shaded SEM region
    fill(x_fill, y_fill, [0 0.5 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % Plot average as a line
    plot(x, avg, 'LineWidth', 1.5, 'Color', [0 0.5 0]);

    % Title and labels for each subplot
    title_str = sprintf('Day %d (Repetitions: %d)', sorted_days(i), num_reps);  % Add repetitions to the title
    title(title_str);
    xlabel('Time from stim onset');
    ylabel('\Delta F/F', 'FontSize', 12);
    xline(0, 'LineWidth', 2)
    ylim([-10^(-3) 2*10^(-3)])
    hold off;
    xlim([-0.1 0.3])
    AP_scalebar(0.1, 10^(-3))
    ap.prettyfig
end


%% Get new pfc ROI from avg Vs

load('D:\matlab_save\ctx_wf');
animals = ctx_wf.animals;

load("D:\matlab_save\new_pfc_ROI.mat");
figure; imagesc(new_pfc_roi_mask); axis image;
title('pfc ROI mask');
% Initialize a structure to hold grouped data
all_new_pfc_avg_grouped = struct();

% Iterate over all animals
for animal_idx = 1:numel(ctx_wf.animals)
    % Get recording data for the current animal
    recording_data = ctx_wf.recording_data{animal_idx};
    
    % Extract the 'days from learning' column
    days_from_learning = recording_data.days_from_learning;
    
    % get U
    U_master = ctx_wf.U_master{animal_idx};

    % Iterate over each day and its corresponding vector
    for day_idx = 1:numel(days_from_learning)
        this_day = days_from_learning{day_idx};

        % get the Vs
        V_avg_contra_stim_align = recording_data.V_avg_contra_stim_align{day_idx};

        % get new ROI
        day_new_pfc_avg_contra_roi_stim_align = ap.wf_roi(U_master,V_avg_contra_stim_align,[],[],new_pfc_roi_mask);
        
        % Create a valid field name for the day
        day_field = ['day_', strrep(num2str(this_day), '-', 'neg')];
        
        % Check if the day already exists in the grouped data
        if ~isfield(all_new_pfc_avg_grouped, day_field)
            all_new_pfc_avg_grouped.(day_field) = {}; % Initialize as a cell array
        end
        
        % Append the vector to the corresponding day
        all_new_pfc_avg_grouped.(day_field){end+1} = day_new_pfc_avg_contra_roi_stim_align';
    end
end

% !!!!!!!!! get wf stim time
wf_stim_time = recording_data.wf_stim_time{1};

% combine and avg for each day
day_fields = fieldnames(all_new_pfc_avg_grouped);
all_new_pfc_avg_grouped_combined = struct();
for day_idx = 1:numel(day_fields)
    day_field = day_fields{day_idx};
    all_new_pfc_avg_grouped_combined.(day_field) = horzcat(all_new_pfc_avg_grouped.(day_field){:}); % Concatenate vectors for this day
end


% Get all day fields from grouped_combined
day_fields = fieldnames(all_new_pfc_avg_grouped_combined);

% Initialize storage for days, averages, and SEMs
sorted_days = zeros(numel(day_fields), 1);
mean_new_pfc_avg_grouped = {};
sems_new_pfc_avg_grouped = {};

% Process each day
for i = 1:numel(day_fields)
    day_field = day_fields{i};
    
    % Extract the numerical day value (handle 'neg' in field names)
    sorted_days(i) = str2double(strrep(strrep(day_field, 'day_neg', '-'), 'day_', ''));
    
    % Extract data for the current day (already a matrix)
    day_matrix = all_new_pfc_avg_grouped_combined.(day_field); % Each column is a variable
    
    % Compute average and SEM across columns
    avg = mean(day_matrix, 2);             
    sem = std(day_matrix, 0, 2) / sqrt(size(day_matrix, 2)); 
    
    % Store results
    mean_new_pfc_avg_grouped{i} = avg;
    sems_new_pfc_avg_grouped{i} = sem;
end

% Sort results by day
[sorted_days, sort_idx] = sort(sorted_days);
mean_new_pfc_avg_grouped = mean_new_pfc_avg_grouped(sort_idx);
sems_new_pfc_avg_grouped = sems_new_pfc_avg_grouped(sort_idx);

% pfc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of days to plot
num_days = numel(sorted_days);

figure;
tiledlayout(4,3);  

% Loop through each day and plot in a separate subplot
for i = 1:num_days

    if sorted_days(i) < 0
        day_field = sprintf('day_neg%d', abs(sorted_days(i)));  % For negative days, use day_neg3, day_neg2, etc.
    else
        day_field = sprintf('day_%d', sorted_days(i));  % For non-negative days, use day_3, day_2, etc.
    end

    % Data for the current day
    x = wf_stim_time;
    avg = mean_new_pfc_avg_grouped{i}';        % Average (1x151)
    sem = sems_new_pfc_avg_grouped{i}';            % SEM (1x151)
    
    % Get the size of the data for the current day (number of repetitions)
    num_reps = size(all_new_pfc_avg_grouped_combined.(day_field), 2);  % Number of rows (repetitions)

    % Create a new tile for each day
    nexttile;
    hold on;
    
    % Create the x and y values for the SEM region
    x_fill = [x, fliplr(x)];  % X values: time points and reversed time points
    y_fill = [avg + sem, fliplr(avg - sem)];  % Y values: upper and lower bounds of SEM
    
    % Plot the shaded SEM region
    fill(x_fill, y_fill, [0 0.5 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % Plot average as a line
    plot(x, avg, 'LineWidth', 1.5, 'Color', [0 0.5 0]);

    % Title and labels for each subplot
    title_str = sprintf('Day %d (Repetitions: %d)', sorted_days(i), num_reps);  % Add repetitions to the title
    title(title_str);
    xlabel('Time from stim onset');
    ylabel('\Delta F/F', 'FontSize', 12);
    xline(0, 'LineWidth', 2)
    ylim([-10^(-3) 2*10^(-3)])
    hold off;
    xlim([-0.1 0.3])
    AP_scalebar(0.1, 10^(-3))
    ap.prettyfig
end

%% SPARE
% figure; imagesc(vis_roi_mask); axis image;
% title('ROI mask');

vis_roi_trace = ap.wf_roi(wf_U,wf_V,wf_avg,[],vis_roi_mask);

figure('Position', [680 460 860 520]);
plot(wf_t, vis_roi_trace, 'color', [0 0.7 0]);
% title('Visual ROI fluorescence');
xline(stimOn_times(contra_good_trials), 'k', 'LineWidth', 2)
ylim([-0.02, 0.02])
xlim([60 120])
xlabel('Time (s)', 'FontSize', 34)
ylabel('{\Delta}F/F', 'FontSize', 34)
box off
AP_scalebar(10, 0.01)

%% - wf avg trials
roi_stim_align = interp1(wf_t,vis_roi_trace,time_stimulus)';
avg_contra_roi_stim_align = mean(roi_stim_align(:, contra_good_trials), 2);
figure('Position', [680 460 860 520]);
plot(timevec, avg_contra_roi_stim_align, 'color', [0 0.7 0])
xline(0, 'k', 'LineWidth', 2)
xline(0.5, 'k', 'LineWidth', 2)
xlabel('Time from stim onset (s)', 'FontSize', 34)
ylabel('{\Delta}F/F', 'FontSize', 34)
% title('Visual ROI fluorescence');
box off
AP_scalebar(0.2, 1*10^-3)
