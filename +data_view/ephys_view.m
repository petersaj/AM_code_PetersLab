%% Manual define info
animal = 'AP006';

% find passive protocol
protocol = 'lcr_passive*';
experiments = general.find_experiments(animal,protocol);
rec_day = experiments(end-3).day;
rec_time = experiments(end-3).experiment{1};

% load experiment
exp_info = general.AM_load_experiment(animal, rec_day, rec_time, bhv=true, ephys=true);

trialStim_values = vertcat(exp_info.trial_events.values.TrialStimX);

general.AM_AP_cellraster(exp_info, exp_info.stimOn_times,trialStim_values);


%% Andy's code
ap.load_experiment

trialStim_values = vertcat(trial_events.values.TrialStimX);
stimOn_times = photodiode_times(1:2:end);

AP_cellraster(stimOn_times,trialStim_values);

% define t
t =timelite.timestamps;

% all moves to left 
tmp_move = [0; diff(wheel_move)];
all_move_left_frames = find(tmp_move==1);
all_move_left_times = t(all_move_left_frames);
AP_cellraster(all_move_left_times);

% all moves to right
all_move_right_frames = find(tmp_move==-1);
all_move_right_times = t(all_move_right_frames);

AP_cellraster(all_move_right_times);

% % check times
% figure; 
% hold on; plot(all_move_left_times, 1, '-or');
% hold on; plot(all_move_right_times, 1, '-ob');
% figure; plot(t, wheel_move)
% hold on; plot(t, wheel_position)


%% move temp stuff
% define t
t =timelite.timestamps;

% all moves
tmp_move = [0; diff(wheel_move)];
all_move_on_frames = find(tmp_move==1);
all_move_on_times = t(all_move_on_frames);

general.AM_AP_cellraster(exp_info, all_move_on_times);
