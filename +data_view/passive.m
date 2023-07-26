%% Initial stuff

animal = 'AP004';
rec_day = '2023-05-11';
rec_time = '1345';
exp_info = general.AM_load_experiment(animal, rec_day, rec_time,"bhv",true, "wf", true);


%% Stimulus movies

% possible stim
possible_stimuli = [-90 0 90];

% find value of stimulus per trial
trialStimulusValue = [exp_info.trial_events.values.TrialStimX];
trialStimulusValue = trialStimulusValue(:);

stimOn_times = exp_info.photodiode_times(1:2:end);
% check where stim on times are
% figure; plot(exp_info.timelite.timestamps, exp_info.timelite.data(:, 5));
% hold on; plot(stimOn_times, 0, 'o');

% create matrix of times for movie
timestep = 0.01;
start_time = -0.5;
end_time = 1;
timevec = [start_time:timestep:end_time];

stim_frame = (-start_time)*(1/timestep)+1;

time_stimulus = stimOn_times+timevec;

% find activity for above time
all_stim_act = interp1(exp_info.wf_times,exp_info.wf_V',time_stimulus);
all_stim_act = permute(all_stim_act, [3,2,1]);

% calculate average act for each stimulus
all_stim_avg_act = nan(size(all_stim_act,1),size(all_stim_act,2),length(possible_stimuli));
% completed_trialStimulusValue = trialStimulusValue(1:n_trials);

% stim aligned wheel move
t = exp_info.timelite.timestamps;
stim_wheel_move = interp1(t,+exp_info.wheel_move,time_stimulus);
no_move_trials = sum(stim_wheel_move(:,stim_frame:end),2)==0;

for stim_idx =1:length(possible_stimuli)
    this_stim_act = all_stim_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
    all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);
end

all_stim_avg_act = all_stim_avg_act - all_stim_avg_act(:,stim_frame,:);

% deconvolved_all_stim_avg_act = AP_deconv_wf(all_stim_avg_act, [], 1/timestep);

% get fluoresence
% all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(Udf,deconvolved_all_stim_avg_act);
all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(exp_info.wf_U, all_stim_avg_act);

% video
AP_imscroll(all_stim_interval_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); colorbar; caxis([-0.01 0.01]); % caxis([10^(-3) 10*10^(-3)])
axis image;

% check wheel trace
% figure;
% plot(exp_info.timelite.timestamps, exp_info.wheel_move);
% hold on;
% plot(stimOn_times, 0, 'o')


