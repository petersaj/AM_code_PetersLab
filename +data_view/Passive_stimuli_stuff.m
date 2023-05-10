%% Paths to stuff

addpath(genpath('C:\Users\Andrada\Documents\GitHub\npy-matlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\widefield'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\Lilrig'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\PupilDetection_DLC'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\AP_scripts_cortexlab'));

%% Initial stuff
mice = {'AP113', 'AP114', 'AP115'};

% retinotopy for each mouse
% mouse_idx = 1;
% 
% for mouse_idx=1:length(mice)
%     animal = mice{mouse_idx};
%     day = '2021-12-07';
%     experiment = 2;
%     verbose = true;
%     AP_load_experiment;
%     
%     lilrig_retinotopy
% end

%% Load experiment
% for day_idx=1:length(training_dates)
%     for mouse_idx=1:length(mice)
%         animal = mice{mouse_idx};
%         day = training_dates{day_idx};
%         experiment = 2;
%         if strcmp(day,'2021-12-13') && strcmp(animal,'AP107')
%             experiment = 3; 
%         end
% %         if day_idx<=3
% %             experiment = 1;
% %         end

for mouse_idx=1:length(mice)
    animal = mice{mouse_idx};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    for day_idx=length(experiments)-1:length(experiments) %1:length(experiments)
        day = experiments(day_idx).day;
        experiment = experiments(day_idx).experiment(end);
        verbose = true;
        AP_load_experiment;
        
        %% Processing stuff
        
        % alignment and deconvolution???
        % test_aligned = AP_align_widefield(Udf,animal,day);
        % deconvolved_fVdf = AP_deconv_wf(fVdf);
        
        %% Stimulus movies
        
        
        % find stimuli
        % trialStimulusValue = trial_conditions(:,2)./90;
        % possible_stimuli = unique(trialStimulusValue);
        
        
        % find possible stimuli
        possible_contrasts = unique(signals_events.stimContrastValues);
        possible_sides = unique(signals_events.stimAzimuthValues)/90;
        possible_stimuli = possible_sides.*possible_contrasts;
        possible_stimuli = unique(possible_stimuli);
        
        % find value of stimulus per trial
        trialStimulusValue = signals_events.stimAzimuthValues/90 .* signals_events.stimContrastValues;
        
        % create matrix of times for movie
        timestep = 0.01;
        start_time = -0.5;
        end_time = 1;
        timevec = [start_time:timestep:end_time];
        
        stim_frame = (-start_time)*(1/timestep)+1;
        
        time_stimulus = stimOn_times+timevec;
        
        % find activity for above time
        all_stim_act = interp1(frame_t,fVdf',time_stimulus);
        all_stim_act = permute(all_stim_act, [3,2,1]);
        
        
        % calculate average act for each stimulus
        all_stim_avg_act = nan(size(all_stim_act,1),size(all_stim_act,2),length(possible_stimuli));
        % completed_trialStimulusValue = trialStimulusValue(1:n_trials);
        
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        
        for stim_idx =1:length(possible_stimuli)
            this_stim_act = all_stim_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
            all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);
        end
        
        all_stim_avg_act = all_stim_avg_act - all_stim_avg_act(:,stim_frame,:);
        
        deconvolved_all_stim_avg_act = AP_deconv_wf(all_stim_avg_act, [], 1/timestep);
 
        % get fluoresence
        all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(Udf,deconvolved_all_stim_avg_act);
        
        % video
        AP_image_scroll(all_stim_interval_avg_fluorescence,timevec);
        axis image;
    end
end
