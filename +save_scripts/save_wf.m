%% SAVE
%% load bhv, save path and dataset
save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Users\Andrada-Maria_Marica\long_str_ctx_data';

% load('D:\matlab_save\swr_bhv.mat');

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

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
% load('D:\matlab_save\ROIs');

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
    bhv_days = {train_rec_passive.day}; % go through this instead
%     wf_days =  bhv_days(arrayfun(@(x) train_rec_passive(x).widefield(end), ...
%         1:length(train_rec_passive))); % add to ephys and other scripts

    wf_animal = table;

    for use_rec=1:length(train_rec_passive)
        rec_day = train_rec_passive(use_rec).day;
        rec_time = train_rec_passive(use_rec).recording{end};
        verbose = true;
        load_parts.behavior = true;
        load_parts.widefield = true;
        ap.load_recording

%         this_day_from_learning = days_from_learning(ismember(bhv_days, rec_day));

        %% -- get no move stuff
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

        %% -- align to master
        [U_master,V_master] = plab.wf.u2master(wf_U,wf_V);

        %% -- interp wf
%         V_interp = interp1(wf_t,wf_V',grab_time)';
        
        wf_stim_time = -0.5:0.01:1;
        wf_stim_interp_time = stimOn_times + wf_stim_time;

        V_stim_align = interp1(wf_t,V_master',wf_stim_interp_time);
        V_no_move_stim_align = V_stim_align(no_move_trials, :, :);

%         px_avg_contra_stim_align = plab.wf.svd2px(U_master,V_avg_contra_stim_align);
%         ap.imscroll(px_avg_contra_stim_align);
%         axis image
%         colormap(ap.colormap('PWG'));
%         title('Example map');
%         clim_val = max(abs(clim)) * 0.7;
%         clim([-clim_val, clim_val]); % Set the color limits
% %         clim(clim/8)


        %% -- save
        wf_animal.animal(use_rec) = {animal};
        wf_animal.rec_day(use_rec) = {rec_day};

        wf_animal.wf_stim_time(use_rec) = {wf_stim_time};
        wf_animal.V_no_move_stim_align(use_rec) = {V_no_move_stim_align};
        wf_animal.trial_stim_values(use_rec) = {trial_stim_values(no_move_trials)};

        disp(['Done day ' num2str(use_rec)])
        
    end
    all_wf_save_cell{animal_idx} = wf_animal;
    disp(['Done ' animal])
end

wf = vertcat(all_wf_save_cell{:});

% ctx_wf.U_master(:) = {U_master};

% ctx_wf.vis_roi_mask(:) = {vis_roi_mask};
% ctx_wf.pfc_roi_mask(:) = {pfc_roi_mask};

save_name = fullfile(save_path, 'ctx_wf');
save(save_name, "wf", "-v7.3");

% copy master
copy_master_path = fullfile(save_path, 'U_master.mat');
copyfile(master_U_fn, copy_master_path);
